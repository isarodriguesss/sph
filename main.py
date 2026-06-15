import csv
import numpy as np
from scipy.spatial import cKDTree
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme

LOG_FILE = "log.csv"
LOG_HEADER = [
    "t",
    "iteration",
    "max_v",
    "mean_v",
    "n_fast",
    "a_marangoni",
    "a_drag",
    "a_pressure",
    "a_flag",
    "a_total",
    "min_cs",
    "max_cs",
    "mean_cs",
    "constrast_cs",
    # Pass M-B — campo de nutriente consumivel c_n (Frente 6)
    "min_c_n",
    "max_c_n",
    "mean_c_n",
    "contrast_c_n",
    "mass_total",  # soma de m — cresce por BiomassGrowth (logistico)
    "pass_n_spawned",  # partículas criadas pelo Pass N desde a última linha de log
    # Instrumentacao do #2 — particao da unidade sigma_a na FRONTIER ATIVA
    # (braços/dendritos, rho_b ∈ [0.1, 0.5]). Mede se o vacuo da frontier e
    # deficit REAL de operador SPH (sigma_a < 0.85, Violeau §3.6 / Liu §3.3.3)
    # ou se a KGC ja o neutralizou. Se mean_sig_arms ≈ 1 e frac_lowsig_arms
    # pequeno → vacuo visual e cosmetico (agar entre dendritos). Se sigma_a
    # baixo e fracao alta → deficit real, decide se C4/PassN/finer-res e preciso.
    "min_sig_arms",  # min sigma_a nos braços (pior deficit de kernel)
    "mean_sig_arms",  # media sigma_a nos braços
    "frac_lowsig_arms",  # fracao de particulas-braço com sigma_a < 0.85
]

x_dim, y_dim = 187, 187  # M-B.10: expandido para preservar dx≈0.054 em dominio 10x10

# Domínio 10×10 centrado na origem (M-B.10: era 8x8 [-4,4]; ampliado para
# evitar colisao de dendritos com paredes invisiveis em t>50s — diagnosticado
# em M-B.9b como mecanismo dominante de fragmentacao pos-50s)
x_min_domain, x_max_domain = -5.0, 5.0
y_min_domain, y_max_domain = -5.0, 5.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)


# ── Parâmetros Calibrados (Coesão Viscosa + Interface Gateada) ──────────
# Meta: v_term ≈ 0.1, acelerações totais 10–50, Marangoni só na interface.
# Coesão vem de: (a) EOS com ramo atrativo (p<0 se rho<rho0),
#                (b) viscosidade maior, (c) Monaghan artificial viscosity forte,
#                (d) kernel com ~35 vizinhos (h_factor=1.8).
mu = 0.020  # K.23: I.2 revert parcial — fortalecer coesão viscosa (I.2 era 0.012, pre-I.2 era 0.025)
gamma = 60.0  # Drag: a_drag = gamma * v_term = 60 * 0.1 = 6
beta = 5.0  # Pass T2d: T1 β=10 → 5 (reduz tracao amplificada na fronteira); razao |F_mar|/B_tension ~4800×
sigma = 10.0  # Pass T2g: T2d σ=20 → 10 — desacelera saturacao cs; reduz mean_cs e tracao no rim (combate colapso de contrast_cs e fragmentacao pos-t25s diagnosticada em Pass N v2.4)
D = 1.5e-3  # Pass I.3: D_int dentro do biofilme — gradiente afiado na interface
D_ext = 0.08  # K.18: 4x — L_D_ext=0.298 (~2x maior); habilita focalizacao Mullins-Sekerka pos-K.17
lambda_ = 0.15  # Decaimento: confina cs mas permite penetracao de ~L_D_ext no exterior
r_growth = 0.02  # M-B.7: 0.05→0.02 — reduz tip pumping (mass cresceu 9.7× em M-B.6 via growth nas pontas com c_n>0.8 a taxa plena)
rho_max = 1.0
alpha_mon = 0.12  # K.23: I.2 revert parcial — previne instabilidade de tração SPH (I.2 era 0.06, pre-I.2 era 0.15)

# Pass M-A (Frente 6): osmolito c_o secretado pelas bacterias.
# Mecanismo: c_o satura nas baias (agar confinado entre dendritos), fica ~0
# nas pontas (agar virgem). |∇c_o| dispara influxo de massa via van't Hoff.
D_o = 0.04  # Pass M-A.2: L_D_o = sqrt(0.04/0.05) = 0.894 ~ d_tip-tip (Mullins-Sekerka)
k_o = 0.5  # taxa de producao por bacteria
lambda_o = 0.05  # decaimento lento (osmolitos persistem mais que cs)
Q0 = 5.0  # Pass M-A.2: forca osmotica Darcy — a_osm_tip ~ Q0*gate*|grad_c_o| ~ 2.8

# Pass M-B (Frente 6): nutriente consumivel c_n.
# c_n inicia em 1.0 (agar virgem). Bacterias consomem na taxa k_n*rho_b.
# Difusao tem que dominar consumo (tau_cons/tau_dif > 4) para evitar morte
# quimica global. Lição da rodada inicial M-B (k_n=1.0, D_n=1e-3): consumo
# dominava → motor cs morria em t<2s. Calibracao corrigida: razao = 25.
D_n = 0.02  # difusao do nutriente no agar (D_n_ext — livre)
D_n_int = 1e-4  # M-B.4: difusao dentro do biofilme (EPS bloqueia transporte)
k_n = 0.5  # taxa de consumo por unidade de biomassa

dt_global = 0.001
total_sim_time = 100.0  # C3.2 validada em t=50s (vacuo profundo preenchido, massa +2.1%, inserções auto-limitantes); estender p/ confirmar sustentabilidade
print_freq = 200

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35  # EOS: B = 1.5²/7 ≈ 0.32 (repulsão suave, atração ~0.1)

use_splitting = False

# ═══ Tratamento do vácuo — OPÇÃO 1 (baseline validado, 2026-06-15) ═══════
# KGC corrige o OPERADOR ∇cs na frontier (pontas σ_a~0.98, medição #2) +
# inserção C3.4 frozen preenche o vácuo ESTRUTURAL (núcleo). Pass N e Shifting
# ABANDONADOS (binds #31/#36/#38). Histórico completo: CLAUDE.md §12.

# Rota A — Particle Shifting (Xu 2009; Lind 2012). ESGOTADA (#31/#32). Preservada.
use_shift = False
SHIFT_COEFF = 0.5
SHIFT_CAP = 0.05
SHIFT_RHO_B_MIN = 0.6

# Rota C — Kernel Gradient Correction (Bonet-Lok 1999; CSPM, Liu §3.3).
# Corrige ∇cs (consistência de 1ª ordem) onde σ_a<1 SEM mover/criar partícula.
# Auto-gateia por det(M): bulk det≈1 (sem correção), frontier det<1 (corrige).
use_kgc = True  # OPÇÃO 1
KGC_DET_MIN = 0.25  # fallback identidade se det(M)<0.25 (|L|≲4×, anti-spike)

# Rota C3.4 — Inserção de filler INERTE no vácuo estrutural (rede dx, h=h0 →
# dt intacto, sem over-pack; Liu §6.5). Frozen-só: filhas não crescem/produzem
# cs e são pinadas (evita runaway #34). Frontier (rho_b<0.5) fica com a KGC.
use_insert = True  # OPÇÃO 1
INSERT_FREQ = 200  # iter entre inserções
INSERT_SIGMA_TRIG = 0.85  # déficit de partição da unidade (Violeau §3.6)
INSERT_RHO_B_MIN = 0.1  # só vácuo estrutural; para frontier: 0.1
INSERT_PROX = 0.7  # insere só em buraco real (>0.7·dx de qualquer vizinho)
INSERT_MAX = 100  # cap de inserções/call (controle de massa)

# Pass N — refinamento hexagonal Vacondio/Feldman (Liu §3.3.3, Violeau §7.4).
# ABANDONADO: bind over-pack↔dt sob dt global (lição #38). Preservado p/ Fase 2
# (timesteps individuais [T9]). Constantes mantidas para reativação.
use_pass_n = False
PASS_N_FREQ = 100
PASS_N_MAX_PARENTS = 100
PASS_N_SIGMA_TRIG = 0.95  # preemptivo (canônico Violeau §3.6 = 0.85)
PASS_N_ALPHA = 0.75  # h_filha=α·h_mãe (grande → preserva dt; #29)
PASS_N_EPSILON = 0.5  # offset=ε·h_mãe (ε/α=0.67; ótimo Feldman=1)
PASS_N_MAX_GEN = 1  # contador de geração explícito (imune ao growth; #37)
PASS_N_RHO_B_MIN = 0.15
PASS_N_RHO_B_MAX = 0.95
PASS_N_PROXIMITY_MIN = (
    0.4  # min dist. filha-vizinho (dx); garante n_d=7 (Violeau §7.4.3)
)


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = (
            None  # snapshot da massa total em t=0 (preenchido em post_step)
        )
        self._pass_n_spawned_since_log = 0  # acumula spawns entre linhas de log

    def create_particles(self):
        fluid_solid = create_initial_state(
            x_dim,
            y_dim,
            rho_max,
            dt_global,
            x_min=x_min_domain,
            x_max=x_max_domain,
            y_min=y_min_domain,
            y_max=y_max_domain,
        )

        for pa in fluid_solid:
            if pa.name == "fluid":
                pa.add_property("noise")
                pa.noise[:] = (
                    1.0
                    + 0.6 * np.sin(8 * np.arctan2(pa.y, pa.x))
                    + 0.01 * np.random.rand(len(pa.x))
                )
                pa.add_property("dt_force")
                pa.add_property("dt_cfl")
                # debug das acelerações (componentes vetoriais + magnitude)
                pa.add_property("au_mar")  # |aceleração Marangoni| (líquida)
                pa.add_property("ax_mar")  # aceleração Marangoni componente x
                pa.add_property("ay_mar")  # aceleração Marangoni componente y
                pa.add_property("au_drag")  # |aceleração Drag|
                # gradiente da densidade
                pa.add_property("grad_rho_b_x")
                pa.add_property("grad_rho_b_y")
                pa.add_property("grad_rho_b_mag")
                # gradiente do surfactante (usado por FlagellarForce)
                pa.add_property("grad_cs_x")
                pa.add_property("grad_cs_y")
                # gradiente do osmolito c_o (Pass M-A — usado por OsmolyteProduction)
                pa.add_property("grad_co_x")
                pa.add_property("grad_co_y")
                # Pass N v2.4 — partição da unidade discreta sigma_a = Σ V_j W_aj
                # (Violeau §3.4-3.6, Liu §3.3.3). Trigger de refinamento.
                pa.add_property("sigma_a")
                pa.sigma_a[:] = 1.0  # inicializa em 1 (consistencia perfeita)
                # Rota A — Fickian shifting: gradiente de concentracao (dC) e
                # vetor de deslocamento (shift). Aplicado a posicao no integrador.
                pa.add_property("shift_dC_x")
                pa.add_property("shift_dC_y")
                pa.add_property("shift_x")
                pa.add_property("shift_y")
                pa.shift_x[:] = 0.0
                pa.shift_y[:] = 0.0
                # Rota C — KGC: matriz de renormalizacao M (acumulada) e sua
                # inversa L (correcao do kernel gradient). L inicia na IDENTIDADE
                # → se KGC desligado, MarangoniForce reduz ao SPH padrao.
                pa.add_property("Mxx")
                pa.add_property("Mxy")
                pa.add_property("Myx")
                pa.add_property("Myy")
                pa.add_property("Lxx")
                pa.add_property("Lxy")
                pa.add_property("Lyx")
                pa.add_property("Lyy")
                pa.Lxx[:] = 1.0
                pa.Lxy[:] = 0.0
                pa.Lyx[:] = 0.0
                pa.Lyy[:] = 1.0
                # C3.3 — marcador de filler inerte (Rota C3): 0 = partícula real
                # (cresce, produz cs, dinâmica normal); 1 = inserida p/ suporte de
                # kernel (não cresce, não produz cs, pinada). Quebra o feedback de
                # nucleus maturation + inflação de cs do runaway t>57s.
                pa.add_property("is_filler")
                pa.is_filler[:] = 0.0
                # Pass N — contador de GERAÇÃO de refinamento. gen=0 = partícula
                # original; cada split incrementa. Substitui o limite por massa
                # (m > m_floor), que era DERROTADO pelo BiomassGrowth: a massa
                # cresce de volta acima do floor → filha re-splita → gen-2,3...
                # → h encolhe → dt colapsa (diagnóstico 2026-06-15, gen-2 com
                # h_min=0.56·h0 + over-pack rho=2.94 → dt 97× pior). gen nunca
                # decresce → imune ao crescimento.
                pa.add_property("gen")
                pa.gen[:] = 0.0
                pa.add_property("ax_drag")
                pa.add_property("ay_drag")
                # flag
                pa.add_property("au_flag")
                # Pass N — posição de referência para rastreamento de deslocamento
                pa.add_property("x_spawn_ref")
                pa.add_property("y_spawn_ref")
                pa.x_spawn_ref[:] = pa.x[:]
                pa.y_spawn_ref[:] = pa.y[:]
                pa.add_output_arrays(
                    [
                        "rho_b_grown",
                        "cs",
                        "c_o",
                        "u",
                        "v",
                        "p",
                        "noise",
                        "au_flag",
                        "au_mar",
                        "sigma_a",  # #2 — particao da unidade p/ 4º painel do plot
                    ]
                )
            elif pa.name == "solid":
                pa.add_property("p")

        return fluid_solid

    def create_scheme(self):
        return MyBiomassScheme(
            fluids=["fluid"],
            solids=["solid"],
            dim=2,
            mu=mu,
            gamma=gamma,
            beta=beta,
            sigma=sigma,
            D=D,
            D_ext=D_ext,
            lambda_=lambda_,
            r_growth=r_growth,
            rho_max=rho_max,
            c0=c0,
            alpha_mon=alpha_mon,
            D_o=D_o,
            k_o=k_o,
            lambda_o=lambda_o,
            Q0=Q0,
            D_n=D_n,
            D_n_int=D_n_int,
            k_n=k_n,
            use_shift=use_shift,
            shift_coeff=SHIFT_COEFF,
            shift_cap=SHIFT_CAP,
            shift_rho_b_min=SHIFT_RHO_B_MIN,
            use_kgc=use_kgc,
            kgc_det_min=KGC_DET_MIN,
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        solver = Solver(
            dim=2,
            integrator=self.scheme.get_integrator(),
            kernel=kernel,
            dt=5e-5,
            adaptive_timestep=True,
            cfl=0.4,
        )
        solver.tf = total_sim_time
        solver.set_print_freq(print_freq)
        return solver

    def post_step(self, solver):
        # Print stats
        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            v_mag = np.sqrt(fluid.u**2 + fluid.v**2)
            max_v = np.max(v_mag)
            mean_v = np.mean(v_mag)
            n_fast = int(np.sum(v_mag > 0.1))

            # Marangoni líquido (magnitude do vetor, não soma de normas)
            a_mar = np.max(np.abs(fluid.au_mar))
            a_drag = np.max(np.abs(fluid.au_drag))
            # Aceleração total (inclui pressão via MomentumEquation + tudo)
            a_total = np.max(np.sqrt(fluid.au**2 + fluid.av**2))
            # ax_pressure = au_total - ax_mar - ax_drag
            ax_p = fluid.au - fluid.ax_mar - fluid.ax_drag
            ay_p = fluid.av - fluid.ay_mar - fluid.ay_drag
            a_pressure = np.max(np.sqrt(ax_p**2 + ay_p**2))
            a_flag = np.max(np.abs(fluid.au_flag))

            # 3. Estatísticas do Surfactante (cs)
            min_cs = np.min(fluid.cs)
            max_cs = np.max(fluid.cs)
            mean_cs = np.mean(fluid.cs)
            contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

            # 4. Pass M-B — Estatísticas do nutriente c_n
            # c_n inicia em 1.0 (agar virgem). Bacterias consomem na taxa k_n*rho_b.
            # mean_c_n esperado: cair de 1.0 mas estabilizar em ~0.3-0.5 (difusao
            # repoe nutriente do exterior). Se cai para <0.1, motor cs vai morrer.
            # contrast_c_n alto = baias esgotadas vs pontas em agar virgem (selecao).
            min_c_n = np.min(fluid.c_n)
            max_c_n = np.max(fluid.c_n)
            mean_c_n = np.mean(fluid.c_n)
            contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

            # 5. Massa total
            # Cresce por BiomassGrowth (logistico). Sem bug osmotico desde M-A.2.
            mass_total = float(np.sum(fluid.m))

            # 6. #2 — particao da unidade sigma_a na FRONTIER ATIVA (braços).
            # Zona motil rho_b ∈ [0.1, 0.5] (lição #35): exclui agar (<0.1),
            # nucleo/junção estrutural (>0.5, ja tratada por C3 filler). Mede
            # se o vacuo dos dendritos e deficit real de kernel (Violeau §3.6).
            arms_mask = (fluid.rho_b_grown >= 0.1) & (fluid.rho_b_grown < 0.5)
            n_arms = int(np.sum(arms_mask))
            if n_arms > 0:
                sig_arms = fluid.sigma_a[arms_mask]
                min_sig_arms = float(np.min(sig_arms))
                mean_sig_arms = float(np.mean(sig_arms))
                frac_lowsig_arms = float(np.mean(sig_arms < 0.85))
            else:
                min_sig_arms = mean_sig_arms = 1.0
                frac_lowsig_arms = 0.0

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
            )
            print(
                f"sigma_a braços (rho_b∈[0.1,0.5]): min={min_sig_arms:.3f} "
                f"mean={mean_sig_arms:.3f} frac<0.85={frac_lowsig_arms:.2%}"
            )
            print("Acelerações:")
            print(f"  > Marangoni (líq): {a_mar:.2f}")
            print(f"  > Drag:            {a_drag:.2f} (Freio)")
            print(f"  > Pressão (est):   {a_pressure:.2f}")
            print(f"  > Total |a|:       {a_total:.2f}")

            with open(LOG_FILE, "a", newline="") as f:
                csv.writer(f).writerow(
                    [
                        f"{solver.t:.4f}",
                        solver.count,
                        f"{max_v:.6f}",
                        f"{mean_v:.6f}",
                        n_fast,
                        f"{a_mar:.4f}",
                        f"{a_drag:.4f}",
                        f"{a_pressure:.4f}",
                        f"{a_flag:.4f}",
                        f"{a_total:.4f}",
                        f"{min_cs:.4f}",
                        f"{max_cs:.4f}",
                        f"{mean_cs:.4f}",
                        f"{contrast_cs:.4f}",
                        f"{min_c_n:.4f}",
                        f"{max_c_n:.4f}",
                        f"{mean_c_n:.4f}",
                        f"{contrast_c_n:.4f}",
                        f"{mass_total:.6e}",
                        self._pass_n_spawned_since_log,
                        f"{min_sig_arms:.4f}",
                        f"{mean_sig_arms:.4f}",
                        f"{frac_lowsig_arms:.4f}",
                    ]
                )
                self._pass_n_spawned_since_log = 0  # reseta após registrar

        if use_splitting:
            print("Iniciando processo de divisão celular...")
            if solver.count > 0 and solver.count % 500 == 0:
                fluid = self.particles[0]
                min_rho_b_for_division = 0.05
                mature_by_mass_indices = np.where(fluid.m > 1.99 * fluid.m0)[0]
                mature_by_rho_b_indices = np.where(
                    fluid.rho_b_grown[mature_by_mass_indices] > min_rho_b_for_division
                )[0]
                indices_to_split = []
                for idx in mature_by_rho_b_indices:
                    original_idx = mature_by_mass_indices[idx]
                    if np.random.rand() < prob_of_splitting:
                        indices_to_split.append(original_idx)

                if len(indices_to_split) > 0:
                    print(
                        f"\n--- Divisão Celular em t={solver.t:.2f}: {len(indices_to_split)} partículas se dividindo. ---"
                    )

                    daughters = fluid.empty_clone()

                    props_to_copy = [
                        "x",
                        "y",
                        "m",
                        "h",
                        "rho",
                        "rho_b_grown",
                        "cs",
                        "u",
                        "v",
                        "au",
                        "av",
                        "a_rho_b_grown",
                        "a_c_s",
                        "m0",
                    ]

                    for parent_idx in indices_to_split:
                        parent_props = {
                            prop: getattr(fluid, prop)[parent_idx]
                            for prop in props_to_copy
                        }

                        for i in range(2):
                            daughter_data = parent_props.copy()
                            daughter_data["m"] = parent_props["m"] / 2.0
                            daughter_data["m0"] = parent_props["m0"]
                            # daughter_data["rho_b_grown"] = parent_props["rho_b_grown"] / 2.0
                            dx_local = parent_props["h"] * 0.5
                            offset_x = dx_local * (np.random.rand() - 0.5)
                            offset_y = dx_local * (np.random.rand() - 0.5)

                            daughter_data = parent_props.copy()
                            daughter_data["x"] = parent_props["x"] + offset_x
                            daughter_data["y"] = parent_props["y"] + offset_y

                            data_to_add = {
                                key: [value] for key, value in daughter_data.items()
                            }

                            daughters.add_particles(**data_to_add)

                    fluid.append_parray(daughters)

                    fluid.remove_particles(indices_to_split)

                    solver.nnps.update()

        # Pass N - Refinamento hexagonal Vacondio/Feldman com consistencia
        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            # Gate v2.3.1: rho_b ∈ [0.3, 0.7] — exclui borda dilute (rho_b<0.3,
            # gap não real, só kernel truncado) E núcleo/shell adjacente ao pin
            # (rho_b>0.7, onde refinamento e irrelevante porque momentum nao e
            # integrado — hard pin K.17).
            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            # Limite de geração via CONTADOR EXPLÍCITO (gen < PASS_N_MAX_GEN),
            # NÃO por massa. O gate antigo (m > m₀/7) assumia que a massa só
            # diminui ao dividir, mas o BiomassGrowth cresce m de volta acima do
            # floor → filha re-splitava → gen-2 espúrio (h_min=0.56·h0) que, junto
            # ao over-pack, colapsou o dt 97% na Fase 1 (diag. 2026-06-15). gen
            # nunca decresce → imune ao crescimento.
            splittable_gen_mask = fluid.gen < PASS_N_MAX_GEN

            # Pass N v2.4 — Trigger por PARTIÇÃO DA UNIDADE (sigma_a < 0.85).
            # sigma_a = Sum_j V_j W_aj e a consistencia de ordem zero do SPH
            # (Violeau §3.6, Liu §3.3.3). Mede diretamente o erro do operador
            # SPH no nó a. INVARIANTE SOB REFINAMENTO — gen 0/1 disparam
            # pelo mesmo limiar. v2.2-v2.3.1 usavam rho_rel < 0.7 que dependia
            # da massa (rho = Sum m_j W) → trigger ficava enviesado em filhas
            # com m_d = m_m/7.
            sigma_a = fluid.sigma_a
            split_mask = (
                colony_mask & (sigma_a < PASS_N_SIGMA_TRIG) & splittable_gen_mask
            )
            split_idx = np.where(split_mask)[0]

            if len(split_idx) > 0:
                # Prioriza partição da unidade mais degradada (menor sigma_a)
                # se exceder limite — refina onde o operador SPH e pior.
                if len(split_idx) > PASS_N_MAX_PARENTS:
                    order = np.argsort(sigma_a[split_idx])
                    split_idx = split_idx[order][:PASS_N_MAX_PARENTS]

                eps = PASS_N_EPSILON
                alpha = PASS_N_ALPHA
                prox_min = PASS_N_PROXIMITY_MIN * dx
                prox_min_sq = prox_min * prox_min

                # Padrão hexagonal: 6 vértices a 60°
                angles_hex = np.arange(6) * (np.pi / 3.0)
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                # KDTree de partículas existentes — proximity guard.
                # Pass N v2.3: excluir as mães do split atual da árvore. Elas serão
                # removidas em remove_particles ao final desta call — manter elas
                # na árvore fazia com que o filtro rejeitasse vértices que caem a
                # r_off < prox_min da própria mãe. Para gen-1 (h=1.08·dx_0) o
                # r_off=0.378·dx_0 é < prox_min=0.4·dx_0, então todos os 6 vértices
                # falhavam → n_d=1 < 4 → split abortado silenciosamente. Resultado
                # v2.2: filhas gen-1 nunca conseguiam se redividir; núcleo
                # esvaziava progressivamente sem reposição. Ver §9 Pass N v2.3.
                tree_mask = np.ones(len(fluid.x), dtype=bool)
                tree_mask[split_idx] = False
                positions = np.column_stack([fluid.x[tree_mask], fluid.y[tree_mask]])
                tree = cKDTree(positions)

                daughters = fluid.empty_clone()
                new_positions = []  # filhas já adicionadas nesta call
                mothers_used = []  # índices de mães que efetivamente splitaram
                n_daughters_total = 0

                # Snapshot dos campos (índices estáveis até o remove)
                x_arr = fluid.x[split_idx]
                y_arr = fluid.y[split_idx]
                h_arr = fluid.h[split_idx]
                m_arr = fluid.m[split_idx]
                u_arr = fluid.u[split_idx]
                v_arr = fluid.v[split_idx]
                rho_arr = fluid.rho[split_idx]
                rhob_arr = fluid.rho_b_grown[split_idx]
                cs_arr = fluid.cs[split_idx]
                co_arr = fluid.c_o[split_idx]
                cn_arr = fluid.c_n[split_idx]
                noise_arr = fluid.noise[split_idx]
                sigma_arr = fluid.sigma_a[split_idx]
                gen_arr = fluid.gen[split_idx]

                for k in range(len(split_idx)):
                    x_m, y_m = float(x_arr[k]), float(y_arr[k])
                    h_m = float(h_arr[k])
                    m_m = float(m_arr[k])
                    r_off = eps * h_m

                    # Filtra vértices que colidem com vizinhos existentes
                    # OU com filhas já adicionadas nesta call.
                    valid_vertices = []
                    for j in range(6):
                        vx = x_m + r_off * cos_a[j]
                        vy = y_m + r_off * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox_min:
                            continue
                        too_close = False
                        for nx, ny in new_positions:
                            if (vx - nx) ** 2 + (vy - ny) ** 2 < prox_min_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        valid_vertices.append((vx, vy))

                    # Pass N v2.4 — EXIGENCIA ESTRITA n_d = 7 (Violeau §7.4.3):
                    # apenas distribuicao hexagonal SIMETRICA preserva (i) erro
                    # de densidade < 5% pos-split (Feldman 2006), (ii) centro de
                    # massa na posicao da mae, (iii) tensor de inercia local.
                    # Splits parciais quebram momento angular e introduzem torque
                    # espurio (Liu §4.2.4) que amplifica instabilidade de tracao.
                    # Se algum vertice falha no proximity guard, split e ADIADO
                    # para proxima call (gap precisa esticar mais).
                    n_d = 1 + len(valid_vertices)
                    if n_d != 7:
                        continue

                    # Massa REDISTRIBUÍDA entre 7 filhas (Liu §3.4 — conservacao
                    # de momentum linear: Sum m_d v_d = 7·(m_m/7)·v_m = m_m·v_m)
                    m_d = m_m / 7.0
                    h_d = alpha * h_m  # v2.4: alpha=0.35 (acoplado a eps=0.35)
                    u_d = float(u_arr[k])
                    v_d = float(v_arr[k])

                    xs = [x_m] + [p[0] for p in valid_vertices]
                    ys = [y_m] + [p[1] for p in valid_vertices]

                    # Atualiza tracking de posições para próximas mães
                    new_positions.append((x_m, y_m))
                    new_positions.extend(valid_vertices)

                    data = {
                        "x": xs,
                        "y": ys,
                        "m": [m_d] * 7,
                        "h": [h_d] * 7,
                        "rho": [float(rho_arr[k])] * 7,
                        "rho_b_grown": [float(rhob_arr[k])] * 7,
                        "cs": [float(cs_arr[k])] * 7,
                        "c_o": [float(co_arr[k])] * 7,
                        "c_n": [float(cn_arr[k])] * 7,
                        "u": [u_d] * 7,
                        "v": [v_d] * 7,
                        "noise": [float(noise_arr[k])] * 7,
                        "sigma_a": [float(sigma_arr[k])] * 7,
                        # is_filler=0: filhas sao BIOMASSA REAL refinada (crescem,
                        # produzem cs, movem). NAO setar isto deixa add_particles
                        # herdar lixo do realloc → filha podia nascer is_filler>0.5
                        # e ser PINADA como filler congelado (scheme.py), matando o
                        # refinamento. Mesmo risco do bug mean_c_n>max_c_n da v2.3.
                        "is_filler": [0.0] * 7,
                        # gen = mãe+1 — contador explícito de geração (imune ao
                        # BiomassGrowth). Garante gen≤MAX_GEN sem depender da massa.
                        "gen": [float(gen_arr[k]) + 1.0] * 7,
                    }
                    daughters.add_particles(**data)
                    mothers_used.append(int(split_idx[k]))
                    n_daughters_total += 7

                if mothers_used:
                    fluid.append_parray(daughters)
                    fluid.remove_particles(np.asarray(mothers_used, dtype=np.uint32))
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_daughters_total
                    print(
                        f"Pass N v2.4 t={solver.t:.1f}s: "
                        f"{len(mothers_used)} mães → {n_daughters_total} filhas "
                        f"(n_d=7 estrito, ε={eps}, α={alpha}, σ_trig={PASS_N_SIGMA_TRIG}, gen≤1)"
                    )

        # Rota C3 — Inserção de partículas no vácuo (em espaçamento dx)
        # Preenche o vácuo central INSERINDO partículas frescas na rede dx
        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx  # massa alvo = volume da rede inicial

            # Vácuo ESTRUTURAL (rho_b>INSERT_RHO_B_MIN=0.5, exclui agar e frontier
            # motil) com deficit de suporte de kernel σ_a < 0.85 (Violeau §3.6 —
            # critério canônico, KernelSum a cada passo). C3.4 frozen-só: filler
            # inerte/congelado no nucleo/junção. C4 (estender p/ frontier ativa +
            # filler ativo) FALHOU — runaway #34, revertido (lição #36).
            void_mask = (fluid.rho_b_grown > INSERT_RHO_B_MIN) & (
                fluid.sigma_a < INSERT_SIGMA_TRIG
            )
            void_idx = np.where(void_mask)[0]

            if len(void_idx) > 0:
                positions = np.column_stack([fluid.x, fluid.y])
                tree = cKDTree(positions)
                prox = INSERT_PROX * dx
                prox_sq = prox * prox

                # Candidatos: 6 vizinhos hexagonais a distância dx de cada
                # partícula-vácuo. Insere onde o spot está VAZIO (rede com buraco).
                angles_hex = np.arange(6) * (np.pi / 3.0)
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                new_x = []
                new_y = []
                parent_idx = []
                added_pts = []  # dedupe entre candidatos da mesma call

                for k in void_idx:
                    xk = float(fluid.x[k])
                    yk = float(fluid.y[k])
                    for j in range(6):
                        vx = xk + dx * cos_a[j]
                        vy = yk + dx * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox:
                            continue  # spot ocupado — não é buraco
                        too_close = False
                        for ax_, ay_ in added_pts:
                            if (vx - ax_) ** 2 + (vy - ay_) ** 2 < prox_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        new_x.append(vx)
                        new_y.append(vy)
                        parent_idx.append(int(k))
                        added_pts.append((vx, vy))
                    if len(new_x) >= INSERT_MAX:
                        break

                if len(new_x) > 0:
                    n_ins = len(new_x)
                    p = np.asarray(parent_idx, dtype=int)
                    inserted = fluid.empty_clone()
                    # Campos herdados da mãe-vácuo (interpolação local simples).
                    # rho recomputado por SummationDensity no proximo passo; L,
                    # sigma_a recomputados por KGC/KernelSum — nao precisam init.
                    data = {
                        "x": new_x,
                        "y": new_y,
                        "m": [m_target] * n_ins,
                        "h": list(fluid.h[p]),
                        "rho": list(fluid.rho[p]),
                        "rho_b_grown": list(fluid.rho_b_grown[p]),
                        "cs": list(fluid.cs[p]),
                        "c_o": list(fluid.c_o[p]),
                        "c_n": list(fluid.c_n[p]),
                        "u": [0.0] * n_ins,
                        "v": [0.0] * n_ins,
                        "noise": list(fluid.noise[p]),
                        "is_filler": [1.0]
                        * n_ins,  # C3.3 — inerte (não cresce/produz, pinada)
                    }
                    inserted.add_particles(**data)
                    fluid.append_parray(inserted)
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_ins
                    print(
                        f"C3 inserção t={solver.t:.1f}s: {n_ins} partículas no "
                        f"vácuo (rede dx, h=h0 → dt intacto)"
                    )


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
