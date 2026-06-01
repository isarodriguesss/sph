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
total_sim_time = 50.0  # Pass N-on validacao curta: medir custo dt real (gen-1) antes de comprometer 100s
print_freq = 200

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35  # EOS: B = 1.5²/7 ≈ 0.32 (repulsão suave, atração ~0.1)

use_splitting = False

# Pass N v2.4 — Refinamento hexagonal Vacondio 2013 / Feldman 2006 com
# consistencia SPH restaurada (Liu §3.3.3, Violeau §3.4-3.6, §7.4-7.5).
# Substitui 1 mãe por 7 filhas (1 centro + 6 vértices a 60°). Massa dividida
# IGUALMENTE (m_filha = m_mãe/7); velocidade IDÊNTICA herdada (Liu §3.4 —
# conservação de momentum linear exata); escalares (rho_b, cs, c_o, c_n,
# noise, sigma_a) copiados. h_filha = α·h_mãe, offset = ε·h_mãe.
#
# Tres mudancas estruturais vs v2.3.1:
#   (a) Trigger MIGRADO de rho_rel<0.7 para sigma_a<0.85 (Violeau §3.6 —
#       particao da unidade discreta, invariante sob refinamento; mede
#       diretamente o erro do operador SPH).
#   (b) α: 0.6 → 0.35 (acoplado a ε=0.35). Feldman 2006 deriva otimo em
#       razao ε/α = 1; manter α=0.6 com ε=0.35 causou over-pack do nucleo
#       em v2.3 (rho_pos-split ≈ 1.7·rho_pre-split → a_pressure=5.015
#       sustentado → dt collapse 12×). Custo: kernel das filhas com ~12
#       vizinhos vs ~35 — degrada acuracia local mas nao viola conservacao
#       (Liu §6.5 — sub-amostragem e menos catastrofica que over-pack).
#   (c) Exigencia ESTRITA n_d=7 (todos os 6 vertices precisam passar no
#       proximity guard). Violeau §7.4.3 — apenas distribuicao hexagonal
#       simetrica preserva (i) erro de densidade < 5%, (ii) centro de massa
#       na posicao da mae, (iii) tensor de inercia local. Splits parciais
#       quebram momento angular e introduzem torque espurio (Liu §4.2.4).
use_pass_n = True  # T2g validado (baseline estavel t=100s); religado p/ refinar braços sobre morfologia que NAO fragmenta mais.
PASS_N_FREQ = 100  # iter entre checks
PASS_N_MAX_PARENTS = 25  # T2g: 100→25 — limita burst inicial (a_pressure=123 em t=3s no run anterior) e acumulo de particulas pequenas
PASS_N_SIGMA_TRIG = 0.85  # trigger v2.4: split se sigma_a < 0.85.
# Violeau §3.6: sigma_a ≈ 0.85 corresponde a ~15% erro nos operadores SPH.
# Invariante sob refinamento — gen 0/1/2 disparam pelo mesmo limiar
# (diferente do trigger rho_rel da v2.2-v2.3.1, que dependia da massa
# da particula porque rho = Sum m_j W). Mede diretamente a quantidade
# que governa a consistencia de ordem zero do SPH.
PASS_N_ALPHA = 0.35  # v2.4: 0.6 → 0.35, acoplado a ε=0.35 (Feldman 2006:
# razao ε/α = 1 minimiza erro de densidade pos-split).
PASS_N_EPSILON = 0.35  # offset filha = ε · h_mãe.
PASS_N_M_FLOOR_RATIO = (
    1.0 / 7.0
)  # T2g: gen ≤ 1 — proibe gen-2 (h=0.22dx → dt_visc 0.015× = killer dos 250×). Cap penalidade dt em ~8×. 7× de resolucao basta p/ braços de 3-4 particulas.
PASS_N_RHO_B_MIN = 0.3  # gate inferior de biomassa — exclui borda dilute da
# Gaussiana inicial (rho_b 0.05-0.3 tem rho SPH naturalmente baixo por kernel
# truncado, não gap real).
PASS_N_RHO_B_MAX = 0.7  # Pass N v2.3.1: gate superior — exclui núcleo
# e shell adjacente ao pin (rho_b ≥ 0.7). Refinamento focado na zona de
# transicao/rim ativo onde dendritos esticam. Mantido em v2.4: hard pin K.17
# torna refinamento no nucleo irrelevante (momentum nao integrado ali —
# Liu §4.5; particao da unidade no nucleo e academic).
PASS_N_PROXIMITY_MIN = 0.4  # min distância filha-vizinho em unidades de dx.
# v2.4: usado em conjunto com exigencia ESTRITA n_d=7 — se QUALQUER vertice
# falhar no proximity guard, o split inteiro e adiado para proxima call.
# Garante simetria hexagonal estrita (Violeau §7.4.3 — centro de massa
# preservado, momento angular conservado).


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
        # dt = solver.dt
        # bact = self.particles[2]

        # bact.x += dt * bact.u
        # bact.y += dt * bact.v

        # domain_width = x_max_domain - x_min_domain
        # domain_height = y_max_domain - y_min_domain

        # bact.x[:] = (bact.x - x_min_domain) % domain_width + x_min_domain
        # bact.y[:] = (bact.y - y_min_domain) % domain_height + y_min_domain

        # if solver.count % trajectory_store_interval == 0:
        #     for i in range(self.n_bact):
        #         self.all_bact_trajectories[i].append([bact.y[i], bact.x[i]])

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

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
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

        # Pass N v2.4 — Refinamento hexagonal Vacondio/Feldman com consistencia
        # SPH restaurada (Liu §3.3.3, Violeau §3.4-3.6, §7.4-7.5):
        # (a) trigger sigma_a < 0.85 (particao da unidade — Violeau §3.6),
        # (b) α = ε = 0.35 (Feldman 2006 — razao otima preserva densidade),
        # (c) gen ≤ 2 (m_floor = m₀/49),
        # (d) gate rho_b ∈ [0.3, 0.7] (exclui borda dilute e nucleo pinado),
        # (e) n_d = 7 ESTRITO (Violeau §7.4.3 — simetria hexagonal).
        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            V_0 = dx * dx  # volume inicial (referência massa)

            # Gate v2.3.1: rho_b ∈ [0.3, 0.7] — exclui borda dilute (rho_b<0.3,
            # gap não real, só kernel truncado) E núcleo/shell adjacente ao pin
            # (rho_b>0.7, onde refinamento e irrelevante porque momentum nao e
            # integrado — hard pin K.17).
            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            # Geração ≤ 2: m_floor = m₀/49 permite 2 splits sucessivos
            m_floor = V_0 * PASS_N_M_FLOOR_RATIO
            splittable_mass_mask = fluid.m > m_floor

            # Pass N v2.4 — Trigger por PARTIÇÃO DA UNIDADE (sigma_a < 0.85).
            # sigma_a = Sum_j V_j W_aj e a consistencia de ordem zero do SPH
            # (Violeau §3.6, Liu §3.3.3). Mede diretamente o erro do operador
            # SPH no nó a. INVARIANTE SOB REFINAMENTO — gen 0/1/2 disparam
            # pelo mesmo limiar. v2.2-v2.3.1 usavam rho_rel < 0.7 que dependia
            # da massa (rho = Sum m_j W) → trigger ficava enviesado em filhas
            # com m_d = m_m/7.
            sigma_a = fluid.sigma_a
            split_mask = (
                colony_mask & (sigma_a < PASS_N_SIGMA_TRIG) & splittable_mass_mask
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
                        f"(n_d=7 estrito, ε=α={eps}, σ_trig={PASS_N_SIGMA_TRIG}, gen≤2)"
                    )


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
