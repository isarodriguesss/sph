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
total_sim_time = 100.0  # C3.2 validada em t=50s (vacuo profundo preenchido, massa +2.1%, inserções auto-limitantes); estender p/ confirmar sustentabilidade
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
# ── Rota A — Fickian Particle Shifting (Xu 2009; Lind 2012) ──────────────
# Alternativa de custo-dt-ZERO ao Pass N para o problema do vacuo: redistribui
# particulas existentes (-∇C) p/ drenar sigma_a → 1 (Violeau §3.6, Liu §6.5),
# sem criar particulas pequenas que cravam o dt (lição #29). Respeita hard pin.
use_shift = False  # Rota A ESGOTADA (v1 congela #31, A.3 deixa vacuo #32). Preservada.
SHIFT_COEFF = 0.5  # D = shift_coeff·h² no δr = -D∇C (relaxacao Fickiana)
SHIFT_CAP = 0.05  # |δr| ≤ 0.05·h por passo (Lind 2012 — estabilidade)
SHIFT_RHO_B_MIN = 0.6  # A.3: gate SO interior [0.6,0.8) — exclui frontier motil

# ── Rota C — Kernel Gradient Correction (Bonet-Lok 1999; CSPM) ───────────
# Corrige o ∇cs da Marangoni (consistencia 1a ordem, Liu §3.3) nos braços
# sub-resolvidos SEM mover particula → nao congela (#31) nem depende de
# separar vacuo↔frontier por rho_b (#32). Auto-gateia pelo det(M): bulk
# det≈1 (sem correcao), rim/vacuo det<1 (corrige). Custo de dt ZERO.
use_kgc = True  # Rota C ativa (KEEP — gradiente da Marangoni válido nas pontas)
KGC_DET_MIN = 0.25  # fallback p/ identidade se det(M)<0.25 (|L|≲4×, anti-spike)

# ── Rota C3 — Inserção de partículas no vácuo (em espaçamento dx) ─────────
# Preenche o vácuo central (evacuação dinâmica núcleo↔rim) INSERINDO partículas
# frescas na rede dx — NÃO dividindo (Vacondio over-packa, #27). h=h0 e
# espaçamento=dx → dt INTACTO (#29) e SEM over-pack (#27). Inseridas herdam
# rho_b alto → pinadas → congelam e preenchem estável. Não move partícula
# existente (escapa #31/#32). Ancorado em particle insertion/packing (Liu §6.5).
use_insert = True  # Rota C3 ativa (lever sob teste; KGC mantido)
INSERT_FREQ = 200  # iter entre inserções (~4s com dt~0.02)
INSERT_RHO_TRIG = 0.6  # C3.2: 0.7→0.6 — fillar so vacuo PROFUNDO (consistencia
# genuinamente quebrada, σ_a severo <0.85 ~ rho/rho0<0.6, Violeau §3.6). A franja
# 0.6-0.7 e under-density leve (nao e buraco real) — inseri-la so inflava massa.
INSERT_RHO_B_MIN = (
    0.5  # só vácuo estrutural (núcleo/junção/braço interior); poupa tips (<0.5)
)
INSERT_PROX = (
    0.7  # spot vazio se nenhum existente a < 0.7·dx (insere só em buraco real)
)
INSERT_MAX = 100  # máx inserções por call (limita crescimento de massa)

use_pass_n = False  # Rota A: Pass N DESLIGADO (preservado) p/ isolar o efeito do shifting. Religar so apos avaliar Rota A.
PASS_N_FREQ = 100  # iter entre checks
PASS_N_MAX_PARENTS = 100  # v2.5: 25→100 — enchimento agressivo. Viavel pq α=0.75 mantem h grande (lição #29: dt e min-h, nao contagem)
PASS_N_SIGMA_TRIG = 0.95  # v2.5: 0.85→0.95 — gatilho PREEMPTIVO (agir ao 1o sinal de estiramento, 5% de perda, nao 15%).
# Violeau §3.6: sigma_a ≈ 0.85 corresponde a ~15% erro nos operadores SPH.
# Invariante sob refinamento — gen 0/1/2 disparam pelo mesmo limiar
# (diferente do trigger rho_rel da v2.2-v2.3.1, que dependia da massa
# da particula porque rho = Sum m_j W). Mede diretamente a quantidade
# que governa a consistencia de ordem zero do SPH.
PASS_N_ALPHA = 0.75  # v2.5: 0.35 → 0.75 — h_filha ~ h_mãe PRESERVA dt-por-h
# (penalidade ~1.5× vs ~85× em v2.4; lição #29 — dt e min-h). ε MANTIDO em 0.35
# (decoupling deliberado de Feldman ε/α=1): aceita-se over-pack inicial, confiando
# na EOS coesiva (tension_ratio=0.30) + Monaghan (alpha_mon=0.12, K.23) p/ relaxar
# (Liu §6.5). RISCO lição #27: ε/α=0.47 < v2.3 (0.58) que travou — critério #2
# (dt avg t>30s ≥ 0.5× inicial) e o teste de relaxacao. Fallback: ε→0.5 ou α→0.6.
PASS_N_EPSILON = 0.35  # offset filha = ε · h_mãe (mantido v2.4 — evita abortar splits no domínio empacotado).
PASS_N_M_FLOOR_RATIO = (
    1.0 / 7.0
)  # T2g: gen ≤ 1 — proibe gen-2 (h=0.22dx → dt_visc 0.015× = killer dos 250×). Cap penalidade dt em ~8×. 7× de resolucao basta p/ braços de 3-4 particulas.
PASS_N_RHO_B_MIN = 0.15  # v2.5: 0.3→0.15 — alarga zona ativa p/ fechar buracos
# em quase toda a colonia (enchimento agressivo).
PASS_N_RHO_B_MAX = 0.95  # v2.5: 0.7→0.95 — INCLUI a junção núcleo-dendrito
# (vácuo crítico). NOTA: 0.8-0.95 cai na zona pinada (hard pin K.17 rho_b≥0.8) —
# refinar lá adiciona vizinhos com força descartada (Liu §4.5, lição #27);
# aposta do usuario: fechar o vácuo da junção vale o risco. Monitorar a_pressure
# na borda do núcleo.
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

        # ── Rota C3 — Inserção de partículas no vácuo (em espaçamento dx) ─────
        # Preenche o vácuo central INSERINDO partículas frescas na rede dx
        # (NÃO split → sem over-pack #27; h=h0 → dt intacto #29). Inseridas
        # herdam campos da mãe-vácuo (pinadas se rho_b alto → congelam estável).
        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx  # massa alvo = volume da rede inicial
            rho_rel = fluid.rho / 1.0  # rho0 = 1.0

            # Vácuo na zona ESTRUTURAL interior (rho_b alto, fora do frontier
            # motil <0.5): núcleo/junção/braço interior com rho/rho0 < trigger.
            void_mask = (fluid.rho_b_grown > INSERT_RHO_B_MIN) & (
                rho_rel < INSERT_RHO_TRIG
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
