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
sigma = 20.0  # Pass T2d: T1 σ=5 → 20 (acelera saturacao em cs_max=0.5); |F_mar| final depende de β·|∇cs|, nao de σ
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
total_sim_time = 50.0  # Validacao T2d ate t=100s antes de aplicar Pass N (pre-requisito morfologia base)
print_freq = 200

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35  # EOS: B = 1.5²/7 ≈ 0.32 (repulsão suave, atração ~0.1)

use_splitting = False

# Pass N v2 — Refinamento hexagonal Vacondio 2013 / Feldman 2006
# (Soleimani 2017 §3.2.6). Substitui 1 mãe por 7 filhas (1 centro + 6 vértices
# hexagonais), massa dividida IGUALMENTE (m_filha = m_mãe/7), velocidade
# IDÊNTICA herdada (conserva momentum linear e angular), escalares (rho_b, cs,
# c_o, c_n, noise) copiados; smoothing length h_filha = α·h_mãe; offset hexagonal
# r = ε·h_mãe; α = ε = 0.6 (Feldman 2006, erro de densidade < 5%).
# Trigger: V_a = m_a / rho_a > V_limit (partícula esticada, não buraco angular).
# Aplica em toda a colônia (rho_b > 0.05) — núcleo + braços + transição.
use_pass_n = True
PASS_N_FREQ = 100  # iter entre checks (mais agressivo que v1=200)
PASS_N_MAX_PARENTS = 100  # máx mães por call (×~7 filhas = ~700 novas/call)
PASS_N_RHO_TRIG = 0.7  # trigger relativo: split se rho_a/rho_0 < 0.7.
# v2.1 usou V_a > 1.5·dx² (absoluto), que para
# gen ≥ 1 requer rho < 0.1·rho_0 (raro) — trigger
# ficava cego a partículas em "gap critico"
# (visualizadas como vermelho/azul no painel 3
# do plot.py). Critério relativo casa exatamente
# com a definição visual de gap (rho < 0.7·rho_0).
PASS_N_ALPHA = 0.6  # h_filha = α · h_mãe (Feldman 2006)
PASS_N_EPSILON = 0.35  # offset filha = ε · h_mãe — REDUZIDO de 0.6 (Feldman puro)
# para 0.35 evita overlap com vizinhos a ~dx (causou
# instabilidade max_v=8.25 em call 3 de v2-inicial).
# Vertices agora a 0.63·dx do centro.
PASS_N_M_FLOOR_RATIO = 1.0 / 49.0  # permite gen ≤ 2: m > m₀/49 ainda splittable.
# v2-inicial usou 1/7 → bloqueava após 1 gen → splits
# cessaram em t=7.6s mesmo com braços alongando.
PASS_N_RHO_B_MIN = 0.3  # gate inferior de biomassa — exclui borda dilute da
# Gaussiana inicial (rho_b 0.05-0.3 tem rho SPH naturalmente baixo por kernel
# truncado, não gap real).
PASS_N_RHO_B_MAX = 0.7  # Pass N v2.3.1 (Opção 2): gate superior — exclui núcleo
# e shell adjacente ao pin (rho_b ≥ 0.7) do refinamento.
# Justificativa: v2.3 desbloqueou gen-1 mas todos os splits
# se concentraram no centro (relaxacao inicial da Gaussiana
# criou gaps no nucleo), causando over-pack do core e
# a_pressure=5.0 sustentado → dt collapse 12× → simulação
# travou em t=13s/iter 7200. Restringindo a rho_b ∈ [0.3, 0.7]
# (zona de transicao/rim ativo), o refinamento so vai disparar
# quando dendritos do rim esticarem (t > 30s). Trade-off:
# nucleo continua perdendo vizinhos com o tempo sem reposicao
# — aceita-se este custo para manter pace temporal saudavel
# e focar refinamento onde a morfologia esta evoluindo.
PASS_N_PROXIMITY_MIN = 0.4  # min distância filha-vizinho em unidades de dx;
# vertices que cairiam < 0.4·dx de partícula existente
# são descartadas (mass redistribuída em N_actual < 7).


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

        # Pass N v2.1 — Refinamento hexagonal Vacondio/Feldman + correções:
        # (a) gen ≤ 2 (m_floor = m₀/49), (b) ε=0.35 (evita overlap c/ vizinhos),
        # (c) gate rho_b > 0.3 (exclui borda da Gaussiana), (d) proximity guard
        # (filhas a < 0.4·dx de vizinho existente são descartadas; massa
        # redistribuída entre N_actual filhas para preservar m_total da mãe).
        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            rho_safe = np.maximum(fluid.rho, 1e-6)  # evita div0 (planktônicas)
            rho0 = 1.0  # densidade de referência SPH

            V_0 = dx * dx  # volume inicial (referência massa)

            # Gate v2.3.1: rho_b ∈ [0.3, 0.7] — exclui borda dilute (rho_b<0.3,
            # gap não real, só kernel truncado) E núcleo/shell adjacente ao pin
            # (rho_b>0.7, onde gaps são gerados pela migração do rim mas o
            # refinamento ali causou over-pack do centro em v2.3 → dt collapse).
            # Refinamento agora focado na zona de transição/rim ativo, onde
            # dendritos esticam e geram gaps morfologicamente significativos.
            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            # Geração ≤ 2: m_floor = m₀/49 permite 2 splits sucessivos
            m_floor = V_0 * PASS_N_M_FLOOR_RATIO
            splittable_mass_mask = fluid.m > m_floor

            # Trigger RELATIVO: rho_a / rho_0 < PASS_N_RHO_TRIG (= 0.7).
            # Casa exatamente com a definição de "gap crítico" visualizada
            # no painel 3 do plot.py. Independente da geração — gen 0, 1, 2
            # disparam pelo mesmo critério (não sofre da restrição artificial
            # do trigger absoluto V_a, que ficava 7× mais estrito para gen 1).
            rho_rel = rho_safe / rho0
            split_mask = (
                colony_mask & (rho_rel < PASS_N_RHO_TRIG) & splittable_mass_mask
            )
            split_idx = np.where(split_mask)[0]

            if len(split_idx) > 0:
                # Prioriza mais esvaziadas (menor rho_rel) se exceder limite
                if len(split_idx) > PASS_N_MAX_PARENTS:
                    order = np.argsort(rho_rel[split_idx])
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

                    # n_d = 1 centro + N vértices válidos
                    n_d = 1 + len(valid_vertices)

                    # Skip se sobraram poucas filhas (ganho de resolução baixo)
                    if n_d < 4:
                        continue

                    # Massa REDISTRIBUÍDA entre n_d filhas (preserva m_total mãe)
                    m_d = m_m / n_d
                    h_d = alpha * h_m
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
                        "m": [m_d] * n_d,
                        "h": [h_d] * n_d,
                        "rho": [float(rho_arr[k])] * n_d,
                        "rho_b_grown": [float(rhob_arr[k])] * n_d,
                        "cs": [float(cs_arr[k])] * n_d,
                        "c_o": [float(co_arr[k])] * n_d,
                        "c_n": [float(cn_arr[k])] * n_d,
                        "u": [u_d] * n_d,
                        "v": [v_d] * n_d,
                        "noise": [float(noise_arr[k])] * n_d,
                    }
                    daughters.add_particles(**data)
                    mothers_used.append(int(split_idx[k]))
                    n_daughters_total += n_d

                if mothers_used:
                    fluid.append_parray(daughters)
                    fluid.remove_particles(np.asarray(mothers_used, dtype=np.uint32))
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_daughters_total
                    avg_d = n_daughters_total / len(mothers_used)
                    print(
                        f"Pass N v2.3 t={solver.t:.1f}s: "
                        f"{len(mothers_used)} mães → {n_daughters_total} filhas "
                        f"(<n_d>={avg_d:.1f}, ε={eps}, ρ_trig={PASS_N_RHO_TRIG}, gen≤2)"
                    )


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
