import csv
import numpy as np
from scipy.spatial import cKDTree
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme


def cubic_spline_w(r, h):
    # W do CubicSpline 2D (PySPH) — usado so na logica de insercao em numpy
    q = np.asarray(r) / h
    fac = 10.0 / (7.0 * np.pi * h * h)
    w = np.zeros_like(q, dtype=float)
    m1 = q <= 1.0
    m2 = (q > 1.0) & (q <= 2.0)
    w[m1] = fac * (1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3)
    w[m2] = fac * 0.25 * (2.0 - q[m2]) ** 3
    return w


def void_fraction(x, y, rho_b, dx, thresholds=(0.7, 1.0, 1.5), n_grid=200):
    """Fracao da AREA da colonia sem nenhuma particula dentro de thr*dx.

    Metrica primaria do Criterio Obrigatorio de Validacao: ao contrario de sigma_a
    (que so existe onde ha particula), esta mede o vacuo geometrico e portanto
    enxerga a quebra de cobertura espacial exigida pela Particao da Unidade
    (Violeau §3.4; consistencia de interpolacao Liu §3.3).
    """
    colony = rho_b > 0.1
    if int(np.sum(colony)) < 10:
        return {t: 0.0 for t in thresholds}
    r = np.hypot(x, y)
    R = float(np.percentile(r[colony], 99))
    if R <= 0:
        return {t: 0.0 for t in thresholds}
    g = np.linspace(-R, R, n_grid)
    GX, GY = np.meshgrid(g, g)
    # Recorta ao DOMINIO: quando a colonia passa da parede, o disco de raio R cobre
    # regiao sem particula por construcao e isso seria contado como vacuo fisico.
    inside = (
        (GX * GX + GY * GY <= R * R)
        & (np.abs(GX) <= x_max_domain)
        & (np.abs(GY) <= y_max_domain)
    )
    if not np.any(inside):
        return {t: 0.0 for t in thresholds}
    d, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX[inside], GY[inside]])
    )
    return {t: float(np.mean(d > t * dx)) for t in thresholds}


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
    "min_c_n",
    "max_c_n",
    "mean_c_n",
    "contrast_c_n",
    "mass_total",
    "pass_n_spawned",
    "min_sig_bio",
    "mean_sig_bio",
    "frac_lowsig_bio",
    "mean_sig_all",
    "frac_lowsig_all",
    "n_bio_arms",
    "n_ins_arms",
    "void_07",
    "void_10",
    "void_15",
    "a_mar_bio_med",
    "a_mar_bio_p95",
    "cs_bio_arms",
    "c_n_bio_arms",
    "biomass_total",
    "biomass_arms",
    "n_pinned",
    "frac_clump",
    "nn_median",
    "n_shift_gate",
]

# Dominio expandido 2026-08-06: em [-5,5] a colonia rompia a parede em t~51s
# (R_p99=6.40, 228 particulas alem de 4.8) — Bloqueio H / licao M-B.10.
# 261 preserva dx: 14/260 = 0.05385 vs 10/186 = 0.05376.
x_dim, y_dim = 261, 261

x_min_domain, x_max_domain = -7.0, 7.0
y_min_domain, y_max_domain = -7.0, 7.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)


mu = 0.020
gamma = 60.0
beta = 5.0
sigma = 10.0
D = 1.5e-3
D_ext = 0.08
lambda_ = 0.15
r_growth = 0.02
rho_max = 1.0
alpha_mon = 0.12

D_o = 0.04
k_o = 0.5
lambda_o = 0.05
Q0 = 5.0

D_n = 0.05
D_n_int = 1e-4
k_n = 0.5

dt_global = 0.001
total_sim_time = 50.0
print_freq = 200

# None = cada rodada tem condicao inicial propria (producao).
# Fixar um inteiro torna a rodada bit-reproduzivel — OBRIGATORIO ao comparar rotas
# (§2.5): sem semente, diferencas de ate ~9 pontos percentuais em metricas de
# amostra pequena (n_bio ~70-90) nao sao atribuiveis ao mecanismo.
SEED = None

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35

use_splitting = False

use_shift = True
SHIFT_COEFF = 0.5
SHIFT_CAP = 0.0006
SHIFT_RHO_B_MIN = 0.1

use_kgc = True
KGC_DET_MIN = 0.25

# S4 = transparencia quimica de cs (sempre ativa em equations.py)
# S5 = + filler nao consome nutriente (alavanca isolada)
FILLER_NUTRIENT_TRANSPARENT = 0

use_insert = True
INSERT_FREQ = 200
INSERT_SIGMA_TRIG = 0.85
INSERT_RHO_B_MIN = (
    0.5  # C3.4 validado: filler frozen so no nucleo estrutural (licao #31/#35)
)
INSERT_PROX = 0.7
INSERT_MAX = 100

use_wake = True
WAKE_FREQ = 100
WAKE_DISP = 1.0
WAKE_PROX = 0.7
WAKE_RHO_B_MIN = 0.05
WAKE_MAX = 150
WAKE_MODE = 2
WAKE_CLUSTER_MAX = 7
WAKE_RING_RATIO = 0.75
WAKE_MASS_BUDGET = 0.12
# Rota C (ABORTADA 2026-08-06) — realocar agar ocioso conservaria massa, MAS o agar
# tem pressao ZERO (BiomassEOS: fade_rep=fade_att=0 para rho_b<0.1), entao o buraco
# deixado pelo doador NAO cicatriza: cada doacao e uma puncao permanente no campo.
# Codigo preservado para referencia; so reativar se o agar ganhar resposta de pressao.
WAKE_RECYCLE = False
WAKE_DONOR_MARGIN = 0.6

use_pass_n = False
PASS_N_FREQ = 100
PASS_N_MAX_PARENTS = 100
PASS_N_SIGMA_TRIG = 0.95
PASS_N_ALPHA = 0.75
PASS_N_EPSILON = 0.5
PASS_N_MAX_GEN = 1
PASS_N_RHO_B_MIN = 0.15
PASS_N_RHO_B_MAX = 0.95
PASS_N_PROXIMITY_MIN = 0.4


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = (
            None  # snapshot da massa total em t=0 (preenchido em post_step)
        )
        self._pass_n_spawned_since_log = 0  # acumula spawns entre linhas de log
        self._wake_mass_added = 0.0  # so o que o WAKE adicionou (nao BiomassGrowth)

    def create_particles(self):
        if SEED is not None:
            np.random.seed(SEED)
        fluid_solid = create_initial_state(
            x_dim,
            y_dim,
            rho_max,
            dt_global,
            x_min=x_min_domain,
            x_max=x_max_domain,
            y_min=y_min_domain,
            y_max=y_max_domain,
            seed=SEED,
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
                pa.add_property("au_mar")
                pa.add_property("ax_mar")
                pa.add_property("ay_mar")
                pa.add_property("au_drag")
                pa.add_property("grad_rho_b_x")
                pa.add_property("grad_rho_b_y")
                pa.add_property("grad_rho_b_mag")
                pa.add_property("grad_cs_x")
                pa.add_property("grad_cs_y")
                pa.add_property("grad_co_x")
                pa.add_property("grad_co_y")
                pa.add_property("sigma_a")
                pa.sigma_a[:] = 1.0
                pa.add_property("shift_dC_x")
                pa.add_property("shift_dC_y")
                pa.add_property("shift_x")
                pa.add_property("shift_y")
                pa.shift_x[:] = 0.0
                pa.shift_y[:] = 0.0
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
                pa.add_property("is_filler")
                pa.is_filler[:] = 0.0
                pa.add_property("is_wake")
                pa.is_wake[:] = 0.0
                pa.add_property("x_dep")
                pa.add_property("y_dep")
                pa.x_dep[:] = pa.x[:]
                pa.y_dep[:] = pa.y[:]
                pa.add_property("gen")
                pa.gen[:] = 0.0
                pa.add_property("ax_drag")
                pa.add_property("ay_drag")
                pa.add_property("au_flag")
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
                        "sigma_a",
                        "is_wake",
                        "is_filler",
                        "c_n",
                        "m",
                        "rho",
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
            filler_nutrient_transparent=FILLER_NUTRIENT_TRANSPARENT,
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
        if self._m_initial is None:
            self._m_initial = float(np.sum(self.particles[0].m))

        # Print stats
        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            v_mag = np.sqrt(fluid.u**2 + fluid.v**2)
            max_v = np.max(v_mag)
            mean_v = np.mean(v_mag)
            n_fast = int(np.sum(v_mag > 0.1))

            a_mar = np.max(np.abs(fluid.au_mar))
            a_drag = np.max(np.abs(fluid.au_drag))
            a_total = np.max(np.sqrt(fluid.au**2 + fluid.av**2))
            ax_p = fluid.au - fluid.ax_mar - fluid.ax_drag
            ay_p = fluid.av - fluid.ay_mar - fluid.ay_drag
            a_pressure = np.max(np.sqrt(ax_p**2 + ay_p**2))
            a_flag = np.max(np.abs(fluid.au_flag))

            # Filler e quimicamente transparente: seu cs fica CONGELADO no valor
            # herdado e nao representa o campo. Incluí-lo infla mean_cs e deprime
            # contrast_cs artificialmente.
            _chem = fluid.is_filler < 0.5
            _cs_real = fluid.cs[_chem] if np.any(_chem) else fluid.cs
            min_cs = np.min(_cs_real)
            max_cs = np.max(_cs_real)
            mean_cs = np.mean(_cs_real)
            contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

            min_c_n = np.min(fluid.c_n)
            max_c_n = np.max(fluid.c_n)
            mean_c_n = np.mean(fluid.c_n)
            contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

            mass_total = float(np.sum(fluid.m))

            arms_mask = (fluid.rho_b_grown >= 0.1) & (fluid.rho_b_grown < 0.5)
            bio_mask = arms_mask & (fluid.is_filler < 0.5)
            n_bio_arms = int(np.sum(bio_mask))
            n_ins_arms = int(np.sum(arms_mask)) - n_bio_arms

            if n_bio_arms > 0:
                sig_bio = fluid.sigma_a[bio_mask]
                min_sig_bio = float(np.min(sig_bio))
                mean_sig_bio = float(np.mean(sig_bio))
                frac_lowsig_bio = float(np.mean(sig_bio < 0.85))
            else:
                min_sig_bio = mean_sig_bio = 1.0
                frac_lowsig_bio = 0.0

            if int(np.sum(arms_mask)) > 0:
                sig_all = fluid.sigma_a[arms_mask]
                mean_sig_all = float(np.mean(sig_all))
                frac_lowsig_all = float(np.mean(sig_all < 0.85))
            else:
                mean_sig_all = 1.0
                frac_lowsig_all = 0.0

            vf = void_fraction(fluid.x, fluid.y, fluid.rho_b_grown, dx)

            # Motor medido SO na biomassa real: 'a_marangoni' acima e um MAXIMO
            # (invariante #1 do §9 — max e cego ao tipico).
            # Biomassa REAL (exclui filler) e contagem de pinados: n_pinned e o
            # teste direto da maturacao dos braços em nucleo (licao #18/#41).
            # Clumping (Liu §6.4): resolvido o vacuo, o defeito remanescente e
            # sobreposicao, nao falta de particula. Mede-se pelo vizinho mais proximo.
            _rr = np.hypot(fluid.x, fluid.y)
            if np.any(fluid.rho_b_grown > 0.1):
                _Rc = float(np.percentile(_rr[fluid.rho_b_grown > 0.1], 99))
            else:
                _Rc = 1.0
            _colony = (fluid.rho_b_grown > 0.05) & (_rr > 0.3 * _Rc)
            if int(np.sum(_colony)) > 10:
                _d, _ = cKDTree(np.column_stack([fluid.x, fluid.y])).query(
                    np.column_stack([fluid.x[_colony], fluid.y[_colony]]), k=2
                )
                _nn = _d[:, 1] / dx
                frac_clump = float(np.mean(_nn < 0.5))
                nn_median = float(np.median(_nn))
            else:
                frac_clump = 0.0
                nn_median = 1.0

            # quantas particulas o ParticleShift efetivamente processa (gate de c_n)
            n_shift_gate = int(
                np.sum(
                    (fluid.rho_b_grown >= SHIFT_RHO_B_MIN)
                    & (fluid.rho_b_grown < 0.8)
                    & (fluid.c_n >= 0.6)
                )
            )

            _real = fluid.is_filler < 0.5
            _vol = fluid.m[_real] / np.maximum(fluid.rho[_real], 1e-9)
            biomass_total = float(np.sum(fluid.rho_b_grown[_real] * _vol))
            n_pinned = int(np.sum(fluid.rho_b_grown >= 0.8))

            _bio = (fluid.rho_b_grown > 0.1) & (fluid.is_filler < 0.5)
            if int(np.sum(_bio)) > 0:
                _r = np.hypot(fluid.x, fluid.y)
                _R = float(np.percentile(_r[fluid.rho_b_grown > 0.1], 99))
                _arm = _bio & (_r > 0.3 * _R)
                a_mar_bio_med = float(np.median(fluid.au_mar[_bio]))
                a_mar_bio_p95 = float(np.percentile(fluid.au_mar[_bio], 95))
                if int(np.sum(_arm)) > 0:
                    cs_bio_arms = float(np.mean(fluid.cs[_arm]))
                    c_n_bio_arms = float(np.mean(fluid.c_n[_arm]))
                    biomass_arms = float(
                        np.sum(
                            fluid.rho_b_grown[_arm]
                            * fluid.m[_arm]
                            / np.maximum(fluid.rho[_arm], 1e-9)
                        )
                    )
                else:
                    cs_bio_arms = c_n_bio_arms = biomass_arms = 0.0
            else:
                a_mar_bio_med = a_mar_bio_p95 = cs_bio_arms = c_n_bio_arms = 0.0
                biomass_arms = 0.0

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
            )
            print(
                f"sigma_a braços biomassa (n={n_bio_arms}): min={min_sig_bio:.3f} "
                f"mean={mean_sig_bio:.3f} frac<0.85={frac_lowsig_bio:.2%} | "
                f"c/ inseridas (n={n_ins_arms}): mean={mean_sig_all:.3f} "
                f"frac<0.85={frac_lowsig_all:.2%}"
            )
            print(
                f"vazio areal: >0.7dx={vf[0.7]:.2%} >1.0dx={vf[1.0]:.2%} "
                f">1.5dx={vf[1.5]:.2%}  (alvo: >1.5dx -> 0)"
            )
            print(
                f"motor na biomassa: a_mar med={a_mar_bio_med:.3f} "
                f"p95={a_mar_bio_p95:.2f} | bracos: cs={cs_bio_arms:.4f} "
                f"c_n={c_n_bio_arms:.3f}   (S0: med=3.38 p95=9.46 cs=0.353 c_n=0.552)"
            )
            print(
                f"biomassa real: total={biomass_total:.4f} bracos={biomass_arms:.4f}"
                f" | n_pinned(rho_b>=0.8)={n_pinned}"
            )
            print(
                f"empacotamento: clump(<0.5dx)={frac_clump:.1%} "
                f"nn_mediana={nn_median:.3f}dx | shift processa {n_shift_gate} part."
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
                        f"{min_sig_bio:.4f}",
                        f"{mean_sig_bio:.4f}",
                        f"{frac_lowsig_bio:.4f}",
                        f"{mean_sig_all:.4f}",
                        f"{frac_lowsig_all:.4f}",
                        n_bio_arms,
                        n_ins_arms,
                        f"{vf[0.7]:.4f}",
                        f"{vf[1.0]:.4f}",
                        f"{vf[1.5]:.4f}",
                        f"{a_mar_bio_med:.4f}",
                        f"{a_mar_bio_p95:.4f}",
                        f"{cs_bio_arms:.5f}",
                        f"{c_n_bio_arms:.4f}",
                        f"{biomass_total:.5f}",
                        f"{biomass_arms:.5f}",
                        n_pinned,
                        f"{frac_clump:.4f}",
                        f"{nn_median:.4f}",
                        n_shift_gate,
                    ]
                )
                self._pass_n_spawned_since_log = 0  # reseta após registrar

        if use_splitting:
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

        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            splittable_gen_mask = fluid.gen < PASS_N_MAX_GEN

            sigma_a = fluid.sigma_a
            split_mask = (
                colony_mask & (sigma_a < PASS_N_SIGMA_TRIG) & splittable_gen_mask
            )
            split_idx = np.where(split_mask)[0]

            if len(split_idx) > 0:
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

                    n_d = 1 + len(valid_vertices)
                    if n_d != 7:
                        continue

                    m_d = m_m / 7.0
                    h_d = alpha * h_m
                    u_d = float(u_arr[k])
                    v_d = float(v_arr[k])

                    xs = [x_m] + [p[0] for p in valid_vertices]
                    ys = [y_m] + [p[1] for p in valid_vertices]

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
                        "is_filler": [0.0] * 7,
                        "gen": [float(gen_arr[k]) + 1.0] * 7,
                        "is_wake": [0.0] * 7,
                        "x_dep": xs,
                        "y_dep": ys,
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

        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx

            void_mask = (fluid.rho_b_grown > INSERT_RHO_B_MIN) & (
                fluid.sigma_a < INSERT_SIGMA_TRIG
            )
            void_idx = np.where(void_mask)[0]

            if len(void_idx) > 0:
                positions = np.column_stack([fluid.x, fluid.y])
                tree = cKDTree(positions)
                prox = INSERT_PROX * dx  # distância de segurança entre duas partículas
                prox_sq = prox * prox  # o quadrado da distância de segurança

                angles_hex = np.arange(6) * (
                    np.pi / 3.0
                )  # angulos hexagonais 0°, 60°, 120°, 180°, 240°, 300°
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
                        # coordenadas polares
                        vx = xk + dx * cos_a[j]
                        vy = yk + dx * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox:
                            continue
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
                        "is_filler": [1.0] * n_ins,
                        "is_wake": [0.0] * n_ins,
                        "x_dep": new_x,
                        "y_dep": new_y,
                    }
                    inserted.add_particles(**data)
                    fluid.append_parray(inserted)
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_ins
                    print(
                        f"C3 inserção t={solver.t:.1f}s: {n_ins} partículas inseridas"
                    )

        if use_wake and solver.count > 0 and solver.count % WAKE_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx

            disp = np.hypot(fluid.x - fluid.x_dep, fluid.y - fluid.y_dep)
            wake_idx = np.where(
                (fluid.rho_b_grown > WAKE_RHO_B_MIN) & (disp >= WAKE_DISP * dx)
            )[0]

            # Conta so a massa que o WAKE adicionou. Antes comparava a massa TOTAL,
            # entao com a biologia ativa o BiomassGrowth sozinho estouraria o teto e
            # o wake morreria por um motivo alheio a ele.
            budget_ok = (
                self._m_initial is None
                or self._wake_mass_added < WAKE_MASS_BUDGET * self._m_initial
            )
            if not budget_ok:
                wake_idx = np.array([], dtype=int)

            if len(wake_idx) > 0:
                tree = cKDTree(np.column_stack([fluid.x, fluid.y]))
                prox = WAKE_PROX * dx
                prox_sq = prox * prox
                h0 = float(fluid.h[0])
                w_self = float(cubic_spline_w(np.array([0.0]), h0)[0])
                r_ring = WAKE_RING_RATIO * dx
                w_ring = float(cubic_spline_w(np.array([r_ring]), h0)[0])

                new_x = []
                new_y = []
                parent_idx = []
                added_pts = []

                for k in wake_idx:
                    sx = float(fluid.x_dep[k])
                    sy = float(fluid.y_dep[k])
                    fluid.x_dep[k] = fluid.x[k]  # reset: deslocamento ja consumido
                    fluid.y_dep[k] = fluid.y[k]

                    d_existing, _ = tree.query([sx, sy])
                    if d_existing < prox:  # rastro ja refluido: nao ha vazio
                        continue
                    too_close = False
                    for ax_, ay_ in added_pts:
                        if (sx - ax_) ** 2 + (sy - ay_) ** 2 < prox_sq:
                            too_close = True
                            break
                    if too_close:
                        continue

                    n_extra = 0
                    if WAKE_MODE == 2:
                        nn = tree.query_ball_point([sx, sy], 2.0 * h0)
                        if len(nn) > 0:
                            nn = np.asarray(nn, dtype=int)
                            d_nn = np.hypot(fluid.x[nn] - sx, fluid.y[nn] - sy)
                            rho_void = float(
                                np.sum(fluid.m[nn] * cubic_spline_w(d_nn, h0))
                            )
                            rho_local = float(np.mean(fluid.rho[nn]))
                            deficit = rho_local - (rho_void + m_target * w_self)
                            if deficit > 0.0:
                                n_extra = int(round(deficit / (m_target * w_ring)))
                            n_extra = max(0, min(n_extra, WAKE_CLUSTER_MAX - 1))

                    new_x.append(sx)
                    new_y.append(sy)
                    parent_idx.append(int(k))
                    added_pts.append((sx, sy))

                    for j in range(n_extra):
                        ang = 2.0 * np.pi * j / max(n_extra, 1)
                        vx = sx + r_ring * np.cos(ang)
                        vy = sy + r_ring * np.sin(ang)
                        d_ex, _ = tree.query([vx, vy])
                        if d_ex < prox:
                            continue
                        bad = False
                        for ax_, ay_ in added_pts:
                            if (vx - ax_) ** 2 + (vy - ay_) ** 2 < prox_sq:
                                bad = True
                                break
                        if bad:
                            continue
                        new_x.append(vx)
                        new_y.append(vy)
                        parent_idx.append(int(k))
                        added_pts.append((vx, vy))

                    if len(new_x) >= WAKE_MAX:
                        break

                if len(new_x) > 0:
                    n_ins = len(new_x)
                    p = np.asarray(parent_idx, dtype=int)
                    if WAKE_RECYCLE:
                        # Rota C: em vez de CRIAR particula (que adiciona massa e
                        # esbarra no teto), REALOCA agar ocioso do campo distante.
                        # Massa exatamente conservada; o buraco deixado fica fora do
                        # disco da colonia, onde C1/C2 nao sao medidos e sigma_a~1.
                        r_all = np.hypot(fluid.x, fluid.y)
                        r_col = float(np.percentile(r_all[fluid.rho_b_grown > 0.1], 99))
                        donor_pool = np.where(
                            (r_all > r_col + WAKE_DONOR_MARGIN)
                            & (fluid.rho_b_grown < 0.05)
                            & (fluid.is_filler < 0.5)
                        )[0]
                        if len(donor_pool) < n_ins:
                            n_ins = len(donor_pool)
                        if n_ins > 0:
                            # amostra ESPALHADA pelo campo distante: tirar sempre os
                            # mais distantes concentraria a depleção nos 4 cantos e
                            # criaria rarefacao (pressao tensil) na fronteira.
                            donors = np.random.choice(
                                donor_pool, size=n_ins, replace=False
                            )
                            p = p[:n_ins]
                            fluid.x[donors] = np.asarray(new_x[:n_ins])
                            fluid.y[donors] = np.asarray(new_y[:n_ins])
                            fluid.u[donors] = 0.0
                            fluid.v[donors] = 0.0
                            fluid.rho_b_grown[donors] = fluid.rho_b_grown[p]
                            fluid.cs[donors] = fluid.cs[p]
                            fluid.c_o[donors] = fluid.c_o[p]
                            fluid.c_n[donors] = fluid.c_n[p]
                            fluid.noise[donors] = fluid.noise[p]
                            fluid.is_filler[donors] = 1.0
                            fluid.is_wake[donors] = 1.0
                            fluid.x_dep[donors] = fluid.x[donors]
                            fluid.y_dep[donors] = fluid.y[donors]
                            solver.nnps.update()
                            self._pass_n_spawned_since_log += n_ins
                            print(
                                f"wake reciclagem t={solver.t:.1f}s: {n_ins} "
                                f"agar realocado (pool={len(donor_pool)})"
                            )
                    else:
                        inserted = fluid.empty_clone()
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
                            "is_filler": [1.0] * n_ins,
                            "is_wake": [1.0] * n_ins,
                            "x_dep": new_x,
                            "y_dep": new_y,
                        }
                        inserted.add_particles(**data)
                        fluid.append_parray(inserted)
                        solver.nnps.update()

                        self._wake_mass_added += n_ins * m_target
                        self._pass_n_spawned_since_log += n_ins
                        print(f"wake inserção t={solver.t:.1f}s: {n_ins} inseridas")


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
