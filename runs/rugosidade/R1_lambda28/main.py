import csv
import json
import os
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline, WendlandQuintic
from pysph.solver.solver import Solver

from src.particles import SEED_MODE, create_initial_state
from src.scheme import MyBiomassScheme


# Baseline P2 (E11 + filler doador). Historico e justificativa de cada valor: CLAUDE.md §7/§9.

x_dim, y_dim = 261, 261
x_min_domain, x_max_domain = -7.0, 7.0
y_min_domain, y_max_domain = -7.0, 7.0
dx = (x_max_domain - x_min_domain) / (x_dim - 1)

KERNEL = "cubic"  # "wendland_c2" com H_FACTOR=1.92 em teste (P2W, licao #88)
H_FACTOR = 1.8

mu = 0.020
gamma = 60.0
alpha_mon = 0.12
c0 = 0.35

beta = 5.0
FLAG_F0 = 3.0
FLAG_GATE_LO = 0.2
FLAG_GATE_HI = 0.6

sigma = 11.1
HILL_K = 0.25
CS_CEILING = 0.5
D = 1.5e-3
D_ext = 0.08
lambda_ = 0.15
LAMBDA_BIO_RATIO = 2.0

r_growth = 0.02
rho_max = 1.0

D_n = 0.05
D_n_int = 1e-4
k_n = 0.5
k_src = 0.3

k_col = 0.03
COL_CS_MIN = 0.3
COL_FILLER_DONOR = 1.0  # P2: filler conta como doador. 0.0 = E11.

total_sim_time = 50.0
print_freq = 200

NOISE_AMP = 0.6
SEED = 20260806

use_shift = True
SHIFT_COEFF = 0.5
SHIFT_CAP = 0.0006
SHIFT_RHO_B_MIN = 0.1

use_kgc = True
KGC_DET_MIN = 0.25

use_insert = True
INSERT_FREQ = 200
INSERT_SIGMA_TRIG = 0.85
INSERT_RHO_B_MIN = 0.5
INSERT_PROX = 0.7
INSERT_MAX = 100

use_wake = True
WAKE_FREQ = 100
WAKE_DISP = 1.0
WAKE_PROX = 0.7
WAKE_RHO_B_MIN = 0.05
WAKE_MAX = 150
WAKE_CLUSTER_MAX = 7
WAKE_RING_RATIO = 0.75
WAKE_MASS_BUDGET = 0.12
WAKE_SEG = False
WAKE_SEG_MAX = 6
RASTRO_W = 5.0
RASTRO_R_MIN = 0.7
RASTRO_COS_FRENTE = 1.0
RASTRO_SEG = True
RASTRO_BAIA = 2.0
RASTRO_BAIA_CHAMADA = False
RASTRO_HIST = 10
RASTRO_R_MOV = 0.4
RASTRO_PASSO_MOV = 0.5
RASTRO_V_MOV = 0.0  # > 0: criterio de VELOCIDADE (u/s) no lugar do passo por chamada
RASTRO_PONTA = 1.0  # > 0: calota eliptica de semi-eixo RASTRO_PONTA*w atras de p1 (1.0 = semicirculo)
RASTRO_PONTA_LINK = 6.0  # > 0 (dx): calota so no lider mais externo de cada grupo; 0 = todos
RASTRO_PONTA_RECOBRE = True  # True: lider com calota converte largura plena ate RASTRO_PONTA*w atras de p0
FILLER_RHO_B_FLOOR = 0.4
FILLER_CS_CONDUZ = 1.0
FILLER_CS_D = 0.0
FILLER_CS_D_INTERNO = 0.0
FILLER_CS_LAMBDA = 0.5
AGAR_CS_LAMBDA = 0.125

# pass-l-aprovado: Pass L (rugosidade) desbloqueado em 2026-09-17 — CLAUDE.md §12 e
# docs/PLANO_RUGOSIDADE.md. Rede triangular de pilares no plano; geometria reancorada no P2R23.
use_pilares = True
PILAR_A = 10.0  # diametro do pilar em dx (~ largura do braco: 9-10.5 dx no P2R23)
PILAR_LAMBDA = 28.0  # espacamento da rede triangular em dx (~1.5 contatos por braco ate t=50)
PILAR_R_EXCL = 1.2  # nenhum centro dentro deste raio (nao perturbar inoculo nem juncao)
PILAR_ROT = 0.0  # rotacao da rede, em graus
PILAR_K = 0.005  # forca de contato (Monaghan & Kajtar 2009); calibrada no plano, secao 4
PILAR_HW = 0.5  # h do kernel de contato, em dx (alcance 2*h_w = 1 dx: zero na rede inicial)


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
    "pass_n_spawned",  # nome historico: particulas inseridas (insert + wake) desde a ultima linha
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
    "c_n_junc",
    "rho_b_junc",
    "rho_b_dip",
    "rho_b_dip_r",
    "n_junc_bio",
    "a_press_max",  # |au - Marangoni - arrasto - flagelo - viscosa|: pressao EOS + visc. artificial
    "a_press_med",  # mediana do mesmo termo na colonia (rho_b >= 0.1 ou filler)
    # pass-l-aprovado: rugosidade (plano secao 6). Zeradas quando use_pilares=False.
    "n_pen",  # fluido com centro a < R_p - 0.5 dx de um centro de pilar — criterio de aborto 2
    # `n_pen` e a definicao PRE-REGISTRADA (docs/CRITERIOS_RUGOSIDADE.md §6, plano §6) e fica
    # como esta. Mas ela le 0 em TODOS os frames do V1, que rodou sem a guarda de deposicao e
    # tinha 6 particulas dentro do raio nominal (5 filler, d/R_p entre 0.971 e 0.998): elas
    # ficam nos vaos do anel de superficie, a mais de R_p - 0.5 dx do centro. `n_pen_sup` e o
    # indicador precoce que teria pego aquele vazamento; `n_pen` segue sendo o aborto.
    "n_pen_sup",  # fluido com centro a < R_p (cruzou a superficie nominal)
    "n_contato",  # fluido a < 1 dx de uma particula de pilar
    "a_rep_max",  # |a| da forca de contato: maximo
    "a_rep_med",  # e mediana, sobre os em contato
]


def kernel_w(r, h):
    q = np.asarray(r) / h
    w = np.zeros_like(q, dtype=float)
    if KERNEL == "wendland_c2":
        m = q < 2.0
        w[m] = (7.0 / (4.0 * np.pi * h * h)) * (1.0 - 0.5 * q[m]) ** 4 * (2.0 * q[m] + 1.0)
        return w
    fac = 10.0 / (7.0 * np.pi * h * h)
    m1 = q <= 1.0
    m2 = (q > 1.0) & (q <= 2.0)
    w[m1] = fac * (1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3)
    w[m2] = fac * 0.25 * (2.0 - q[m2]) ** 3
    return w


# pass-l-aprovado: mascara de pilar na metrica areal — Pass L desbloqueado em 2026-09-17
# (CLAUDE.md §1/§12); a area do pilar e solida e nao pode contar como vazio da colonia.
def void_fraction(x, y, rho_b, dx, thresholds=(0.7, 1.0, 1.5), n_grid=200,
                  pilar_c=None, pilar_rp=0.0):
    """Fracao da area da colonia sem particula a menos de thr*dx (criterio C1, §2.5).

    Sem pilares (`pilar_c=None`) o resultado e identico ao anterior.
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
    # recorta ao dominio: fora dele nao ha particula por construcao
    inside = (
        (GX * GX + GY * GY <= R * R)
        & (np.abs(GX) <= x_max_domain)
        & (np.abs(GY) <= y_max_domain)
    )
    if pilar_c is not None and len(pilar_c) > 0:
        d_pil, _ = cKDTree(pilar_c).query(np.column_stack([GX.ravel(), GY.ravel()]))
        inside &= d_pil.reshape(GX.shape) > pilar_rp
    if not np.any(inside):
        return {t: 0.0 for t in thresholds}
    d, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX[inside], GY[inside]])
    )
    return {t: float(np.mean(d > t * dx)) for t in thresholds}


def filler_data(fluid, new_x, new_y, parent, is_wake):
    n = len(new_x)
    p = np.asarray(parent, dtype=int)
    return {
        "x": new_x,
        "y": new_y,
        "m": [dx * dx] * n,
        "h": list(fluid.h[p]),
        "rho": list(fluid.rho[p]),
        "rho_b_grown": list(np.maximum(fluid.rho_b_grown[p], FILLER_RHO_B_FLOOR)),
        "cs": list(fluid.cs[p]),
        "c_n": list(fluid.c_n[p]),
        "u": [0.0] * n,
        "v": [0.0] * n,
        "noise": list(fluid.noise[p]),
        "is_filler": [1.0] * n,
        "is_wake": [is_wake] * n,
        "x_dep": new_x,
        "y_dep": new_y,
    }


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = None
        self._inseridas_desde_log = 0
        self._wake_mass_added = 0.0
        self._rastro_x = None
        self._rastro_t = None
        self._rastro_y = None
        self._rastro_hist = []
        self._col_tree = None
        self._rastro_ponta = None
        # pass-l-aprovado: geometria dos pilares, usada pelas mascaras de metrica e pelas
        # guardas de deposicao. Fica None quando use_pilares=False (comportamento anterior).
        self._pilar_c = None
        self._pilar_rp = 0.0
        self._pilar_pts = None

    def create_particles(self):
        fluid_solid = create_initial_state(
            x_dim,
            y_dim,
            rho_max,
            x_min=x_min_domain,
            x_max=x_max_domain,
            y_min=y_min_domain,
            y_max=y_max_domain,
            seed=SEED,
            h_factor=H_FACTOR,
        )

        for pa in fluid_solid:
            if pa.name == "fluid":
                pa.add_property("noise")
                pa.noise[:] = (
                    1.0
                    + NOISE_AMP * np.sin(SEED_MODE * np.arctan2(pa.y, pa.x))
                    + 0.01 * np.random.rand(len(pa.x))
                )
                for prop in (
                    "dt_force",
                    "dt_cfl",
                    "au_mar",
                    "ax_mar",
                    "ay_mar",
                    "au_drag",
                    "ax_drag",
                    "ay_drag",
                    "au_flag",
                    "ax_flag",
                    "ay_flag",
                    "ax_vis",
                    "ay_vis",
                    # pass-l-aprovado: instrumentacao da forca de contato (plano secao 6).
                    # Zeradas por ForcaContornoPilar.initialize; com use_pilares=False a
                    # equacao nao entra no scheme e os arrays ficam identicamente nulos.
                    "ax_rep",
                    "ay_rep",
                    "grad_rho_b_x",
                    "grad_rho_b_y",
                    "grad_rho_b_mag",
                    "grad_cs_x",
                    "grad_cs_y",
                    "rho_b_smooth",
                    "rho_b_w2",
                    "shift_dC_x",
                    "shift_dC_y",
                    "shift_x",
                    "shift_y",
                    "Mxx",
                    "Mxy",
                    "Myx",
                    "Myy",
                    "Lxx",
                    "Lxy",
                    "Lyx",
                    "Lyy",
                    "sigma_a",
                    "is_filler",
                    "is_wake",
                    "x_dep",
                    "y_dep",
                ):
                    pa.add_property(prop)
                pa.sigma_a[:] = 1.0
                pa.Lxx[:] = 1.0
                pa.Lyy[:] = 1.0
                pa.x_dep[:] = pa.x[:]
                pa.y_dep[:] = pa.y[:]
                pa.add_output_arrays(
                    [
                        "rho_b_grown",
                        "cs",
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

        if use_pilares:
            fluid_solid.append(self._faz_pilares(fluid_solid[0]))

        return fluid_solid

    def _faz_pilares(self, fluid):
        """Rede triangular de pilares: os proprios pontos da grade viram parede estatica.

        pass-l-aprovado: Pass L desbloqueado (CLAUDE.md §12). Mesma rede, mesma massa dx^2,
        entao a densidade vista pelos vizinhos nao muda (plano, secao 2).
        """
        lam = PILAR_LAMBDA * dx
        rp = 0.5 * PILAR_A * dx
        ang = np.radians(PILAR_ROT)
        n = int(np.ceil(2.0 * max(abs(x_min_domain), x_max_domain) / lam)) + 2
        i, j = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1))
        cx = lam * (i + 0.5 * j).ravel()
        cy = lam * (np.sqrt(3.0) / 2.0) * j.ravel()
        cx, cy = (cx * np.cos(ang) - cy * np.sin(ang), cx * np.sin(ang) + cy * np.cos(ang))
        borda = 2.0 * dx
        dentro_dom = (
            (cx > x_min_domain + rp + borda)
            & (cx < x_max_domain - rp - borda)
            & (cy > y_min_domain + rp + borda)
            & (cy < y_max_domain - rp - borda)
        )
        fora_inoculo = np.hypot(cx, cy) >= PILAR_R_EXCL
        cx, cy = cx[dentro_dom & fora_inoculo], cy[dentro_dom & fora_inoculo]

        tree = cKDTree(np.column_stack([fluid.x, fluid.y]))
        idx = np.unique(np.concatenate(tree.query_ball_point(np.column_stack([cx, cy]), rp)))
        idx = np.asarray(idx, dtype=int)
        pilar = fluid.extract_particles(idx)
        pilar.set_name("pilar")
        pilar.u[:] = 0.0
        pilar.v[:] = 0.0
        pilar.add_output_arrays(["m", "rho"])
        # pass-l-aprovado: duas representacoes da mesma geometria — o disco analitico
        # (centro, raio) para as mascaras de metrica, e o conjunto de particulas para as
        # guardas de deposicao, que reusam o mesmo teste de proximidade do fluido.
        self._pilar_c = np.column_stack([cx, cy])
        self._pilar_rp = rp
        self._pilar_pts = cKDTree(
            np.column_stack([np.asarray(pilar.x), np.asarray(pilar.y)])
        )
        fluid.remove_particles(idx)

        out = os.path.join(self.output_dir, "pilares.json")
        os.makedirs(self.output_dir, exist_ok=True)
        with open(out, "w") as fp:
            json.dump(
                {
                    "centros": np.column_stack([cx, cy]).tolist(),
                    "raio": rp,
                    "lambda": lam,
                    "rotacao_graus": PILAR_ROT,
                    "r_exclusao": PILAR_R_EXCL,
                    "n_particulas": int(len(idx)),
                },
                fp,
            )
        print(f"pilares: {len(cx)} centros, {len(idx)} particulas, raio {rp / dx:.1f} dx, "
              f"lambda {lam / dx:.0f} dx, area {100 * len(idx) / (len(idx) + len(fluid.x)):.1f}%")
        return pilar

    def create_scheme(self):
        return MyBiomassScheme(
            fluids=["fluid"],
            solids=["solid"],
            pilares=["pilar"] if use_pilares else (),
            pilar_k=PILAR_K,
            pilar_hw=PILAR_HW * dx,
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
            D_n=D_n,
            D_n_int=D_n_int,
            k_n=k_n,
            k_src=k_src,
            k_col=k_col,
            cs_max=CS_CEILING,
            col_cs_min=COL_CS_MIN,
            col_filler_donor=COL_FILLER_DONOR,
            hill_k=HILL_K,
            lambda_bio_ratio=LAMBDA_BIO_RATIO,
            flag_gate_lo=FLAG_GATE_LO,
            flag_gate_hi=FLAG_GATE_HI,
            flag_f0=FLAG_F0,
            use_shift=use_shift,
            shift_coeff=SHIFT_COEFF,
            shift_cap=SHIFT_CAP,
            shift_rho_b_min=SHIFT_RHO_B_MIN,
            use_kgc=use_kgc,
            kgc_det_min=KGC_DET_MIN,
            filler_cs_conduz=FILLER_CS_CONDUZ,
            filler_cs_D=FILLER_CS_D,
            filler_cs_D_interno=FILLER_CS_D_INTERNO,
            filler_cs_lambda=FILLER_CS_LAMBDA,
            agar_cs_lambda=AGAR_CS_LAMBDA,
        )

    def create_solver(self):
        kernel = WendlandQuintic(dim=2) if KERNEL == "wendland_c2" else CubicSpline(dim=2)
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
        if solver.count % print_freq == 0:
            self._registra(solver)
        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            self._insere_vacuo(solver)
        if use_wake and solver.count > 0 and solver.count % WAKE_FREQ == 0:
            self._deposita_rastro(solver)
            if RASTRO_W > 0.0:
                self._alarga_rastro(solver)

    def _registra(self, solver):
        fluid = self.particles[0]
        rb = fluid.rho_b_grown
        viva = fluid.is_filler < 0.5

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
        # `a_pressure` (historica) inclui o flagelo e a viscosa; esta e so a pressao da EOS
        a_press = np.hypot(
            ax_p - fluid.ax_flag - fluid.ax_vis - fluid.ax_rep,
            ay_p - fluid.ay_flag - fluid.ay_vis - fluid.ay_rep,
        )
        corpo = (rb >= 0.1) | (fluid.is_filler > 0.5)
        a_press_max = float(np.max(a_press))
        a_press_med = float(np.median(a_press[corpo])) if np.any(corpo) else 0.0

        # pass-l-aprovado: instrumentacao da rugosidade (plano secao 6). `n_pen` e o criterio
        # de aborto 2 dos criterios pre-registrados: fluido dentro do corpo solido do pilar.
        n_pen = 0
        n_pen_sup = 0
        n_contato = 0
        a_rep_max = 0.0
        a_rep_med = 0.0
        if self._pilar_c is not None and self._pilar_pts is not None:
            xy = np.column_stack([fluid.x, fluid.y])
            d_centro, _ = cKDTree(self._pilar_c).query(xy)
            n_pen = int(np.sum(d_centro < self._pilar_rp - 0.5 * dx))
            n_pen_sup = int(np.sum(d_centro < self._pilar_rp))
            d_sup, _ = self._pilar_pts.query(xy)
            contato = d_sup < dx
            n_contato = int(np.sum(contato))
            a_rep = np.hypot(fluid.ax_rep, fluid.ay_rep)
            a_rep_max = float(np.max(a_rep))
            a_rep_med = float(np.median(a_rep[contato])) if n_contato else 0.0

        # filler tem cs congelado: fora das estatisticas de cs (licao #39)
        cs_viva = fluid.cs[viva] if np.any(viva) else fluid.cs
        min_cs = np.min(cs_viva)
        max_cs = np.max(cs_viva)
        mean_cs = np.mean(cs_viva)
        contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

        min_c_n = np.min(fluid.c_n)
        max_c_n = np.max(fluid.c_n)
        mean_c_n = np.mean(fluid.c_n)
        contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

        mass_total = float(np.sum(fluid.m))

        arms_mask = (rb >= 0.1) & (rb < 0.5)
        bio_mask = arms_mask & viva
        n_bio_arms = int(np.sum(bio_mask))
        n_ins_arms = int(np.sum(arms_mask & (fluid.is_filler > 0.5)))

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

        vf = void_fraction(
            fluid.x, fluid.y, rb, dx, pilar_c=self._pilar_c, pilar_rp=self._pilar_rp
        )

        rr = np.hypot(fluid.x, fluid.y)
        Rc = float(np.percentile(rr[rb > 0.1], 99)) if np.any(rb > 0.1) else 1.0
        colony = (rb > 0.05) & (rr > 0.3 * Rc)
        if int(np.sum(colony)) > 10:
            d, _ = cKDTree(np.column_stack([fluid.x, fluid.y])).query(
                np.column_stack([fluid.x[colony], fluid.y[colony]]), k=2
            )
            nn = d[:, 1] / dx
            frac_clump = float(np.mean(nn < 0.5))
            nn_median = float(np.median(nn))
        else:
            frac_clump = 0.0
            nn_median = 1.0

        # juncao nucleo-braco medida so nas portadoras vivas (licao #62)
        carrier = viva & (rb > 1e-6)
        ju = (rr >= 0.4) & (rr < 1.2)
        c_n_junc = float(np.mean(fluid.c_n[ju])) if np.any(ju) else 0.0
        juc = ju & carrier
        rho_b_junc = float(np.percentile(rb[juc], 90)) if np.any(juc) else 0.0
        n_junc_bio = int(np.sum(ju & viva & (rb > 0.1)))

        th = np.arctan2(fluid.y, fluid.x)
        arm = (rr > 0.55 * Rc) & (rr < 0.75 * Rc) & (rb > 0.1)
        rho_b_dip = 1.0
        rho_b_dip_r = 0.0
        if int(np.sum(arm)) > 20:
            hist, edges = np.histogram(th[arm], bins=72, range=(-np.pi, np.pi))
            step = max(0.1, 0.05 * Rc)
            for b in np.argsort(hist)[-4:]:
                c = 0.5 * (edges[b] + edges[b + 1])
                cone = np.abs(((th - c + np.pi) % (2 * np.pi)) - np.pi) < np.deg2rad(9)
                for r0 in np.arange(0.15 * Rc, 0.75 * Rc, step):
                    anel = cone & (np.abs(rr - r0) < step)
                    if not np.any(anel):
                        continue
                    ac = anel & carrier
                    v = float(np.max(rb[ac])) if np.any(ac) else 0.0
                    if v < rho_b_dip:
                        rho_b_dip = v
                        rho_b_dip_r = float(r0)

        n_shift_gate = int(np.sum((rb >= SHIFT_RHO_B_MIN) & (rb < 0.8) & (fluid.c_n >= 0.6)))

        vol = fluid.m[viva] / np.maximum(fluid.rho[viva], 1e-9)
        biomass_total = float(np.sum(rb[viva] * vol))
        n_pinned = int(np.sum(rb >= 0.8))

        bio = (rb > 0.1) & viva
        if int(np.sum(bio)) > 0:
            arm_bio = bio & (rr > 0.3 * Rc)
            a_mar_bio_med = float(np.median(fluid.au_mar[bio]))
            a_mar_bio_p95 = float(np.percentile(fluid.au_mar[bio], 95))
            if int(np.sum(arm_bio)) > 0:
                cs_bio_arms = float(np.mean(fluid.cs[arm_bio]))
                c_n_bio_arms = float(np.mean(fluid.c_n[arm_bio]))
                biomass_arms = float(
                    np.sum(rb[arm_bio] * fluid.m[arm_bio] / np.maximum(fluid.rho[arm_bio], 1e-9))
                )
            else:
                cs_bio_arms = c_n_bio_arms = biomass_arms = 0.0
        else:
            a_mar_bio_med = a_mar_bio_p95 = cs_bio_arms = c_n_bio_arms = 0.0
            biomass_arms = 0.0

        print("-" * 50)
        print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
        print(f"Velocidade Máx: {max_v:.4f} | Contraste CS: {contrast_cs:.4f}")
        print(f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} | massa: {mass_total:.2f}")
        print(
            f"sigma_a braços (n={n_bio_arms}): min={min_sig_bio:.3f} mean={mean_sig_bio:.3f} "
            f"frac<0.85={frac_lowsig_bio:.2%} | c/ filler: mean={mean_sig_all:.3f}"
        )
        print(f"vazio areal: >0.7dx={vf[0.7]:.2%} >1.0dx={vf[1.0]:.2%} >1.5dx={vf[1.5]:.2%}")
        print(
            f"motor na biomassa: a_mar med={a_mar_bio_med:.3f} p95={a_mar_bio_p95:.2f} | "
            f"biomassa={biomass_total:.4f} bracos={biomass_arms:.4f} | pinadas={n_pinned}"
        )
        print(f"clump(<0.5dx)={frac_clump:.1%} nn_mediana={nn_median:.3f}dx")
        print(
            f"Acelerações: Marangoni {a_mar:.2f} | Drag {a_drag:.2f} | Flagelo {a_flag:.2f} | "
            f"Pressão EOS máx {a_press_max:.2e} med {a_press_med:.2e} | Total {a_total:.2f}"
        )
        if self._pilar_c is not None:
            alerta = "  <<< PENETRACAO (criterio de aborto 2)" if n_pen else ""
            print(
                f"pilares: n_pen={n_pen} (sup {n_pen_sup}) n_contato={n_contato} | "
                f"a_rep máx={a_rep_max:.2e} med={a_rep_med:.2e}{alerta}"
            )

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
                    self._inseridas_desde_log,
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
                    f"{c_n_junc:.4f}",
                    f"{rho_b_junc:.4f}",
                    f"{rho_b_dip:.4f}",
                    f"{rho_b_dip_r:.3f}",
                    n_junc_bio,
                    f"{a_press_max:.4e}",
                    f"{a_press_med:.4e}",
                    n_pen,
                    n_pen_sup,
                    n_contato,
                    f"{a_rep_max:.4e}",
                    f"{a_rep_med:.4e}",
                ]
            )
        self._inseridas_desde_log = 0

    def _adiciona(self, solver, data):
        fluid = self.particles[0]
        novas = fluid.empty_clone()
        novas.add_particles(**data)
        fluid.append_parray(novas)
        solver.nnps.update()
        self._inseridas_desde_log += len(data["x"])

    # pass-l-aprovado: guarda geometrica das duas vias de deposicao (insert e wake). A
    # cKDTree do fluido nao enxerga os pilares, que vivem em array proprio — sem isto o
    # filler nasce dentro do solido (docs/PLANO_RUGOSIDADE.md secao 5).
    def _livre_pilar(self, px, py, prox):
        """False se (px, py) esta a menos de `prox` de uma particula de pilar."""
        if self._pilar_pts is None:
            return True
        d, _ = self._pilar_pts.query([px, py])
        return d >= prox

    def _insere_vacuo(self, solver):
        fluid = self.particles[0]
        void_idx = np.where(
            (fluid.rho_b_grown > INSERT_RHO_B_MIN) & (fluid.sigma_a < INSERT_SIGMA_TRIG)
        )[0]
        if len(void_idx) == 0:
            return

        tree = cKDTree(np.column_stack([fluid.x, fluid.y]))
        prox = INSERT_PROX * dx
        prox_sq = prox * prox
        ang = np.arange(6) * (np.pi / 3.0)
        cos_a = np.cos(ang)
        sin_a = np.sin(ang)

        new_x, new_y, parent = [], [], []
        for k in void_idx:
            xk = float(fluid.x[k])
            yk = float(fluid.y[k])
            for j in range(6):
                vx = xk + dx * cos_a[j]
                vy = yk + dx * sin_a[j]
                if not self._livre_pilar(vx, vy, prox):
                    continue
                d_existing, _ = tree.query([vx, vy])
                if d_existing < prox:
                    continue
                if any((vx - ax) ** 2 + (vy - ay) ** 2 < prox_sq for ax, ay in zip(new_x, new_y)):
                    continue
                new_x.append(vx)
                new_y.append(vy)
                parent.append(int(k))
            if len(new_x) >= INSERT_MAX:
                break

        if new_x:
            self._adiciona(solver, filler_data(fluid, new_x, new_y, parent, 0.0))
            print(f"inserção t={solver.t:.1f}s: {len(new_x)} partículas inseridas")

    def _deposita_rastro(self, solver):
        fluid = self.particles[0]
        m_target = dx * dx

        disp = np.hypot(fluid.x - fluid.x_dep, fluid.y - fluid.y_dep)
        wake_idx = np.where((fluid.rho_b_grown > WAKE_RHO_B_MIN) & (disp >= WAKE_DISP * dx))[0]
        # orcamento conta so a massa que o wake adicionou, nao o crescimento
        if self._wake_mass_added >= WAKE_MASS_BUDGET * self._m_initial:
            return
        if len(wake_idx) == 0:
            return

        tree = cKDTree(np.column_stack([fluid.x, fluid.y]))
        prox = WAKE_PROX * dx
        prox_sq = prox * prox
        h0 = float(fluid.h[0])
        w_self = float(kernel_w(np.array([0.0]), h0)[0])
        r_ring = WAKE_RING_RATIO * dx
        w_ring = float(kernel_w(np.array([r_ring]), h0)[0])

        new_x, new_y, parent = [], [], []

        def livre(px, py):
            if not self._livre_pilar(px, py, prox):
                return False
            d_ex, _ = tree.query([px, py])
            if d_ex < prox:
                return False
            return not any((px - ax) ** 2 + (py - ay) ** 2 < prox_sq for ax, ay in zip(new_x, new_y))

        for k in wake_idx:
            sx = float(fluid.x_dep[k])
            sy = float(fluid.y_dep[k])
            fluid.x_dep[k] = fluid.x[k]
            fluid.y_dep[k] = fluid.y[k]

            if WAKE_SEG:
                ex = float(fluid.x[k])
                ey = float(fluid.y[k])
                n_seg = min(int(np.hypot(ex - sx, ey - sy) / dx), WAKE_SEG_MAX)
                for j in range(1, n_seg + 1):
                    f = j / (n_seg + 1.0)
                    px = sx + f * (ex - sx)
                    py = sy + f * (ey - sy)
                    if livre(px, py):
                        new_x.append(px)
                        new_y.append(py)
                        parent.append(int(k))

            if not livre(sx, sy):
                continue

            # completa a densidade do vazio com um anel de ate WAKE_CLUSTER_MAX-1 particulas
            n_extra = 0
            nn = tree.query_ball_point([sx, sy], 2.0 * h0)
            if len(nn) > 0:
                nn = np.asarray(nn, dtype=int)
                d_nn = np.hypot(fluid.x[nn] - sx, fluid.y[nn] - sy)
                rho_void = float(np.sum(fluid.m[nn] * kernel_w(d_nn, h0)))
                rho_local = float(np.mean(fluid.rho[nn]))
                deficit = rho_local - (rho_void + m_target * w_self)
                if deficit > 0.0:
                    n_extra = int(round(deficit / (m_target * w_ring)))
                n_extra = max(0, min(n_extra, WAKE_CLUSTER_MAX - 1))

            new_x.append(sx)
            new_y.append(sy)
            parent.append(int(k))

            for j in range(n_extra):
                ang = 2.0 * np.pi * j / max(n_extra, 1)
                vx = sx + r_ring * np.cos(ang)
                vy = sy + r_ring * np.sin(ang)
                if not livre(vx, vy):
                    continue
                new_x.append(vx)
                new_y.append(vy)
                parent.append(int(k))

            if len(new_x) >= WAKE_MAX:
                break

        if new_x:
            self._adiciona(solver, filler_data(fluid, new_x, new_y, parent, 1.0))
            self._wake_mass_added += len(new_x) * m_target
            print(f"wake inserção t={solver.t:.1f}s: {len(new_x)} inseridas")

    def _alarga_rastro(self, solver):
        fluid = self.particles[0]
        rb = fluid.rho_b_grown
        fil = fluid.is_filler > 0.5
        r = np.hypot(fluid.x, fluid.y)
        r99 = float(np.percentile(r[(rb >= 0.1) | fil], 99))
        r_min = RASTRO_R_MIN * r99
        lider = np.where((~fil) & (rb >= 0.1) & (r >= r_min))[0]
        if RASTRO_R_MOV > 0.0 and self._rastro_x is not None:
            n_prev = len(self._rastro_x)
            cand = np.where((~fil[:n_prev]) & (rb[:n_prev] >= 0.1)
                            & (r[:n_prev] >= RASTRO_R_MOV * r99) & (r[:n_prev] < r_min))[0]
            ddx = fluid.x[cand] - self._rastro_x[cand]
            ddy = fluid.y[cand] - self._rastro_y[cand]
            passo = np.hypot(ddx, ddy)
            radial = (ddx * fluid.x[cand] + ddy * fluid.y[cand]) / np.maximum(r[cand], 1e-12)
            if RASTRO_V_MOV > 0.0 and self._rastro_t is not None:
                dt_call = max(float(solver.t) - self._rastro_t, 1e-9)
                ok = (passo / dt_call) >= RASTRO_V_MOV
            else:
                ok = passo >= RASTRO_PASSO_MOV * dx
            mov = cand[ok & (radial >= 0.5 * passo)]
            lider = np.concatenate([lider, mov])
        agar = np.where((~fil) & (rb <= 1e-12))[0]
        if len(lider) == 0 or len(agar) == 0:
            return

        self._rastro_ponta = None
        if RASTRO_PONTA > 0.0 and RASTRO_PONTA_LINK > 0.0:
            tl = cKDTree(np.column_stack([fluid.x[lider], fluid.y[lider]]))
            n_g, grupo = connected_components(
                tl.sparse_distance_matrix(tl, RASTRO_PONTA_LINK * dx), directed=False
            )
            rl = r[lider]
            self._rastro_ponta = set()
            for g in range(n_g):
                m = np.where(grupo == g)[0]
                self._rastro_ponta.add(int(lider[m[np.argmax(rl[m])]]))
        tree = cKDTree(np.column_stack([fluid.x[agar], fluid.y[agar]]))
        if RASTRO_SEG:
            if RASTRO_BAIA > 0.0:
                col = np.where(fil | (rb > 1e-12))[0]
                self._col_tree = (cKDTree(np.column_stack([fluid.x[col], fluid.y[col]])), col)
            mae = self._rastro_segmento(fluid, lider, agar, tree)
            self._rastro_x = fluid.x.copy()
            self._rastro_t = float(solver.t)
            self._rastro_y = fluid.y.copy()
            if RASTRO_BAIA > 0.0:
                self._rastro_hist.append((self._rastro_x, self._rastro_y))
                self._rastro_hist = self._rastro_hist[-RASTRO_HIST:]
            if not mae:
                return
            alvo = np.fromiter(mae.keys(), dtype=int)
            origem = np.fromiter(mae.values(), dtype=int)
            self._converte_rastro(fluid, alvo, origem, rb)
            print(f"rastro t={solver.t:.1f}s: {len(alvo)} convertidas ({len(lider)} lideres)")
            return
        viz = tree.query_ball_point(
            np.column_stack([fluid.x[lider], fluid.y[lider]]), RASTRO_W * dx
        )
        mae = {}
        for k, v in zip(lider, viz):
            if RASTRO_COS_FRENTE < 1.0:
                vk = np.hypot(fluid.u[k], fluid.v[k])
                if vk < 1e-12:
                    continue
                v = np.asarray(v, dtype=int)
                ddx = fluid.x[agar[v]] - fluid.x[k]
                ddy = fluid.y[agar[v]] - fluid.y[k]
                cos = (ddx * fluid.u[k] + ddy * fluid.v[k]) / (np.hypot(ddx, ddy) * vk + 1e-12)
                v = v[cos <= RASTRO_COS_FRENTE]
            for j in v:
                mae.setdefault(int(agar[j]), int(k))
        if not mae:
            return

        alvo = np.fromiter(mae.keys(), dtype=int)
        origem = np.fromiter(mae.values(), dtype=int)
        self._converte_rastro(fluid, alvo, origem, rb)
        print(f"rastro t={solver.t:.1f}s: {len(alvo)} convertidas ({len(lider)} lideres)")

    def _converte_rastro(self, fluid, alvo, origem, rb):
        fluid.is_filler[alvo] = 1.0
        fluid.is_wake[alvo] = 1.0
        fluid.rho_b_grown[alvo] = np.maximum(rb[origem], FILLER_RHO_B_FLOOR)
        fluid.x_dep[alvo] = fluid.x[alvo]
        fluid.y_dep[alvo] = fluid.y[alvo]

    def _rastro_segmento(self, fluid, lider, agar, tree):
        mae = {}
        if self._rastro_x is None:
            return mae
        w = RASTRO_W * dx
        n_prev = len(self._rastro_x)
        novos = []
        for k in lider:
            if k >= n_prev:
                continue
            p0 = np.array([self._rastro_x[k], self._rastro_y[k]])
            p1 = np.array([fluid.x[k], fluid.y[k]])
            d = p1 - p0
            L2 = float(d @ d)
            if L2 < (0.1 * dx) ** 2:
                continue
            tapa = RASTRO_PONTA > 0.0 and (self._rastro_ponta is None or k in self._rastro_ponta)
            atras = RASTRO_PONTA * w if (tapa and RASTRO_PONTA_RECOBRE) else 0.0
            if atras > 0.0:
                raio = 0.5 * np.sqrt(L2) + w + atras
            else:
                raio = 0.5 * np.sqrt(L2) + w
            cand = np.asarray(tree.query_ball_point(0.5 * (p0 + p1), raio), dtype=int)
            if len(cand) == 0:
                continue
            q = np.column_stack([fluid.x[agar[cand]], fluid.y[agar[cand]]])
            s = ((q - p0) @ d) / L2
            if atras > 0.0:
                t = np.clip(s, -atras / np.sqrt(L2), 1.0)
            else:
                t = np.clip(s, 0.0, 1.0)
            dist = np.hypot(q[:, 0] - (p0[0] + t * d[0]), q[:, 1] - (p0[1] + t * d[1]))
            if tapa:
                a = RASTRO_PONTA * w
                u = np.clip((1.0 - s) * np.sqrt(L2), 0.0, a)
                lim = w * np.sqrt(u * (2.0 * a - u)) / a
            else:
                lim = w
            sel = cand[(s <= 1.0) & (dist <= lim)]
            if RASTRO_BAIA > 0.0 and len(sel):
                sel = self._guarda_baia(fluid, k, p1, agar, sel, w, novos)
                novos.extend(zip(fluid.x[agar[sel]], fluid.y[agar[sel]]))
            for j in sel:
                mae.setdefault(int(agar[j]), int(k))
        return mae

    def _guarda_baia(self, fluid, k, p1, agar, sel, w, novos):
        ctree, col = self._col_tree
        g = RASTRO_BAIA * dx
        q = np.column_stack([fluid.x[agar[sel]], fluid.y[agar[sel]]])
        perto = np.unique(
            np.concatenate([np.asarray(v, dtype=int) for v in ctree.query_ball_point(q, g)] + [np.zeros(0, int)])
        )
        c = np.column_stack([fluid.x[col[perto]], fluid.y[col[perto]]])
        if novos and RASTRO_BAIA_CHAMADA:
            c = np.vstack([c, np.asarray(novos)])
        if len(c) == 0:
            return sel
        pts = [np.array([hx[k], hy[k]]) for hx, hy in self._rastro_hist if k < len(hx)] + [p1]
        d_poly = np.hypot(c[:, 0] - pts[0][0], c[:, 1] - pts[0][1])
        for a, b in zip(pts[:-1], pts[1:]):
            e = b - a
            e2 = float(e @ e)
            tt = np.clip(((c - a) @ e) / e2, 0.0, 1.0) if e2 > 0.0 else np.zeros(len(c))
            d_poly = np.minimum(d_poly, np.hypot(c[:, 0] - a[0] - tt * e[0], c[:, 1] - a[1] - tt * e[1]))
        viz = c[d_poly > w + 0.5 * dx]
        if len(viz) == 0:
            return sel
        d_viz, _ = cKDTree(viz).query(q)
        return sel[d_viz > g]


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
