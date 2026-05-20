from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
from pysph.sph.wc.basic import MomentumEquation
from .equations import (
    BiomassEOS,
    BiomassGrowth,
    FlagellarForce,
    MarangoniForce,
    OxigenConsumption,
    SurfactantEquation,
    LinearDrag,
    ViscousForce,
    BiomassGradient,
)


class CustomEulerStep(EulerStep):
    def stage1(
        self,
        d_idx,
        d_m,
        d_am,
        d_u,
        d_v,
        d_au,
        d_av,
        d_x,
        d_y,
        d_rho_b_grown,
        d_a_rho_b_grown,
        d_cs,
        d_a_c_s,
        d_c_o,
        d_a_c_o,
        d_c_n,
        d_a_c_n,
        dt,
    ):
        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_cs[d_idx] += dt * d_a_c_s[d_idx]
        d_c_o[d_idx] += dt * d_a_c_o[d_idx]
        d_c_n[d_idx] += dt * d_a_c_n[d_idx]

        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0))
        d_cs[d_idx] = max(1e-9, d_cs[d_idx])
        d_c_o[d_idx] = max(0.0, min(d_c_o[d_idx], 1.0))
        d_c_n[d_idx] = max(1e-9, min(d_c_n[d_idx], 1.0))

        # M-B.9b: revert ao pin quimico hard de M-B.8 (c_n<0.4) — M-B.9a
        # demonstrou que drag suave + coesao baixa gera fragmentacao
        # distribuida (lição §22). Esta versao testa se reforcar coesao
        # (tension_ratio 0.02→0.08 em BiomassEOS) reduz a fragmentacao
        # interna pos-t=50s observada com pin hard. Mantem K.17 mecanico.
        if d_rho_b_grown[d_idx] >= 0.8 or d_c_n[d_idx] < 0.6:
            d_u[d_idx] = 0.0
            d_v[d_idx] = 0.0
        else:
            d_u[d_idx] += dt * d_au[d_idx]
            d_v[d_idx] += dt * d_av[d_idx]

        d_x[d_idx] += dt * d_u[d_idx]  # x avança com u (que é 0 se pinnado)
        d_y[d_idx] += dt * d_v[d_idx]
        d_m[d_idx] += dt * d_am[d_idx]

        # vmax = 5.0

        # v = np.sqrt(d_u[d_idx]**2 + d_v[d_idx]**2)
        # if v > vmax:
        #     scale = vmax / v
        #     d_u[d_idx] *= scale
        #     d_v[d_idx] *= scale


class MyBiomassScheme(Scheme):
    def __init__(
        self,
        fluids,
        solids,
        dim,
        mu,
        gamma,
        beta,
        sigma,
        D,
        D_ext,
        lambda_,
        r_growth,
        rho_max,
        c0=10.0,
        alpha_mon=0.5,
        p0=0.0,
        D_o=1e-3,
        k_o=0.5,
        lambda_o=0.05,
        Q0=5e-4,
        D_n=0.02,
        D_n_int=1e-4,
        k_n=0.5,
    ):
        self.mu = mu
        self.gamma = gamma
        self.beta = beta
        self.sigma = sigma
        self.D = D
        self.D_ext = D_ext
        self.lambda_ = lambda_
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.c0 = c0
        self.alpha_mon = alpha_mon
        self.p0 = p0
        self.D_o = D_o
        self.k_o = k_o
        self.lambda_o = lambda_o
        self.Q0 = Q0
        self.D_n = D_n
        self.D_n_int = D_n_int
        self.k_n = k_n
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid"]),
                BiomassEOS(
                    dest="fluid",
                    sources=None,
                    rho0=1.0,
                    c0=self.c0,
                    tension_ratio=0.30,  # M-B.9b: 0.02 → 0.08 (4x) — coesao reforcada
                    # contra fragmentacao interna pos-t=50s (lição §15/§22). B_tension
                    # sobe de 0.00069 → 0.00275. Risco: re-aparecer halo nas baias.
                ),
            ],
            real=False,
        )

        equations_main = Group(
            equations=[
                # MomentumEquation com Monaghan artificial viscosity forte
                # (alpha=0.5) para manter continuidade no braço dendrítico.
                # Pressão do EOS já faz repulsão E coesão (p<0 → atração).
                MomentumEquation(
                    dest="fluid",
                    sources=["fluid"],
                    c0=self.c0,
                    alpha=self.alpha_mon,
                    beta=0.0,
                ),
                BiomassGrowth(
                    dest="fluid",
                    sources=None,
                    r_growth=self.r_growth,
                    rho_max=self.rho_max,
                ),
                # BiomassGradient DEVE vir antes de MarangoniForce
                # (Marangoni usa grad_rho_b_mag como gate de interface)
                BiomassGradient(dest="fluid", sources=["fluid"]),
                MarangoniForce(dest="fluid", sources=["fluid"], beta=self.beta),
                ViscousForce(dest="fluid", sources=["fluid"], mu=self.mu),
                LinearDrag(
                    dest="fluid",
                    sources=None,
                    gamma_base=self.gamma,
                    gamma_mature=self.gamma
                    * 1.5,  # K.16c: revertido a baseline K.15; pinning via edge_fade invertido
                ),
                SurfactantEquation(
                    dest="fluid",
                    sources=["fluid"],
                    D=self.D,
                    D_ext=self.D_ext,
                    sigma=self.sigma,
                    lambda_=self.lambda_,
                    lambda_ext_ratio=5.0,  # decaimento no agar mantem halo finito (Trinschek-like)
                    k_consume=0.0,  # Pass T1: sumidouro removido — saturacao agora via (1-cs/cs_max)
                    cs_max=0.5,  # Pass T1: Γ_max do painel (b) Trinschek 2018 — alvo de saturacao
                ),
                OxigenConsumption(
                    dest="fluid",
                    sources=["fluid"],
                    D_n=self.D_n,
                    D_n_int=self.D_n_int,
                    k_n=self.k_n,
                ),
                # Pass M-A (Frente 6): osmolitos secretados pelas bacterias.
                # |∇c_o| dispara influxo de massa (van't Hoff) — pontas incham,
                # baias estagnam. DEVE preceder FlagellarForce (que usa cs).
                # OsmolyteProduction(
                #     dest="fluid",
                #     sources=["fluid"],
                #     D_o=self.D_o,
                #     k_o=self.k_o,
                #     lambda_o=self.lambda_o,
                #     Q0=self.Q0,
                # ),
                FlagellarForce(
                    dest="fluid",
                    sources=["fluid"],
                    f0=3.0,  # K.22: 0.5→3.0 — restaura ignição nas pontas (alvo §8: f0~γ·v_term/2=3)
                ),
            ],
        )

        return [equations_pre, equations_main]

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
