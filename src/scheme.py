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
        dt,
    ):
        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_cs[d_idx] += dt * d_a_c_s[d_idx]

        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0))
        d_cs[d_idx] = max(1e-9, d_cs[d_idx])

        # K.17 — hard core pinning: rho_b >= 0.8 representa matriz EPS madura
        # (gel solido imovel). Velocidade zerada; posicao congelada.
        # cs e rho_b continuam evoluindo (producao/difusao ativas no nucleo).
        if d_rho_b_grown[d_idx] < 0.8:
            d_u[d_idx] += dt * d_au[d_idx]
            d_v[d_idx] += dt * d_av[d_idx]
            d_x[d_idx] += dt * d_u[d_idx]
            d_y[d_idx] += dt * d_v[d_idx]
            d_m[d_idx] += dt * d_am[d_idx]
        else:
            d_u[d_idx] = 0.0
            d_v[d_idx] = 0.0

        d_u[d_idx] += dt * d_au[d_idx]
        d_v[d_idx] += dt * d_av[d_idx]
        d_x[d_idx] += dt * d_u[d_idx]
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
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid"]),
                BiomassEOS(dest="fluid", sources=None, rho0=1.0, c0=self.c0),
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
                ),
                FlagellarForce(
                    dest="fluid",
                    sources=["fluid"],
                    f0=0.5,
                ),
            ],
        )

        return [equations_pre, equations_main]

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
