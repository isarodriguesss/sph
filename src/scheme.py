from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
from pysph.sph.wc.basic import MomentumEquation

from .equations import (
    BiomassGrowth,
    SurfactantProductionDecay,
    SurfactantDiffusion,
    MarangoniForce,
    LinearDrag,
)


class CustomEulerStep(EulerStep):
    def stage1(
        self,
        d_idx,
        d_u,
        d_v,
        d_au,
        d_av,
        d_x,
        d_y,
        d_rho,
        d_rho_b_grown,
        d_a_rho_b_grown,
        d_cs,
        d_a_c_s,
        dt,
    ):
        d_u[d_idx] += dt * d_au[d_idx]
        d_v[d_idx] += dt * d_av[d_idx]
        d_x[d_idx] += dt * d_u[d_idx]
        d_y[d_idx] += dt * d_v[d_idx]

        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_cs[d_idx] += dt * d_a_c_s[d_idx]

        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0))
        d_cs[d_idx] = max(1e-9, d_cs[d_idx])


class MyBiomassScheme(Scheme):
    def __init__(
        self, fluids, solids, dim, mu, gamma, beta, sigma, D, lambda_, r_growth, rho_max
    ):
        self.mu = mu
        self.gamma = gamma
        self.beta = beta
        self.sigma = sigma
        self.D = D
        self.lambda_ = lambda_
        self.r_growth = r_growth
        self.rho_max = rho_max
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid", "solid"]),
            ]
        )

        equations_fluid_solid = Group(
            equations=[
                MomentumEquation(
                    dest="fluid",
                    sources=["fluid", "solid"],
                    c0=10.0,
                    alpha=self.mu,
                    beta=0.0,
                )
            ]
        )

        equation_fluid_fluid = Group(
            equations=[
                MarangoniForce(dest="fluid", sources=["fluid"], beta=self.beta),
                SurfactantDiffusion(dest="fluid", sources=["fluid"], D=self.D),
            ]
        )

        equations_pointwise = Group(
            equations=[
                BiomassGrowth(
                    dest="fluid",
                    sources=None,
                    r_growth=self.r_growth,
                    rho_max=self.rho_max,
                ),
                SurfactantProductionDecay(
                    dest="fluid", sources=None, sigma=self.sigma, lambda_=self.lambda_
                ),
                LinearDrag(dest="fluid", sources=None, gamma=self.gamma),
            ]
        )

        return [
            equations_pre,
            equations_fluid_solid,
            equation_fluid_fluid,
            equations_pointwise,
        ]

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
