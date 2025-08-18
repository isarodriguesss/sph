from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
from .equations import (
    BiomassGrowth,
    SurfactantEquation,
    LinearDrag,
    InterpolateVelocity,
    ViscousForce,
    MarangoniForce,
    ActiveSpreadingForce,
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
        self,
        fluids,
        solids,
        others,
        dim,
        mu,
        gamma,
        beta,
        sigma,
        D,
        lambda_,
        r_growth,
        rho_max,
        rho0,
        c0,
        P_active
    ):
        self.mu = mu
        self.gamma = gamma
        self.beta = beta
        self.sigma = sigma
        self.D = D
        self.lambda_ = lambda_
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.rho0 = rho0
        self.c0 = c0
        self.P_active = P_active
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid", "solid"]),
            ],
            real=False,
        )

        equations_main = Group(
            equations=[
                MarangoniForce(dest="fluid", sources=["fluid"], beta=self.beta),
                ViscousForce(dest="fluid", sources=["fluid"], mu=self.mu),
                LinearDrag(dest="fluid", sources=None, gamma=self.gamma),
                ActiveSpreadingForce(dest="fluid", sources=["fluid"], P_active=self.P_active),
                SurfactantEquation(
                    dest="fluid",
                    sources=["fluid"],
                    D=self.D,
                    sigma=self.sigma,
                    lambda_=self.lambda_,
                ),
                BiomassGrowth(
                    dest="fluid",
                    sources=None,
                    r_growth=self.r_growth,
                    rho_max=self.rho_max,
                ),
            ],
        )

        equations_interp = Group(
            equations=[InterpolateVelocity(dest="bact", sources=["fluid"])]
        )

        return [equations_pre, equations_main, equations_interp]

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
