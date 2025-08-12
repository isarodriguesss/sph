# src/scheme.py

from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

# Importações
from pysph.sph.basic_equations import SummationDensity
from pysph.sph.wc.basic import TaitEOS, MomentumEquation
from .equations import (
    BiomassGrowth,
    SurfactantProductionDecay,
    SurfactantDiffusion,
    MarangoniForce,
    LinearDrag,
    InterpolateVelocity,
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
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        # Grupo 1: Pré-cálculos para densidade e pressão.
        # A opção `real=False` informa ao PySPH que este grupo não calcula
        # as acelerações finais, então ele não deve zerar as propriedades de taxa.
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid", "solid"]),
                TaitEOS(
                    dest="fluid", sources=None, rho0=self.rho0, c0=self.c0, gamma=7.0
                ),
            ],
            real=False,
        )

        # Grupo 2: O grupo principal. TODAS as equações que calculam
        # forças e taxas de mudança de escalares vão aqui.
        equations_main = Group(
            equations=[
                # Equações de Força (contribuem para au, av)
                MomentumEquation(
                    dest="fluid",
                    sources=["fluid", "solid"],
                    c0=self.c0,
                    alpha=self.mu,
                    beta=0.0,
                ),
                MarangoniForce(dest="fluid", sources=["fluid"], beta=self.beta),
                LinearDrag(dest="fluid", sources=None, gamma=self.gamma),
                # Equações de Taxa de Escalares (contribuem para a_cs, a_rho_b_grown)
                SurfactantDiffusion(dest="fluid", sources=["fluid"], D=self.D),
                SurfactantProductionDecay(
                    dest="fluid", sources=None, sigma=self.sigma, lambda_=self.lambda_
                ),
                BiomassGrowth(
                    dest="fluid",
                    sources=None,
                    r_growth=self.r_growth,
                    rho_max=self.rho_max,
                ),
            ]
        )

        # Grupo 3: Interpolação para as bactérias (separado, pois tem um destino diferente)
        equations_interp = Group(
            equations=[InterpolateVelocity(dest="bact", sources=["fluid"])]
        )

        return [equations_pre, equations_main, equations_interp]

    def get_integrator(self):
        # Seu integrador customizado está correto e é necessário.
        return EulerIntegrator(fluid=CustomEulerStep())
