from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.basic_equations import SummationDensity
from .equations import BiomassGrowth, SurfactantProductionDecay, SurfactantForceAndDrag
from pysph.sph.equation import Group
import numpy as np

class MyBiomassScheme(Scheme):
    def __init__(self, fluids, solids, dim, rho_max, r_growth, sigma, lambda_, beta, gamma, D, mu, periodic_domain=None, **kw):
        super(MyBiomassScheme, self).__init__(fluids, solids, dim, **kw)
        self.rho_max = rho_max
        self.r_growth = r_growth
        self.sigma = sigma
        self.lambda_ = lambda_
        self.beta = beta
        self.gamma = gamma
        self.D = D
        self.mu = mu
        self.periodic_domain = periodic_domain

    def get_equations(self):
        equations = [
            Group(equations=[
                SummationDensity(dest='fluid', sources=['fluid']),
            ]),
            Group(equations=[
                BiomassGrowth(dest='fluid', sources=None, r_growth=self.r_growth, rho_max=self.rho_max)
            ]),
            Group(equations=[
                SurfactantProductionDecay(dest='fluid', sources=None, sigma=self.sigma, lambda_=self.lambda_)
            ]),
            Group(equations=[
                SurfactantForceAndDrag(dest='fluid', sources=None, beta=self.beta, gamma=self.gamma)
            ]),
        ]
        return equations

    def get_integrator(self):
        return EulerIntegrator(fluid=EulerStep())