from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.basic_equations import SummationDensity
from .equations import BiomassGrowth, SurfactantProductionDecay, SurfactantForceAndDrag
from pysph.sph.equation import Group
import numpy as np

# *************************************************************
# CORREÇÃO: stage1 DEVE SER UM MÉTODO PYTHON NORMAL!
# O PySPH se encarrega da Cythonização internamente.
# Remova 'cdef inline void' e 'noexcept:'.
# *************************************************************
class CustomEulerStep(EulerStep):
    # O PySPH vai inferir os tipos e Cythonizar este método.
    def stage1(self, d_idx,
               d_u, d_v, d_w,
               d_au, d_av, d_aw,
               d_x, d_y, d_z,
               d_rho, d_arho,
               dt,
               # Suas propriedades personalizadas e suas derivadas
               d_rho_b_grown, d_d_rho_b_grown,
               d_c_s, d_d_c_s):

        # Implemente as integrações padrão do EulerStep diretamente
        d_u[d_idx] += dt * d_au[d_idx]
        d_v[d_idx] += dt * d_av[d_idx]
        d_w[d_idx] += dt * d_aw[d_idx]

        d_x[d_idx] += dt * d_u[d_idx]
        d_y[d_idx] += dt * d_v[d_idx]
        d_z[d_idx] += dt * d_w[d_idx]

        d_rho[d_idx] += dt * d_arho[d_idx]

        # Agora, adicione suas integrações de propriedades customizadas
        d_rho_b_grown[d_idx] += dt * d_d_rho_b_grown[d_idx]
        d_c_s[d_idx] += dt * d_d_c_s[d_idx]

        # Opcional: Adicionar prints de depuração se ainda precisar.
        # if d_idx == 0:
        #    print(f"DEBUG CustomEulerStep.stage1: rho_b_grown[{d_idx}] AFTER_INTEGRATION: {d_rho_b_grown[d_idx]:.4f}")
        #    print(f"DEBUG CustomEulerStep.stage1: c_s[{d_idx}] AFTER_INTEGRATION: {d_c_s[d_idx]:.4f}")


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

    def initialize_properties(self, particle_arrays, clean=True):
        super(MyBiomassScheme, self).initialize_properties(particle_arrays, clean)

        for pa in particle_arrays:
            if pa.name == 'fluid':
                if not hasattr(pa, 'rho_b_grown'):
                    print("WARNING: rho_b_grown not found in particle array at scheme initialize_properties")
                    pa.add_property('rho_b_grown', default=0.0)
                if not hasattr(pa, 'c_s'):
                    print("WARNING: c_s not found in particle array at scheme initialize_properties")
                    pa.add_property('c_s', default=0.0)

                print(f"DEBUG scheme.py initialize_properties (BEFORE SOLVER): Max rho_b_grown: {np.max(pa.rho_b_grown):.4f}")
                print(f"DEBUG scheme.py initialize_properties (BEFORE SOLVER): Max c_s: {np.max(pa.c_s):.4f}")

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
        return EulerIntegrator(fluid=CustomEulerStep())