from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.basic_equations import SummationDensity
from .equations import BiomassGrowth, SurfactantProductionDecay, SurfactantForceAndDrag
from pysph.sph.equation import Group

# O integrador customizado agora integra TODAS as propriedades (padrão e customizadas)
class CustomEulerStep(EulerStep):
    def stage1(self, d_idx,
               d_u, d_v, d_w, d_au, d_av, d_aw,
               d_x, d_y, d_z, d_rho, d_arho,
               dt,
               # Adicione aqui o estado e a taxa de suas propriedades
               d_rho_b_grown, d_a_rho_b_grown,
               d_c_s, d_a_c_s):
        
        # Integração padrão de Euler
        d_u[d_idx] += dt * d_au[d_idx]
        d_v[d_idx] += dt * d_av[d_idx]
        d_w[d_idx] += dt * d_aw[d_idx]

        d_x[d_idx] += dt * d_u[d_idx]
        d_y[d_idx] += dt * d_v[d_idx]
        d_z[d_idx] += dt * d_w[d_idx]

        d_rho[d_idx] += dt * d_arho[d_idx]

        # Integração das suas propriedades: ESTADO += dt * TAXA
        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_c_s[d_idx] += dt * d_a_c_s[d_idx]

        # É uma boa prática aplicar limites/clipagem aqui, após a atualização
        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0)) # rho_max é 1.0
        d_c_s[d_idx] = max(1e-9, d_c_s[d_idx])

class MyBiomassScheme(Scheme):
    def __init__(self, fluids, solids, dim, rho_max, r_growth, sigma, lambda_, beta, gamma, D, mu):
        super(MyBiomassScheme, self).__init__(fluids, solids, dim)
        self.rho_max = rho_max
        self.r_growth = r_growth
        self.sigma = sigma
        self.lambda_ = lambda_
        self.beta = beta
        self.gamma = gamma
        self.D = D
        self.mu = mu

    def get_equations(self):
        equations_interaction = Group(equations=[
            SummationDensity(dest='fluid', sources=['fluid']),
        ])

        equations_pointwise = Group(equations=[
            BiomassGrowth(dest='fluid', sources=None, r_growth=self.r_growth, rho_max=self.rho_max),
            SurfactantProductionDecay(dest='fluid', sources=None, sigma=self.sigma, lambda_=self.lambda_),
            SurfactantForceAndDrag(dest='fluid', sources=None, beta=self.beta, gamma=self.gamma),
        ])

        return [equations_interaction, equations_pointwise]

    def get_integrator(self):
        # O integrador usa o nosso stepper customizado
        return EulerIntegrator(fluid=CustomEulerStep())