from pysph.sph.equation import Equation
import numpy as np

class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, **kw):
        super(BiomassGrowth, self).__init__(dest, sources, **kw)
        self.r_growth = r_growth
        self.rho_max = rho_max

    # O loop lê o estado (d_rho_b_grown) e escreve na taxa (d_a_rho_b_grown)
    def loop(self, d_idx, d_rho_b_grown, d_a_rho_b_grown):
        current_rho_b = d_rho_b_grown[d_idx]
        
        # A lógica da equação logística está correta
        rate_of_change = self.r_growth * current_rho_b * (1 - current_rho_b / self.rho_max)

        # Armazena o resultado na propriedade de TAXA ('a_rho_b_grown')
        d_a_rho_b_grown[d_idx] = rate_of_change

        if d_idx == 0:
            print(f"DEBUG_GROWTH_LOOP (Partícula 0): lendo rho_b={current_rho_b:.4f}, calculando taxa={rate_of_change:.6e}")

class SurfactantProductionDecay(Equation):
    def __init__(self, dest, sources, sigma, lambda_, **kw):
        super(SurfactantProductionDecay, self).__init__(dest, sources, **kw)
        self.sigma = sigma
        self.lambda_ = lambda_
    
    # Aplicando o mesmo padrão para a produção de surfactante
    def loop(self, d_idx, d_cs, d_rho_b_grown, d_a_c_s):
        dcs_dt = self.sigma * d_rho_b_grown[d_idx] - self.lambda_ * d_cs[d_idx]
        
        # Armazena a taxa de mudança de c_s na sua propriedade de taxa ('a_c_s')
        d_a_c_s[d_idx] = dcs_dt

class SurfactantForceAndDrag(Equation):
    def __init__(self, dest, sources, beta, gamma, **kw):
        super(SurfactantForceAndDrag, self).__init__(dest, sources, **kw)
        self.beta = beta
        self.gamma = gamma

    # Esta equação afeta as acelerações padrão 'ax' e 'ay'
    def loop(self, d_idx, d_u, d_v, d_ax, d_ay):
        d_ax[d_idx] += -self.gamma * d_u[d_idx]
        d_ay[d_idx] += -self.gamma * d_v[d_idx]