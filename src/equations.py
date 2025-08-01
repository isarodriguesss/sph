from pysph.sph.equation import Equation
import numpy as np

class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, **kw):
        super(BiomassGrowth, self).__init__(dest, sources, **kw)
        self.r_growth = r_growth
        self.rho_max = rho_max

    def loop(self, d_idx, d_rho_b_grown):
        current_rho_b = d_rho_b_grown[d_idx]
        rate_of_change = self.r_growth * current_rho_b * (1 - current_rho_b / self.rho_max)

        if d_idx == 0:
            print(f"DEBUG_GROWTH_LOOP: r_growth={self.r_growth:.4f}, rho_max={self.rho_max:.4f}")
            print(f"DEBUG_GROWTH_LOOP: Partícula {d_idx}: rho_b={current_rho_b:.4f}, Taxa_de_mudanca={rate_of_change:.6e}")
            if abs(rate_of_change) < 1e-12:
                print(f"DEBUG_GROWTH_LOOP: Taxa de mudança muito pequena! Termo (1 - rho_b/rho_max) = {(1 - current_rho_b / self.rho_max):.6e}")

        return rate_of_change

    def post_loop(self, d_idx, d_rho_b_grown, d_rho, dt):
        current_val = d_rho_b_grown[d_idx]
        clipped_val = max(0.0, min(current_val, self.rho_max))

        if d_idx == 0:
            if abs(clipped_val - current_val) > 1e-9:
                print(f"DEBUG_POST_LOOP: Partícula {d_idx}: Valor antes da clipagem={current_val:.4f}, Valor clipado={clipped_val:.4f}")
            else:
                print(f"DEBUG_POST_LOOP: Partícula {d_idx}: Valor após cálculo={current_val:.4f}, Não clipado.")

        d_rho_b_grown[d_idx] = clipped_val

class SurfactantProductionDecay(Equation):
    def __init__(self, dest, sources, sigma, lambda_, **kw):
        super(SurfactantProductionDecay, self).__init__(dest, sources, **kw)
        self.sigma = sigma
        self.lambda_ = lambda_

    def initialize(self, d_idx, d_cs):
        if d_cs[d_idx] <= 1e-9:
            d_cs[d_idx] = 1e-9

    def post_step(self, d_idx, d_cs, d_rho_b_grown, dt):
        dcs_dt = self.sigma * d_rho_b_grown[d_idx] - self.lambda_ * d_cs[d_idx]
        d_cs[d_idx] += dcs_dt * dt
        d_cs[d_idx] = max(d_cs[d_idx], 1e-9)

class SurfactantForceAndDrag(Equation):
    def __init__(self, dest, sources, beta, gamma, **kw):
        super(SurfactantForceAndDrag, self).__init__(dest, sources, **kw)
        self.beta = beta
        self.gamma = gamma

    def post_step(self, d_idx, d_u, d_v, d_ax, d_ay, dt):
        d_ax[d_idx] += -self.gamma * d_u[d_idx]
        d_ay[d_idx] += -self.gamma * d_v[d_idx]