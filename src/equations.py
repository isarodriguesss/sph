from pysph.sph.equation import Equation

class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, **kw):
        super(BiomassGrowth, self).__init__(dest, sources, **kw)
        self.r_growth = r_growth
        self.rho_max = rho_max

    def initialize(self, d_idx, d_rho, d_rho_b_grown):
        d_rho_b_grown[d_idx] = d_rho[d_idx]

    def post_step(self, d_idx, d_rho, d_rho_b_grown, dt):
        d_rho_dt = self.r_growth * d_rho_b_grown[d_idx] * (1 - d_rho_b_grown[d_idx] / self.rho_max)
        d_rho_b_grown[d_idx] += d_rho_dt * dt
        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], self.rho_max))

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