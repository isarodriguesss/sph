from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max):
        self.r_growth = r_growth
        self.rho_max = rho_max
        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho_b_grown, d_a_rho_b_grown):
        if d_rho_b_grown[d_idx] > 1e-12:
            rate = (
                self.r_growth
                * d_rho_b_grown[d_idx]
                * (1.0 - d_rho_b_grown[d_idx] / self.rho_max)
            )
            d_a_rho_b_grown[d_idx] = rate
        else:
            d_a_rho_b_grown[d_idx] = 0.0


class SurfactantProductionDecay(Equation):
    def __init__(self, dest, sources, sigma, lambda_):
        self.sigma = sigma
        self.lambda_ = lambda_
        super(SurfactantProductionDecay, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho_b_grown, d_cs, d_a_c_s):
        rate = self.sigma * d_rho_b_grown[d_idx] - self.lambda_ * d_cs[d_idx]
        d_a_c_s[d_idx] += rate


class SurfactantDiffusion(Equation):
    def __init__(self, dest, sources, D):
        self.D = D
        super(SurfactantDiffusion, self).__init__(dest, sources)

    def loop(self, d_idx, s_idx, s_m, s_rho, d_cs, d_rho, s_cs, d_a_c_s, RIJ, DWIJ):
        cs_ij = d_cs[d_idx] - s_cs[s_idx]
        rho_ij = (d_rho[d_idx] + s_rho[s_idx]) * 0.5  

        dw_x = DWIJ[0]
        dw_y = DWIJ[1]

        term = (s_m[s_idx] / s_rho[s_idx]) * (cs_ij / RIJ) * (dw_x + dw_y)

        d_a_c_s[d_idx] += 2.0 * self.D * term


class MarangoniForce(Equation):
    def __init__(self, dest, sources, beta):
        self.beta = -beta
        super(MarangoniForce, self).__init__(dest, sources)

    def loop(self, d_idx, s_idx, s_m, d_rho, s_rho, d_cs, s_cs, d_au, d_av, DWIJ):
        cs_i = d_cs[d_idx]
        cs_j = s_cs[s_idx]
        rho_i = d_rho[d_idx]
        rho_j = s_rho[s_idx]

        factor = (cs_i / (rho_i**2)) + (cs_j / (rho_j**2))

        d_au[d_idx] += self.beta * s_m[s_idx] * factor * DWIJ[0]
        d_av[d_idx] += self.beta * s_m[s_idx] * factor * DWIJ[1]


class LinearDrag(Equation):
    def __init__(self, dest, sources, gamma):
        self.gamma = -gamma
        super(LinearDrag, self).__init__(dest, sources)

    def loop(self, d_idx, d_u, d_v, d_au, d_av):
        d_au[d_idx] += self.gamma * d_u[d_idx]
        d_av[d_idx] += self.gamma * d_v[d_idx]


class InterpolateVelocity(Equation):
    def initialize(self, d_idx, d_u, d_v):
        d_u[d_idx] = 0.0
        d_v[d_idx] = 0.0

    def loop(self, d_idx, s_idx, d_u, d_v, s_u, s_v, s_m, s_rho, WIJ):
        vol_j = s_m[s_idx] / s_rho[s_idx]

        wij = WIJ

        d_u[d_idx] += s_u[s_idx] * vol_j * wij
        d_v[d_idx] += s_v[s_idx] * vol_j * wij
