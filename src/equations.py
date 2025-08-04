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

    def loop(self, d_idx, s_idx, d_rho, s_rho, d_cs, s_cs, s_m, d_a_c_s, RIJ, DWIJ):
        cs_i = d_cs[d_idx]
        cs_j = s_cs[s_idx]
        rho_j = s_rho[s_idx]

        factor = 2.0 * s_m[s_idx] / rho_j

        term = (
            (cs_i - cs_j)
            / (RIJ * RIJ + 1e-6 * self.D)
            * (DWIJ[0] * RIJ + DWIJ[1] * RIJ)
        )

        d_a_c_s[d_idx] += self.D * factor * term


class MarangoniForce(Equation):
    def __init__(self, dest, sources, beta):
        self.beta = -beta
        super(MarangoniForce, self).__init__(dest, sources)

    def loop(self, d_idx, s_idx, d_rho, s_rho, d_cs, s_cs, s_m, d_au, d_av, DWIJ):
        cs_i = d_cs[d_idx]
        cs_j = s_cs[s_idx]

        grad_cs_x = (cs_j - cs_i) * DWIJ[0]
        grad_cs_y = (cs_j - cs_i) * DWIJ[1]

        factor = s_m[s_idx] / (d_rho[d_idx] * s_rho[s_idx])

        d_au[d_idx] += self.beta * factor * grad_cs_x
        d_av[d_idx] += self.beta * factor * grad_cs_y


class LinearDrag(Equation):
    def __init__(self, dest, sources, gamma):
        self.gamma = -gamma
        super(LinearDrag, self).__init__(dest, sources)

    def loop(self, d_idx, d_u, d_v, d_au, d_av):
        d_au[d_idx] += self.gamma * d_u[d_idx]
        d_av[d_idx] += self.gamma * d_v[d_idx]
