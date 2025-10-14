from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max):
        self.r_growth = r_growth
        self.rho_max = rho_max

        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho_b_grown, d_a_rho_b_grown, d_m, d_am):
        rate = 0.0
        if d_rho_b_grown[d_idx] > 1e-12:
            rate = (
                self.r_growth
                * d_rho_b_grown[d_idx]
                * (1.0 - d_rho_b_grown[d_idx] / self.rho_max)
            )
            d_a_rho_b_grown[d_idx] = rate

        d_am[d_idx] = rate * d_m[d_idx]


class MassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max):
        self.r_growth = r_growth
        self.rho_max = rho_max
        super(MassGrowth, self).__init__(dest, sources)

    def loop(self, d_idx, d_m, d_rho_b_grown, d_am):
        if d_rho_b_grown[d_idx] > 1e-12:
            growth_rate_factor = self.r_growth * (
                1.0 - d_rho_b_grown[d_idx] / self.rho_max
            )
            dm_dt = growth_rate_factor * d_m[d_idx]
            d_am[d_idx] = dm_dt


class SurfactantEquation(Equation):
    def __init__(self, dest, sources, sigma, lambda_, D):
        self.D = D
        self.sigma = sigma
        self.lambda_ = lambda_
        super(SurfactantEquation, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_s):
        d_a_c_s[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_rho, d_cs, s_cs, s_m, RIJ, XIJ, DWIJ, d_a_c_s, d_h):
        cs_ij = d_cs[d_idx] - s_cs[s_idx]
        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2
        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        term = (s_m[s_idx] / s_rho[s_idx]) * (cs_ij / rij_sq) * dot_product
        d_a_c_s[d_idx] += 2.0 * self.D * term

    def post_loop(self, d_idx, d_rho_b_grown, d_a_c_s, d_cs):
        reaction_rate = self.sigma * d_rho_b_grown[d_idx] - self.lambda_ * d_cs[d_idx]

        d_a_c_s[d_idx] += reaction_rate


class MarangoniForce(Equation):
    def __init__(self, dest, sources, beta, acc_limit=10.0):
        self.beta = -beta
        self.acc_limit = acc_limit
        self.acc_limit_sq = acc_limit * acc_limit
        super(MarangoniForce, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au, d_av):
        d_au[d_idx] = 0.0
        d_av[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, d_rho, s_rho, d_cs, s_cs, d_au, d_av, DWIJ):
        vol_j = s_m[s_idx] / s_rho[s_idx]
        cs_ij = s_cs[s_idx] - d_cs[d_idx]

        acc_x = self.beta * vol_j * cs_ij * DWIJ[0]
        acc_y = self.beta * vol_j * cs_ij * DWIJ[1]

        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y

    def post_loop(self, d_idx, d_au, d_av):
        acc_sq = d_au[d_idx] ** 2 + d_av[d_idx] ** 2

        if acc_sq > self.acc_limit_sq:
            scale = self.acc_limit / (acc_sq**0.5)
            d_au[d_idx] *= scale
            d_av[d_idx] *= scale


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


class ViscousForce(Equation):
    def __init__(self, dest, sources, mu):
        self.mu = mu
        super(ViscousForce, self).__init__(dest, sources)

    def loop(
        self,
        d_idx,
        s_idx,
        s_rho,
        d_u,
        d_v,
        s_u,
        s_v,
        s_m,
        RIJ,
        XIJ,
        DWIJ,
        d_au,
        d_av,
        d_h,
    ):
        rho_j = s_rho[s_idx]

        u_ij = d_u[d_idx] - s_u[s_idx]
        v_ij = d_v[d_idx] - s_v[s_idx]

        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2

        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        common_term = (s_m[s_idx] / rho_j) * (dot_product / rij_sq)

        acc_x = 2.0 * self.mu * common_term * u_ij
        acc_y = 2.0 * self.mu * common_term * v_ij

        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y
