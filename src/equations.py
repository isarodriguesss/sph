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


class SurfactantEquation(Equation):
    def __init__(self, dest, sources, sigma, lambda_, D):
        self.D = D
        self.sigma = sigma
        self.lambda_ = lambda_
        super(SurfactantEquation, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_s):
        d_a_c_s[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_rho, d_cs, s_cs, s_m, RIJ, XIJ, DWIJ, d_a_c_s, d_h):
        cs_i = d_cs[d_idx]
        cs_j = s_cs[s_idx]
        rho_j = s_rho[s_idx]

        cs_ij = cs_i - cs_j

        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2

        xij_dot_dwij = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        term = (s_m[s_idx] / rho_j) * (cs_ij / rij_sq) * xij_dot_dwij

        d_a_c_s[d_idx] += 2.0 * self.D * term

    def post_loop(self, d_idx, d_rho_b_grown, d_a_c_s, d_cs):
        diffusion_rate = d_a_c_s[d_idx]

        prodution_rate = self.sigma * d_rho_b_grown[d_idx]
        decay_rate = self.lambda_ * d_cs[d_idx]

        d_a_c_s[d_idx] = diffusion_rate + prodution_rate - decay_rate


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


class ViscousForce(Equation):
    """
    Calcula a força de viscosidade como mu * Laplaciano(v).
    Usa a mesma forma estável do Laplaciano para evitar instabilidades.
    """

    def __init__(self, dest, sources, mu):
        self.mu = mu
        super(ViscousForce, self).__init__(dest, sources)

    def loop(
        self,
        d_idx,
        s_idx,
        d_rho,
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

        # Diferença de velocidade
        u_ij = d_u[d_idx] - s_u[s_idx]
        v_ij = d_v[d_idx] - s_v[s_idx]

        # Termo de suavização para estabilidade
        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2

        # Produto escalar (∇W ⋅ r)
        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        # Aceleração devido à viscosidade
        # O termo completo é (2 * mu / rho_i) * sum(...)
        # O PySPH já divide pela densidade no integrador, então aqui calculamos a força por unidade de volume.
        # Mas para ser consistente, adicionamos à aceleração. A densidade rho_i será considerada.
        common_term = (s_m[s_idx] / rho_j) * (dot_product / rij_sq)

        # Fator 2*mu vem da formulação do Laplaciano para vetores
        acc_x = 2.0 * self.mu * common_term * u_ij
        acc_y = 2.0 * self.mu * common_term * v_ij

        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y
