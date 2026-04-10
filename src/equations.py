from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, density_limit=1.1):
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.density_limit = density_limit

        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho_b_grown, d_a_rho_b_grown, d_m, d_am, d_rho):
        d_a_rho_b_grown[d_idx] = 0.0
        d_am[d_idx] = 0.0
        if d_rho_b_grown[d_idx] > 1e-12:
            rate = self.r_growth * (1.0 - d_rho_b_grown[d_idx] / self.rho_max)
            d_a_rho_b_grown[d_idx] = rate * d_rho_b_grown[d_idx]
            d_am[d_idx] = rate * d_m[d_idx]


class BiomassGradient(Equation):
    def initialize(self, d_idx, d_grad_rho_b_x, d_grad_rho_b_y):
        d_grad_rho_b_x[d_idx] = 0.0
        d_grad_rho_b_y[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        s_rho,
        d_rho_b_grown,
        s_rho_b_grown,
        DWIJ,
        d_grad_rho_b_x,
        d_grad_rho_b_y,
    ):
        vol_j = s_m[s_idx] / s_rho[s_idx]

        diff = s_rho_b_grown[s_idx] - d_rho_b_grown[d_idx]

        d_grad_rho_b_x[d_idx] += vol_j * diff * DWIJ[0]
        d_grad_rho_b_y[d_idx] += vol_j * diff * DWIJ[1]

    def post_loop(self, d_idx, d_grad_rho_b_x, d_grad_rho_b_y, d_grad_rho_b_mag):
        gx = d_grad_rho_b_x[d_idx]
        gy = d_grad_rho_b_y[d_idx]

        d_grad_rho_b_mag[d_idx] = (gx * gx + gy * gy) ** 0.5


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

    def post_loop(self, d_idx, d_rho_b_grown, d_a_c_s, d_cs, d_noise, d_grad_rho_b_mag):
        grad = min(d_grad_rho_b_mag[d_idx], 10.0)

        production = self.sigma * d_rho_b_grown[d_idx] * (0.3 + grad) * d_noise[d_idx]

        # Lambda constante: o lambda espacialmente dependente foi removido porque ele matava
        # o cs na interface (onde rho_b < 0.1), destruindo o gradiente antes que o SPH
        # pudesse calculá-lo. O confinamento do surfactante é garantido pelo D pequeno.
        d_a_c_s[d_idx] += production - self.lambda_ * d_cs[d_idx]


class MarangoniForce(Equation):
    """
    Força de Marangoni aplicada APENAS na interface da colônia.

    Detecção de interface: usa d_grad_rho_b_mag (calculado por BiomassGradient).
    No bulk interior grad_rho_b ≈ 0 → força zero. No border do dendrito
    grad_rho_b é grande → força máxima. Isso previne empuxo espúrio por
    ruído de cs no interior.

    gate = smoothstep(grad_rho_b_mag, grad_low, grad_high)
    F = gate · (-β) · ∇cs
    """

    def __init__(self, dest, sources, beta, grad_low=0.1, grad_high=1.0):
        self.beta = -beta
        self.grad_low = grad_low
        self.grad_high = grad_high
        super(MarangoniForce, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au_mar, d_ax_mar, d_ay_mar):
        d_au_mar[d_idx] = 0.0
        d_ax_mar[d_idx] = 0.0
        d_ay_mar[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        d_rho,
        s_rho,
        d_cs,
        s_cs,
        d_au,
        d_av,
        d_ax_mar,
        d_ay_mar,
        d_grad_rho_b_mag,
        DWIJ,
    ):
        # Gate de interface: só atua onde |∇rho_b| é significativo
        grad_mag = d_grad_rho_b_mag[d_idx]
        if grad_mag >= self.grad_low:
            gate_raw = (grad_mag - self.grad_low) / (self.grad_high - self.grad_low)
            if gate_raw > 1.0:
                gate = 1.0
            else:
                gate = gate_raw * gate_raw * (3.0 - 2.0 * gate_raw)

            vol_j = s_m[s_idx] / s_rho[s_idx]
            cs_ij = s_cs[s_idx] - d_cs[d_idx]

            acc_x = gate * self.beta * vol_j * cs_ij * DWIJ[0]
            acc_y = gate * self.beta * vol_j * cs_ij * DWIJ[1]

            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_ax_mar[d_idx] += acc_x
            d_ay_mar[d_idx] += acc_y

    def post_loop(self, d_idx, d_au_mar, d_ax_mar, d_ay_mar):
        # Registra a magnitude da aceleração vetorial LÍQUIDA de Marangoni
        # (após cancelamento de contribuições de vizinhos).
        d_au_mar[d_idx] = (
            d_ax_mar[d_idx] * d_ax_mar[d_idx] + d_ay_mar[d_idx] * d_ay_mar[d_idx]
        ) ** 0.5


class LinearDrag(Equation):
    """
    Drag linear com GRADIENTE DE MOBILIDADE biológico.

    Modelo: γ_eff = γ_base + γ_mature · rho_b²

    - Partículas com rho_b ≈ 0 (fluido livre): γ_eff = γ_base (mobilidade total)
    - Partículas com rho_b ≈ 0.3 (interface/edge swarm): γ_eff ≈ γ_base + 0.09·γ_mature
    - Partículas com rho_b ≈ 1.0 (núcleo maduro / biofilme): γ_eff = γ_base + γ_mature

    Justificativa biológica: em P. aeruginosa swarming, células no núcleo
    estão embebidas em matriz EPS (biofilme) e são essencialmente imóveis,
    enquanto células no edge swarm são motoras (flagelo + surfactante).
    Sem esse gradiente, o transiente inicial varre todo o material para
    fora junto com a interface, criando o anel oco que vimos.

    Use γ_mature ~5–10× γ_base para diferenciar fortemente os regimes.
    """

    def __init__(self, dest, sources, gamma_base, gamma_mature=0.0):
        self.gamma_base = -gamma_base
        self.gamma_mature = -gamma_mature
        super(LinearDrag, self).__init__(dest, sources)

    def loop(self, d_idx, d_u, d_v, d_au, d_av, d_au_drag, d_rho_b_grown):
        rho_b = d_rho_b_grown[d_idx]
        # Quadrático em rho_b: edge mobile, núcleo congelado
        gamma_eff = self.gamma_base + self.gamma_mature * rho_b * rho_b
        acc_x = gamma_eff * d_u[d_idx]
        acc_y = gamma_eff * d_v[d_idx]
        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y
        d_au_drag[d_idx] = (acc_x * acc_x + acc_y * acc_y) ** 0.5


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

        u_ij = d_u[d_idx] - s_u[s_idx]
        v_ij = d_v[d_idx] - s_v[s_idx]

        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2

        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        common_term = (s_m[s_idx] / rho_j) * (dot_product / rij_sq)

        # Viscosidade ponderada pela densidade: interface (rho < rho0) é menos viscosa,
        # permitindo que as pontas dos dendritos se deformem livremente. O interior denso
        # (rho ~ rho0) mantém a coesão da colônia com viscosidade plena.
        rho_avg = 0.5 * (d_rho[d_idx] + rho_j)
        mu_eff = self.mu * min(rho_avg, 1.0)

        acc_x = 2.0 * mu_eff * common_term * u_ij
        acc_y = 2.0 * mu_eff * common_term * v_ij

        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y


class BiomassEOS(Equation):
    """
    EOS "Soft Interior, Cohesive Edge" com tensão superficial intrínseca.

    Dois regimes:
    - Compressão (ρ > ρ₀): pressão repulsiva quadrática (evita overlap)
    - Rarefação (ρ < ρ₀): pressão levemente NEGATIVA (tensão superficial)

    Via MomentumEquation do PySPH, p < 0 cria força atrativa entre vizinhos
    → substitui a OsmoticForce como mecanismo de coesão. O regime atrativo
    é modulado pelo edge_fade: zero na borda livre (rho_b < 0.1), máximo
    no interior (rho_b > 0.5).

    P = B · excess² · edge_fade          se ρ > ρ₀
    P = -B_tension · deficit · edge_fade  se ρ < ρ₀

    - B = ρ₀·c₀²/γ — repulsão forte para impedir colapso
    - B_tension = B · tension_ratio — atração fraca para coesão (ratio ~0.3)
    - edge_fade: smoothstep em rho_b ∈ [0.1, 0.5]
    """

    def __init__(self, dest, sources, rho0, c0, gamma_eos=7.0, tension_ratio=0.1):
        self.rho0 = rho0
        self.c0 = c0
        self.B = rho0 * c0 * c0 / gamma_eos
        self.B_tension = self.B * tension_ratio
        super(BiomassEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_rho_b_grown):
        ratio = d_rho[d_idx] / self.rho0
        rho_b = d_rho_b_grown[d_idx]

        # Edge fade: smoothstep — pressão zero onde rho_b < 0.1,
        # transição suave, full onde rho_b > 0.5. Garante borda livre
        # (sem pressão) para as pontas dos dendritos.
        if rho_b < 0.1:
            edge_fade = 0.0
        elif rho_b < 0.5:
            t = (rho_b - 0.1) / 0.4
            edge_fade = t * t * (3.0 - 2.0 * t)
        else:
            edge_fade = 1.0

        if ratio > 1.0:
            # Compressão: repulsão quadrática suave
            excess = ratio - 1.0
            d_p[d_idx] = self.B * excess * excess * edge_fade
        else:
            # Rarefação: tensão superficial (atração leve)
            # Limitamos o deficit a 0.3 para evitar atração excessiva em gaps
            deficit = 1.0 - ratio
            if deficit > 0.3:
                deficit = 0.3
            d_p[d_idx] = -self.B_tension * deficit * edge_fade


class OsmoticForce(Equation):
    """
    Força osmótica anisotrópica: preenche buracos no interior sem circularizar.

    Princípio: A força aponta na direção ∇(rho_SPH) — detecta rarefação
    mecânica (gaps no fluido). Modulada por quão "interior" a partícula é.
    Na borda da colônia (rho_b baixo) a força é ZERO.

    F = k_osm · interior_gate(rho_b) · ∇(rho_SPH)

    - interior_gate: smoothstep 0→1 para rho_b ∈ [0.3, 0.7]
    - ∇(rho_SPH) calculado par-a-par no loop SPH
    - k_osm pequeno (ratio Marangoni/Osmótica ≈ 20:1)
    """

    def __init__(self, dest, sources, k_osm=0.5, rho_b_threshold=0.3):
        self.k_osm = k_osm
        self.rho_b_threshold = rho_b_threshold
        super(OsmoticForce, self).__init__(dest, sources)

    def initialize(self, d_idx, d_au_osm):
        d_au_osm[d_idx] = 0.0

    def loop(
        self, d_idx, s_idx, s_m, d_rho, s_rho, d_rho_b_grown, d_au, d_av, d_au_osm, DWIJ
    ):
        rho_b = d_rho_b_grown[d_idx]

        # Gate: só no interior (rho_b > 0.3), smoothstep até 0.7
        if rho_b >= self.rho_b_threshold:
            gate_raw = (rho_b - self.rho_b_threshold) / 0.4
            if gate_raw > 1.0:
                gate = 1.0
            else:
                gate = gate_raw * gate_raw * (3.0 - 2.0 * gate_raw)

            # Gradiente da densidade SPH — detecta gaps no fluido
            vol_j = s_m[s_idx] / s_rho[s_idx]
            rho_diff = s_rho[s_idx] - d_rho[d_idx]

            acc_x = self.k_osm * gate * vol_j * rho_diff * DWIJ[0]
            acc_y = self.k_osm * gate * vol_j * rho_diff * DWIJ[1]

            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_au_osm[d_idx] += (acc_x * acc_x + acc_y * acc_y) ** 0.5
