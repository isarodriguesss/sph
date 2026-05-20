from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, density_limit=1.1):
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.density_limit = density_limit

        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(
        self, d_idx, d_rho_b_grown, d_a_rho_b_grown, d_m, d_am, d_rho, d_c_n, d_u, d_v
    ):
        d_a_rho_b_grown[d_idx] = 0.0
        d_am[d_idx] = 0.0
        if d_rho_b_grown[d_idx] > 1e-12 and d_rho_b_grown[d_idx] < 0.8:
            c_n = d_c_n[d_idx]
            if c_n < 0.6:
                c_n_factor = 0.0
            elif c_n > 0.9:
                c_n_factor = 1.0
            else:
                t = (c_n - 0.6) / 0.3
                c_n_factor = t * t * (3.0 - 2.0 * t)

            # Gate por motilidade: so swarmers ativos (|v| ~ v_term=f0/gamma=0.05)
            # crescem rho_b. Rim parado da casca (|v|~0.001 por Marangoni anulado
            # radialmente) nao cresce. Discrimina tip avancando em agar virgem
            # vs rim mature em torno do nucleo — ambos tem c_n~0.9, mas so o tip
            # se move. Smoothstep [0.01, 0.05] = [v_pinned, v_term].
            v_mag = (d_u[d_idx] * d_u[d_idx] + d_v[d_idx] * d_v[d_idx]) ** 0.5
            if v_mag < 0.01:
                motility_gate = 0.0
            elif v_mag > 0.05:
                motility_gate = 1.0
            else:
                tv = (v_mag - 0.01) / 0.04
                motility_gate = tv * tv * (3.0 - 2.0 * tv)

            rate = (
                self.r_growth * (1.0 - d_rho_b_grown[d_idx] / self.rho_max) * c_n_factor
                # * motility_gate
            )
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
    def __init__(
        self,
        dest,
        sources,
        sigma,
        lambda_,
        D,
        D_ext=0.03,
        lambda_ext_ratio=5.0,
        k_consume=0.0,
        cs_max=0.5,
    ):
        self.D = D
        self.D_ext = D_ext
        self.sigma = sigma
        self.lambda_ = lambda_
        self.lambda_ext = lambda_ * lambda_ext_ratio
        self.k_consume = k_consume
        self.cs_max = cs_max
        super(SurfactantEquation, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_s):
        d_a_c_s[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        s_rho,
        d_cs,
        s_cs,
        s_m,
        RIJ,
        XIJ,
        DWIJ,
        d_a_c_s,
        d_h,
        d_rho_b_grown,
        s_rho_b_grown,
    ):
        # Difusao bi-escala: D_int dentro da colonia, D_ext no agar exterior.
        # Ramnolipideo difunde rapido no agar livre e lento dentro do biofilme (EPS).
        # Usa a media geometrica dos rho_b do par para transicao suave.
        rho_b_avg = 0.5 * (d_rho_b_grown[d_idx] + s_rho_b_grown[s_idx])
        # Smoothstep gate: D_ext para rho_b<0.1, D_int para rho_b>0.5
        if rho_b_avg < 0.1:
            D_eff = self.D_ext
        elif rho_b_avg < 0.5:
            t = (rho_b_avg - 0.1) / 0.4
            gate = t * t * (3.0 - 2.0 * t)
            D_eff = self.D_ext + (self.D - self.D_ext) * gate
        else:
            D_eff = self.D

        cs_ij = d_cs[d_idx] - s_cs[s_idx]
        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2
        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        term = (s_m[s_idx] / s_rho[s_idx]) * (cs_ij / rij_sq) * dot_product
        d_a_c_s[d_idx] += 2.0 * D_eff * term

    def post_loop(
        self,
        d_idx,
        d_rho_b_grown,
        d_a_c_s,
        d_cs,
        d_noise,
        d_c_n,
    ):
        rho_b = d_rho_b_grown[d_idx]
        qs = rho_b * rho_b / (rho_b * rho_b + 0.01)

        # Pass T1 (Trinschek-like, §3.0 [T1] painel (b)): producao satura
        # localmente em cs_max via (1 - cs/cs_max). Substitui a combinacao
        # growth_headroom + tip_boost + motile_boost + c_n_factor + k_consume
        # por uma unica forma analoga a Trinschek 2018 eq. para Γ.
        # Resultado esperado (validado §3.3.6): cs uniforme = cs_max em rho_b>0.5,
        # decai monotonicamente no rim+agar — perfil §3.3-compliant.
        saturation = 1.0 - d_cs[d_idx] / self.cs_max
        if saturation < 0.0:
            saturation = 0.0

        c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + 0.1)

        # Pontas rápidas produzem mais surfactante
        # Problema: os dendritos correm mais rápido do que a química conseguiu acompanhar.
        # Splitting resolveria o problema
        # v_mag = sqrt(d_u[d_idx] * d_u[d_idx] + d_v[d_idx] * d_v[d_idx])
        # if v_mag > 0.1:
        #     v_mag = 0.1
        # motile_boost = 1.0 + 100.0 * v_mag  # 1x em estatico, 3.5x em v=0.1

        # revisar: qual diferença com e sem c_n_factor?
        production = self.sigma * qs * saturation * d_noise[d_idx] * c_n_factor

        # Decaimento somente no agar (rho_b < 0.1). No biofilme a saturacao
        # via (1-cs/cs_max) ja limita cs em cs_max — nao precisa de sumidouro.
        # Isto reproduz a forma de Γ uniforme em Trinschek (b): biomassa cheia,
        # halo no agar com decaimento finito.
        if rho_b < 0.1:
            lambda_eff = self.lambda_ * 2.0  # Decaimento LENTO no ágar (gera o halo)
        else:
            lambda_eff = self.lambda_ * 1.0

        d_a_c_s[d_idx] += production - lambda_eff * d_cs[d_idx]


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

    def __init__(self, dest, sources, beta, grad_low=0.15, grad_high=0.5):
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

    def loop(
        self,
        d_idx,
        d_u,
        d_v,
        d_au,
        d_av,
        d_au_drag,
        d_rho_b_grown,
        d_ax_drag,
        d_ay_drag,
    ):
        rho_b = d_rho_b_grown[d_idx]
        # Quadrático em rho_b: edge mobile, núcleo congelado
        gamma_eff = self.gamma_base + self.gamma_mature * rho_b * rho_b
        acc_x = gamma_eff * d_u[d_idx]
        acc_y = gamma_eff * d_v[d_idx]
        d_ax_drag[d_idx] = acc_x
        d_ay_drag[d_idx] = acc_y
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

    def __init__(self, dest, sources, rho0, c0, gamma_eos=7.0, tension_ratio=0.02):
        self.rho0 = rho0
        self.c0 = c0
        self.B = rho0 * c0 * c0 / gamma_eos
        self.B_tension = self.B * tension_ratio
        super(BiomassEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_rho_b_grown):
        ratio = d_rho[d_idx] / self.rho0
        rho_b = d_rho_b_grown[d_idx]

        # K.16e: edge_fade ASSIMETRICO — repulsao e atracao desacopladas.
        # Repulsao (compressao): zerada no agar (rho_b<0.1) E no nucleo (rho_b>=0.8)
        #   → elimina a_pressure no core saturado (alavanca §2.4).
        # Atracao (rarefacao): zerada so no agar (rho_b<0.1), preservada no core
        #   → mantem coesao estrutural do nucleo EPS, evita colapso K.16c.
        if rho_b < 0.1:
            fade_rep = 0.0
            fade_att = 0.0
        elif rho_b < 0.5:
            t = (rho_b - 0.1) / 0.4
            s = t * t * (3.0 - 2.0 * t)
            fade_rep = s
            fade_att = s
        elif rho_b < 0.8:
            t = (rho_b - 0.5) / 0.3
            s = t * t * (3.0 - 2.0 * t)
            fade_rep = 1.0 - s
            fade_att = 1.0
        else:
            fade_rep = 0.0
            fade_att = 1.0

        if ratio > 1.0:
            # Compressão: repulsão quadrática suave (zerada no core)
            excess = ratio - 1.0
            d_p[d_idx] = self.B * excess * excess * fade_rep
        else:
            # Rarefação: tensão superficial (atração leve, preservada no core)
            # Limitamos o deficit a 0.3 para evitar atração excessiva em gaps
            deficit = 1.0 - ratio
            if deficit > 0.3:
                deficit = 0.3
            d_p[d_idx] = -self.B_tension * deficit * fade_att


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


class OsmolyteProduction(Equation):
    """
    Pass M-A — Campo de osmolito c_o produzido pelas bacterias (Frente 6).

    Mecanismo (Srinivasan 2019, Bru 2023): bacterias secretam LPS/EPS que
    funcionam como osmolitos. c_o satura nas baias (agar confinado entre
    dois dendritos) e fica ~0 nas pontas (agar virgem). |∇c_o| pequeno
    nas baias → influxo zero. |∇c_o| grande nas pontas → influxo de massa
    via van't Hoff (V0 ∝ |∇c_o|). Resultado: pontas incham, baias estagnam.

    Loop: difusao SPH (Brookshaw) + acumulo de gradiente simetrico.
    Post_loop: producao biomassa-dependente, decaimento, e influxo de massa
    proporcional a |∇c_o| · gate(rho_b).
    """

    def __init__(self, dest, sources, D_o=1e-3, k_o=0.5, lambda_o=0.05, Q0=5e-4):
        self.D_o = D_o
        self.k_o = k_o
        self.lambda_o = lambda_o
        self.Q0 = Q0
        super(OsmolyteProduction, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_o, d_grad_co_x, d_grad_co_y):
        d_a_c_o[d_idx] = 0.0
        d_grad_co_x[d_idx] = 0.0
        d_grad_co_y[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        d_c_o,
        s_c_o,
        d_a_c_o,
        d_grad_co_x,
        d_grad_co_y,
        s_rho,
        s_m,
        DWIJ,
        XIJ,
        RIJ,
    ):
        co_ij = s_c_o[s_idx] - d_c_o[d_idx]
        Vj = s_m[s_idx] / s_rho[s_idx]

        d_grad_co_x[d_idx] += Vj * co_ij * DWIJ[0]
        d_grad_co_y[d_idx] += Vj * co_ij * DWIJ[1]

        if RIJ > 1e-12:
            eij_dot_dwij = (XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]) / (RIJ * RIJ)
            d_a_c_o[d_idx] += 2.0 * self.D_o * Vj * co_ij * eij_dot_dwij

    def post_loop(
        self,
        d_idx,
        d_a_c_o,
        d_c_o,
        d_rho_b_grown,
        d_grad_co_x,
        d_grad_co_y,
        d_au,
        d_av,
    ):
        rho_b = d_rho_b_grown[d_idx]

        # Producao: bacterias secretam osmolito (proporcional a densidade)
        # qs = Hill QS, satura em c_o=1 via (1 - c_o)
        qs = rho_b * rho_b / (rho_b * rho_b + 0.01)
        production = self.k_o * qs * (1.0 - d_c_o[d_idx])

        # Decaimento lento (osmolitos persistem mais que cs)
        decay = self.lambda_o * d_c_o[d_idx]

        d_a_c_o[d_idx] += production - decay

        # Forca osmotica de Darcy: F = -Q0 * gate * nabla(c_o)
        # c_o alto dentro, baixo fora -> nabla(c_o) aponta para dentro.
        # -nabla(c_o) aponta para fora (expansao do rim), mesmo sinal de FlagellarForce.
        # Forma proporcional ao gradiente (nao normalizada): preserva diferenciacao
        # tip/baia — pontas tem |nabla(c_o)| maior que baias confinadas (2x).
        # Gate smoothstep em rho_b in [0.1, 0.6] com pico em rho_b = 0.35
        if rho_b < 0.1:
            gate = 0.0
        elif rho_b > 0.6:
            gate = 0.0
        else:
            t = (rho_b - 0.1) / 0.5
            gate = t * t * (3.0 - 2.0 * t)
            t2 = (rho_b - 0.35) / 0.25
            envelope = 1.0 - t2 * t2
            if envelope < 0.0:
                envelope = 0.0
            gate *= envelope

        d_au[d_idx] -= self.Q0 * gate * d_grad_co_x[d_idx]
        d_av[d_idx] -= self.Q0 * gate * d_grad_co_y[d_idx]


class FlagellarForce(Equation):
    """
    Pass K — Motilidade flagelar orientada por gradiente (Frente 5, CLAUDE.md).

    F_flag = -f0 * gate(rho_b) * n̂(∇cs)

    Forca propulsiva constante em magnitude (f0), alinhada a -∇cs (mesma direcao
    de Marangoni — para agar fresco). Sobrevive a saturacao do reservatorio de cs
    porque depende apenas da DIRECAO do gradiente, nao da magnitude. Modela o
    flagelo polar de P. aeruginosa como motor ativo quimiotactico.

    Gate rho_b em [0.1, 0.6] com pico em 0.35: so swarmers de borda recebem a forca,
    nucleo maduro (rho_b~1) e agar livre (rho_b~0) sao imunes.

    Loop: acumula gradiente SPH completo de cs em grad_cs_x/y.
    Post_loop: normaliza o vetor total e aplica com sinal correto.
    """

    def __init__(self, dest, sources, f0):
        self.f0 = f0
        super().__init__(dest, sources)

    def initialize(self, d_idx, d_au_flag, d_grad_cs_x, d_grad_cs_y):
        d_au_flag[d_idx] = 0.0
        d_grad_cs_x[d_idx] = 0.0
        d_grad_cs_y[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        s_rho,
        s_cs,
        d_cs,
        d_grad_cs_x,
        d_grad_cs_y,
        DWIJ,
    ):
        vol_j = s_m[s_idx] / s_rho[s_idx]
        cs_ij = s_cs[s_idx] - d_cs[d_idx]
        d_grad_cs_x[d_idx] += vol_j * cs_ij * DWIJ[0]
        d_grad_cs_y[d_idx] += vol_j * cs_ij * DWIJ[1]

    def post_loop(
        self,
        d_idx,
        d_rho_b_grown,
        d_grad_cs_x,
        d_grad_cs_y,
        d_au,
        d_av,
        d_au_flag,
    ):
        rho_b = d_rho_b_grown[d_idx]
        # Gate: swarmers na borda apenas, pico em rho_b=0.35
        if rho_b >= 0.2 and rho_b <= 0.6:
            if rho_b < 0.4:
                t = (rho_b - 0.2) / 0.2
            else:
                t = (0.6 - rho_b) / 0.2
            gate = t * t * (3.0 - 2.0 * t)

            gx = d_grad_cs_x[d_idx]
            gy = d_grad_cs_y[d_idx]
            mag = (gx * gx + gy * gy) ** 0.5 + 1e-9

            # Sinal negativo: forca aponta na direcao de -∇cs (mesma de Marangoni,
            # para o agar fresco). Magnitude constante = f0 * gate, independente
            # de |∇cs|, por isso sobrevive a saturacao do reservatorio.
            acc_x = -self.f0 * gate * gx / mag
            acc_y = -self.f0 * gate * gy / mag
            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_au_flag[d_idx] = (acc_x * acc_x + acc_y * acc_y) ** 0.5


class OxigenConsumption(Equation):
    """Pass M-B.4 (Frente 6): campo de nutriente/oxigenio consumivel c_n com difusao bi-escala.

    dc_n/dt = D_n_eff * Laplacian(c_n) - k_n * rho_b * c_n

    D_n_eff = D_n_int dentro do biofilme (EPS bloqueia transporte)
            = D_n_ext no agar livre
    Smoothstep em rho_b_avg [0.1, 0.5].

    A difusao bi-escala e o equivalente de Pass I.6 para c_n:
    - D_n_int << D_n_ext cria L_D_n_int = sqrt(D_n_int/(k_n*rho_b)) ~ 0.018
      (vs L_D_n = 0.258 com D_n uniforme).
    - Transicao entre nucleo e baia e depletada em <2s: c_n_factor ~ 0 no interior.
    - Baias geometricamente fechadas pelo biofilme nao recebem reabastecimento
      de c_n do agar exterior — producao cs e crescimento rho_b parados.
    - Tips em agar virgem (D_n_ext ativo) recebem c_n fresco — crescem normalmente.
    """

    def __init__(self, dest, sources, D_n=0.02, D_n_int=1e-4, k_n=0.5):
        self.D_n = D_n  # D_n_ext: difusao no agar livre
        self.D_n_int = D_n_int  # difusao dentro do biofilme (EPS bloqueia)
        self.k_n = k_n
        super(OxigenConsumption, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_n):
        d_a_c_n[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        d_c_n,
        s_c_n,
        d_a_c_n,
        s_rho,
        s_m,
        DWIJ,
        XIJ,
        RIJ,
        d_h,
        d_rho_b_grown,
        s_rho_b_grown,
    ):
        # M-B.6: barreira EPS [0.05, 0.3]. Em rho_b=0.3 (baia tipica) D_eff = D_n_int.
        # Transicao um pouco mais larga que M-B.5 [0.05, 0.25] — preserva blocking
        # nas baias mas suaviza descontinuidade na fronteira agar/biofilme.
        rho_b_avg = 0.5 * (d_rho_b_grown[d_idx] + s_rho_b_grown[s_idx])
        if rho_b_avg < 0.05:
            D_eff = self.D_n
        elif rho_b_avg < 0.3:
            t = (rho_b_avg - 0.05) / 0.25
            gate = t * t * (3.0 - 2.0 * t)
            D_eff = self.D_n + (self.D_n_int - self.D_n) * gate
        else:
            D_eff = self.D_n_int

        cn_ij = d_c_n[d_idx] - s_c_n[s_idx]
        Vj = s_m[s_idx] / s_rho[s_idx]
        rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2
        dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]

        d_a_c_n[d_idx] += 2.0 * D_eff * Vj * (cn_ij / rij_sq) * dot_product

    def post_loop(self, d_idx, d_a_c_n, d_c_n, d_rho_b_grown):
        # Consumo biomassa-dependente
        consumption = self.k_n * d_rho_b_grown[d_idx] * d_c_n[d_idx]
        d_a_c_n[d_idx] -= consumption
