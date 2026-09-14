from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max, density_limit=1.1):
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.density_limit = density_limit

        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(
        self,
        d_idx,
        d_rho_b_grown,
        d_a_rho_b_grown,
        d_m,
        d_am,
        d_c_n,
        d_is_filler,
        d_is_env,
    ):
        d_a_rho_b_grown[d_idx] = 0.0
        d_am[d_idx] = 0.0
        # limiar 0.01 e nao 1e-12: a banda [0, 0.01) e o PISO DE ENVELOPE, marcacao de
        # representacao e nao celula. Com o gate antigo as 7365 particulas do piso
        # entravam no crescimento (biomassa +1.8% no E6).
        if (
            d_rho_b_grown[d_idx] > 1e-12
            and d_rho_b_grown[d_idx] < 0.8
            and d_is_filler[d_idx] < 0.5
            and d_is_env[d_idx] < 0.5
        ):
            c_n = d_c_n[d_idx]
            if c_n < 0.4:
                c_n_factor = 0.0
            elif c_n > 0.8:
                c_n_factor = 1.0
            else:
                t = (c_n - 0.4) / 0.4
                c_n_factor = t * t * (3.0 - 2.0 * t)

            rate = (
                self.r_growth * (1.0 - d_rho_b_grown[d_idx] / self.rho_max) * c_n_factor
            )
            d_a_rho_b_grown[d_idx] = rate * d_rho_b_grown[d_idx]
            # licao #50: `rate*m` nao e proporcional a rho_b, entao particula de traco
            # ganha massa a taxa cheia. Inofensivo em r_growth=0.02 (por isso o C4
            # funciona), runaway assim que a taxa sobe: D3b mediu +154%.
            d_am[d_idx] = d_m[d_idx] * d_a_rho_b_grown[d_idx] / self.rho_max


class BiomassColonization(Equation):
    """Invasao de espaco vazio por biomassa vizinha (K3 — docs/PLANO_K3_JUNCAO.md).

    `BiomassGrowth` e MULTIPLICATIVO pela propria biomassa, entao rho_b = 0 e
    estado absorvente: a particula nunca se torna viva, por mais biomassa que
    exista em volta. Este e o unico termo ADITIVO. Biologicamente e a divisao
    celular — a filha ocupa o espaco vizinho.

    Relaxacao UNILATERAL em direcao a densidade do DOADOR,

        d_rho_b/dt += k_col * c_n_factor * max(0, rho_b_doador - rho_b)
        rho_b_doador = sum_j V_j rho_b_j^2 W_ij / sum_j V_j rho_b_j W_ij

    O doador e a media de rho_b PONDERADA POR rho_b, nao a media Shepard: num
    campo 92% vazio a Shepard e a media entre a mae e o vacuo (0.006-0.025 medido
    no C4) e recruta ABAIXO do quorum 0.1, onde a particula e mecanicamente
    invisivel (licao #53). A filha nasce com a densidade da MAE (licao #48).

    Quatro guardas:

    1. UNILATERAL (so soma) — o nucleo nunca drena, ao contrario da difusao (D1),
       que o esvaziou de 1.0 para 0.48 e zerou o hard pin.
    2. FILLER FORA da soma de doadores — o filler carrega rho_b herdado e
       congelado; sem esta exclusao 93% dos doadores sao fantasmas (licao #39).
    3. GATE `cs > cs_min` — recruta so onde o surfactante ja esta perto do teto,
       onde `(1 - cs/cs_max)` limita por construcao o que uma fonte nova
       acrescenta. Sem ele, 45% do recrutamento cai na frente e afoga o gradiente
       de Marangoni (licao #48).
    4. AUTO-GATEADA pela vizinhanca — numa baia as duas somas vao a zero.
    """

    def __init__(self, dest, sources, k_col, rho_max, cs_min):
        self.k_col = k_col
        self.rho_max = rho_max
        self.cs_min = cs_min
        super(BiomassColonization, self).__init__(dest, sources)

    def initialize(self, d_idx, d_rho_b_smooth, d_rho_b_w2):
        d_rho_b_smooth[d_idx] = 0.0
        d_rho_b_w2[d_idx] = 0.0

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        s_rho,
        s_rho_b_grown,
        s_is_filler,
        s_is_env,
        WIJ,
        d_rho_b_w2,
        d_rho_b_smooth,
    ):
        if s_is_filler[s_idx] < 0.5 and s_is_env[s_idx] < 0.5:
            vw = (s_m[s_idx] / s_rho[s_idx]) * s_rho_b_grown[s_idx] * WIJ
            d_rho_b_smooth[d_idx] += vw
            d_rho_b_w2[d_idx] += vw * s_rho_b_grown[s_idx]

    def post_loop(
        self,
        d_idx,
        d_rho_b_smooth,
        d_rho_b_w2,
        d_rho_b_grown,
        d_a_rho_b_grown,
        d_am,
        d_m,
        d_c_n,
        d_cs,
        d_is_filler,
    ):
        if (
            d_rho_b_grown[d_idx] < 0.8
            and d_is_filler[d_idx] < 0.5
            and d_cs[d_idx] > self.cs_min
        ):
            den = d_rho_b_smooth[d_idx]
            if den > 1e-9:  # sem doador vivo na vizinhanca: evita 0/0
                local = d_rho_b_w2[d_idx] / den
                deficit = local - d_rho_b_grown[d_idx]
                if deficit > 0.0:
                    c_n = d_c_n[d_idx]
                    if c_n < 0.4:
                        c_n_factor = 0.0
                    elif c_n > 0.8:
                        c_n_factor = 1.0
                    else:
                        t = (c_n - 0.4) / 0.4
                        c_n_factor = t * t * (3.0 - 2.0 * t)

                    rate = self.k_col * c_n_factor * deficit
                    d_a_rho_b_grown[d_idx] += rate
                    d_am[d_idx] += d_m[d_idx] * rate / self.rho_max


class BiomassDiffusion(Equation):
    """Redistribuicao de biomassa DENTRO da fase densa (Rota D1).

    O degrau de rho_b na juncao nucleo-braco e um artefato de historia: o nucleo
    esta pinado e nao cresce, os bracos levaram a biomassa para fora, e rho_b nao
    tem nenhum termo de transporte alem da adveccao. Nada reequilibra o campo.

    Fundamentacao: [T2] Srinivasan 2019 trata a colonia como duas fases (ativa +
    passiva) com balanco de massa/momento, o que inclui FLUXO da fase ativa. Um
    Laplaciano de Brookshaw e a forma discreta desse fluxo.

    O gate `rho_b_i * rho_b_j` (produto, nao media) faz o coeficiente cair a zero
    na interface colonia-agar: so ha fluxo onde AMBAS as particulas sao densas.
    Sem isso a difusao borraria os bracos para dentro das baias e destruiria a
    morfologia dendritica (o mesmo modo de falha da licao #31).
    """

    def __init__(self, dest, sources, D_b):
        self.D_b = D_b
        super(BiomassDiffusion, self).__init__(dest, sources)

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        s_rho,
        d_rho_b_grown,
        s_rho_b_grown,
        d_a_rho_b_grown,
        d_h,
        RIJ,
        XIJ,
        DWIJ,
    ):
        gate = d_rho_b_grown[d_idx] * s_rho_b_grown[s_idx]
        if gate > 0.01:
            rb_ij = d_rho_b_grown[d_idx] - s_rho_b_grown[s_idx]
            rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2
            dot = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]
            vol_j = s_m[s_idx] / s_rho[s_idx]
            d_a_rho_b_grown[d_idx] += (
                2.0 * self.D_b * gate * vol_j * (rb_ij / rij_sq) * dot
            )


class ChemotacticFlux(Equation):
    """Fluxo quimiotatico de BIOMASSA — a metade que faltava da quimiotaxia.

    [T3] Giverso, Verani & Ciarletta 2016 compara crescimento volumetrico contra
    fluxo quimiotatico `m = χ·ρ·∇n` e mostra que o segundo produz "padroes mais
    simetricos com multiplos dendritos". O modelo tinha quimiotaxia so como FORCA
    sobre particulas (`FlagellarForce`); como `rho_b` e escalar passivo advectado,
    ela movia sempre os mesmos ~150 portadores e nunca levava biomassa a territorio
    novo. Dai o estado absorvente da licao #48.

        u = χ·(−∇cs)        deriva quimiotatica (agar fresco, §3.2 Frente 5)
        dρ_b/dt = −∇·(ρ_b·u)

    Forma antissimetrica par-a-par: o termo do par (i,j) em i e o negativo do termo
    em j, porque `∇_j W_ji = −∇_i W_ij`. Conserva `Σ V·dρ_b` (Liu §3.4; Violeau §5.3
    — atencao: o invariante NAO e `Σ ρ_b·V`, que muda com a deformacao do fluido).
    Particula em `rho_b = 0` tem fluxo de SAIDA zero e recebe do vizinho — e assim
    que o estado absorvente cai, sem termo aditivo.

    UPWIND obrigatorio. A forma simetrica `(ρ_i u_i + ρ_j u_j)` conserva no contínuo
    mas nao preserva POSITIVIDADE: num passo de Euler `rho_b` fica negativo onde o
    fluxo diverge, e o clamp `max(0, rho_b)` do integrador vira termo-FONTE, disparando
    a cada passo. Medido no teste de conservacao com a forma simetrica: vazamento de
    **+14% em 5 s**, com 23 739 particulas sentadas em zero recebendo fluxo liquido
    negativo. Com upwind o doador perde no maximo o que tem (`dρ_b/dt ≥ −C·ρ_b`,
    decaimento exponencial), `rho_b` nunca fica negativo e o clamp nunca dispara.

    Gate SIMETRICO em `rho_b < rho_b_pin`: assimetrico quebraria a conservacao (o
    vizinho ganharia o que o nucleo nao perde). Exclui o nucleo maduro, que nao
    quimiotaxa (§2.4) e cujo fluxo o dreno em τ≈41 s — o pin cinematico nao protege,
    porque zera a velocidade da particula, nao o fluxo do escalar.

    Roda em Group SEPARADO apos `equations_main`: `grad_cs_x/y` sao acumulados por
    `FlagellarForce.loop` e so estao finalizados no fim daquele grupo.
    """

    def __init__(self, dest, sources, chi, rho_b_pin=0.8, rho_target=0.4):
        self.chi = chi
        self.rho_b_pin = rho_b_pin
        self.rho_target = rho_target
        super().__init__(dest, sources)

    def loop(
        self,
        d_idx,
        s_idx,
        s_m,
        s_rho,
        d_rho_b_grown,
        s_rho_b_grown,
        d_grad_cs_x,
        d_grad_cs_y,
        s_grad_cs_x,
        s_grad_cs_y,
        d_is_filler,
        s_is_filler,
        d_a_rho_b_grown,
        XIJ,
        DWIJ,
    ):
        if (
            d_rho_b_grown[d_idx] < self.rho_b_pin
            and s_rho_b_grown[s_idx] < self.rho_b_pin
            and d_is_filler[d_idx] < 0.5
            and s_is_filler[s_idx] < 0.5
        ):
            ux = -self.chi * (d_grad_cs_x[d_idx] + s_grad_cs_x[s_idx])
            uy = -self.chi * (d_grad_cs_y[d_idx] + s_grad_cs_y[s_idx])
            # UPWIND: rho_b de quem o fluxo SAI. XIJ = x_i - x_j, entao o fluxo vai
            # de i para j quando u aponta de i para j, i.e. u.XIJ < 0. A escolha e
            # simetrica no par (para o destino j a condicao inverte junto com XIJ),
            # o que preserva a antissimetria e portanto a conservacao.
            if ux * XIJ[0] + uy * XIJ[1] < 0.0:
                rb_up = d_rho_b_grown[d_idx]
                cap = self.rho_target - s_rho_b_grown[s_idx]
            else:
                rb_up = s_rho_b_grown[s_idx]
                cap = self.rho_target - d_rho_b_grown[d_idx]
            # CAPACIDADE: o receptor so aceita ate `rho_target`. Sem isso o fluxo
            # espalha sem limite e dilui abaixo dos gates (licao #54: banda flagelar
            # 2041->34, colonia congelou). Com o limite a frente avanca como DEGRAU
            # (onda de preenchimento), que e como frentes de colonia real avancam.
            # O limite e propriedade do PAR (doador e receptor sao os mesmos qualquer
            # que seja o destino avaliado), entao a antissimetria — e a conservacao —
            # sobrevivem.
            if cap < 0.0:
                cap = 0.0
            if rb_up > cap:
                rb_up = cap
            vol_j = s_m[s_idx] / s_rho[s_idx]
            d_a_rho_b_grown[d_idx] -= vol_j * rb_up * (ux * DWIJ[0] + uy * DWIJ[1])


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
        d_is_env,
        s_is_env,
        DWIJ,
        d_grad_rho_b_x,
        d_grad_rho_b_y,
    ):
        # o PISO DE ENVELOPE ([0, 0.01)) e lido como ZERO aqui: ele nao pode mover o
        # `grad_rho_b_mag`, que e o gate de interface da MarangoniForce. Sem isto o piso
        # nao e inerte — foi o que custou 9% de AR no E6.
        rb_s = s_rho_b_grown[s_idx]
        rb_d = d_rho_b_grown[d_idx]
        if s_is_env[s_idx] > 0.5:
            rb_s = 0.0
        if d_is_env[d_idx] > 0.5:
            rb_d = 0.0

        vol_j = s_m[s_idx] / s_rho[s_idx]

        diff = rb_s - rb_d

        d_grad_rho_b_x[d_idx] += vol_j * diff * DWIJ[0]
        d_grad_rho_b_y[d_idx] += vol_j * diff * DWIJ[1]

    def post_loop(self, d_idx, d_grad_rho_b_x, d_grad_rho_b_y, d_grad_rho_b_mag):
        gx = d_grad_rho_b_x[d_idx]
        gy = d_grad_rho_b_y[d_idx]

        d_grad_rho_b_mag[d_idx] = (gx * gx + gy * gy) ** 0.5


class KernelSum(Equation):
    # Soma dos volumes dos vizinhos ponderados pelo kernel
    def initialize(self, d_idx, d_sigma_a):
        d_sigma_a[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho, WIJ, d_sigma_a):
        rho_safe = s_rho[s_idx]  # densidade da particula vizinha
        if rho_safe < 1e-6:  # minimo para evitar NaN
            rho_safe = 1e-6
        d_sigma_a[d_idx] += (s_m[s_idx] / rho_safe) * WIJ


class KernelGradientCorrection(Equation):
    def __init__(self, dest, sources, det_min=0.25):
        self.det_min = det_min
        super().__init__(dest, sources)

    def initialize(self, d_idx, d_Mxx, d_Mxy, d_Myx, d_Myy):
        d_Mxx[d_idx] = 0.0
        d_Mxy[d_idx] = 0.0
        d_Myx[d_idx] = 0.0
        d_Myy[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho, XIJ, DWIJ, d_Mxx, d_Mxy, d_Myx, d_Myy):
        vol_j = s_m[s_idx] / s_rho[s_idx]
        dxx = -XIJ[0]
        dyy = -XIJ[1]
        d_Mxx[d_idx] += vol_j * dxx * DWIJ[0]
        d_Mxy[d_idx] += vol_j * dxx * DWIJ[1]
        d_Myx[d_idx] += vol_j * dyy * DWIJ[0]
        d_Myy[d_idx] += vol_j * dyy * DWIJ[1]

    def post_loop(self, d_idx, d_Mxx, d_Mxy, d_Myx, d_Myy, d_Lxx, d_Lxy, d_Lyx, d_Lyy):
        det = d_Mxx[d_idx] * d_Myy[d_idx] - d_Mxy[d_idx] * d_Myx[d_idx]
        if det < self.det_min:
            d_Lxx[d_idx] = 1.0
            d_Lxy[d_idx] = 0.0
            d_Lyx[d_idx] = 0.0
            d_Lyy[d_idx] = 1.0
        else:
            inv_det = 1.0 / det
            d_Lxx[d_idx] = d_Myy[d_idx] * inv_det
            d_Lxy[d_idx] = -d_Mxy[d_idx] * inv_det
            d_Lyx[d_idx] = -d_Myx[d_idx] * inv_det
            d_Lyy[d_idx] = d_Mxx[d_idx] * inv_det


class ParticleShift(Equation):
    def __init__(
        self,
        dest,
        sources,
        shift_coeff,
        shift_cap,
        rho_b_min=0.6,
        rho_b_pin=0.8,
        c_n_pin=0.6,
    ):
        self.shift_coeff = shift_coeff
        self.shift_cap = shift_cap
        self.rho_b_min = rho_b_min
        self.rho_b_pin = rho_b_pin
        self.c_n_pin = c_n_pin
        super().__init__(dest, sources)

    def initialize(self, d_idx, d_shift_dC_x, d_shift_dC_y):
        d_shift_dC_x[d_idx] = 0.0
        d_shift_dC_y[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho, DWIJ, d_shift_dC_x, d_shift_dC_y):
        # gradiente da "concentração de partículas"
        vol_j = s_m[s_idx] / s_rho[s_idx]
        d_shift_dC_x[d_idx] += vol_j * DWIJ[0]
        d_shift_dC_y[d_idx] += vol_j * DWIJ[1]

    def post_loop(
        self,
        d_idx,
        d_shift_dC_x,
        d_shift_dC_y,
        d_shift_x,
        d_shift_y,
        d_rho_b_grown,
        d_is_env,
        d_c_n,
        d_h,
    ):
        rho_b = d_rho_b_grown[d_idx]
        if d_is_env[d_idx] > 0.5:
            rho_b = 0.0
        sx = 0.0
        sy = 0.0
        # C2: cláusula c_n removida. Com D_n=0.05 o gate ja abria sozinho em ~25% da
        # banda; o clumping (Liu §6.4) esta no restante, que continuava sem shifting.
        if rho_b >= self.rho_b_min and rho_b < self.rho_b_pin:
            h = d_h[d_idx]
            D = self.shift_coeff * h * h
            sx = (
                -D * d_shift_dC_x[d_idx]
            )  # O sinal negativo manda a partícula para longe do aglomerado, em direção ao buraco
            sy = -D * d_shift_dC_y[d_idx]

            mag = (sx * sx + sy * sy) ** 0.5
            cap = self.shift_cap * h
            if mag > cap:
                sx = sx * cap / mag
                sy = sy * cap / mag

        d_shift_x[d_idx] = sx
        d_shift_y[d_idx] = sy


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
        hill_k=0.1,
        lambda_bio_ratio=2.0,
    ):
        self.hill_k2 = hill_k * hill_k
        self.lambda_bio = lambda_ * lambda_bio_ratio
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
        d_is_filler,
        s_is_filler,
        d_is_matrix,
        s_is_matrix,
    ):
        if (d_is_filler[d_idx] < 0.5 or d_is_matrix[d_idx] > 0.5) and (
            s_is_filler[s_idx] < 0.5 or s_is_matrix[s_idx] > 0.5
        ):
            rho_b_avg = 0.5 * (d_rho_b_grown[d_idx] + s_rho_b_grown[s_idx])
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
        d_is_filler,
        d_is_matrix,
        d_is_env,
    ):
        if d_is_env[d_idx] > 0.5:
            # marcador de envelope: conduz cs (fica no loop) mas nao produz. Decai como
            # agar (0.5*lambda) — nao tem celula para degradar ramnolipideo.
            d_a_c_s[d_idx] += -0.5 * self.lambda_ * d_cs[d_idx]
        elif d_is_matrix[d_idx] > 0.5:
            d_a_c_s[d_idx] += -0.5 * self.lambda_ * d_cs[d_idx]
        elif d_is_filler[d_idx] > 0.5:
            d_a_c_s[d_idx] = 0.0
        else:
            rho_b = d_rho_b_grown[d_idx]
            qs = rho_b * rho_b / (rho_b * rho_b + self.hill_k2)

            saturation = 1.0 - d_cs[d_idx] / self.cs_max
            if saturation < 0.0:
                saturation = 0.0

            c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + 0.1)

            production = self.sigma * qs * saturation * d_noise[d_idx] * c_n_factor

            # CLAUDE.md §7 / Pass T1: decay LENTO no agar (gera o halo), rapido no
            # biofilme. O codigo tinha o oposto.
            if rho_b < 0.1:
                lambda_eff = self.lambda_ * 0.5
            else:
                lambda_eff = self.lambda_bio

            d_a_c_s[d_idx] += production - lambda_eff * d_cs[d_idx]


class MarangoniForce(Equation):
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
        d_is_filler,
        s_is_filler,
        d_is_matrix,
        s_is_matrix,
        d_Lxx,
        d_Lxy,
        d_Lyx,
        d_Lyy,
        DWIJ,
    ):
        grad_mag = d_grad_rho_b_mag[d_idx]
        if (
            grad_mag >= self.grad_low
            and (d_is_filler[d_idx] < 0.5 or d_is_matrix[d_idx] > 0.5)
            and (s_is_filler[s_idx] < 0.5 or s_is_matrix[s_idx] > 0.5)
        ):
            gate_raw = (grad_mag - self.grad_low) / (self.grad_high - self.grad_low)
            if gate_raw > 1.0:
                gate = 1.0
            else:
                gate = gate_raw * gate_raw * (3.0 - 2.0 * gate_raw)

            vol_j = s_m[s_idx] / s_rho[s_idx]
            cs_ij = s_cs[s_idx] - d_cs[d_idx]

            cdwij_x = d_Lxx[d_idx] * DWIJ[0] + d_Lxy[d_idx] * DWIJ[1]
            cdwij_y = d_Lyx[d_idx] * DWIJ[0] + d_Lyy[d_idx] * DWIJ[1]

            acc_x = gate * self.beta * vol_j * cs_ij * cdwij_x
            acc_y = gate * self.beta * vol_j * cs_ij * cdwij_y

            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_ax_mar[d_idx] += acc_x
            d_ay_mar[d_idx] += acc_y

    def post_loop(self, d_idx, d_au_mar, d_ax_mar, d_ay_mar):
        d_au_mar[d_idx] = (
            d_ax_mar[d_idx] * d_ax_mar[d_idx] + d_ay_mar[d_idx] * d_ay_mar[d_idx]
        ) ** 0.5


class LinearDrag(Equation):
    """Arrasto linear. `agar_ratio` reduz o arrasto onde nao ha biomassa.

    Diagnostico que motivou (2026-08-14): o agar engolido pela colonia nao e
    empurrado NEM arrastado. Nao e empurrado porque `BiomassEOS` zera `fade` abaixo
    de `rho_b=0.1` e ele fica com `|p| = 0` EXATO mesmo comprimido a `rho/rho0` 7.1.
    Nao e arrastado porque `gamma_base=60` se aplica a ele tambem: a forca viscosa
    que a colonia exerce (`~mu*v/h^2 = 1.1e-3`) da velocidade terminal `1.8e-5`,
    contra `~5e-4` da colonia — 4%. Medido: `|v|` mediano do agar = 4.9e-29.

    Fisicamente `gamma` e friccao flagelo-substrato, propriedade da BACTERIA. Aplica-la
    ao meio e o que o transforma em fundo rigido. Com `agar_ratio=0.1` a velocidade
    terminal do agar sobe ~10x e ele passa a acompanhar a frente em vez de ser engolido.
    """

    def __init__(self, dest, sources, gamma_base, gamma_mature=0.0, agar_ratio=1.0):
        self.gamma_base = -gamma_base
        self.gamma_mature = -gamma_mature
        self.agar_ratio = agar_ratio
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
        # smoothstep em [0.1, 0.3] para nao criar descontinuidade de arrasto na
        # interface colonia-agar, que geraria cisalhamento numerico ali.
        if rho_b < 0.1:
            g_base = self.gamma_base * self.agar_ratio
        elif rho_b < 0.3:
            t = (rho_b - 0.1) / 0.2
            s = t * t * (3.0 - 2.0 * t)
            g_base = self.gamma_base * (self.agar_ratio + (1.0 - self.agar_ratio) * s)
        else:
            g_base = self.gamma_base
        gamma_eff = g_base + self.gamma_mature * rho_b * rho_b
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

        rho_avg = 0.5 * (d_rho[d_idx] + rho_j)
        mu_eff = self.mu * min(rho_avg, 1.0)

        acc_x = 2.0 * mu_eff * common_term * u_ij
        acc_y = 2.0 * mu_eff * common_term * v_ij

        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y


class BiomassEOS(Equation):
    def __init__(
        self, dest, sources, rho0, c0, gamma_eos=7.0, tension_ratio=0.02, agar_fade=0.0
    ):
        self.rho0 = rho0
        self.c0 = c0
        self.B = rho0 * c0 * c0 / gamma_eos
        self.B_tension = self.B * tension_ratio
        self.agar_fade = agar_fade
        super(BiomassEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_rho_b_grown):
        ratio = d_rho[d_idx] / self.rho0
        rho_b = d_rho_b_grown[d_idx]

        # `agar_fade` da REPULSAO ao meio sem biomassa (atracao continua zero: agar nao
        # e coeso). Com agar_fade=0 o agar tem `|p| = 0` EXATO mesmo comprimido a
        # `rho/rho0 = 4.8`, entao nao transmite empurrao: a colonia o ATRAVESSA em vez
        # de desloca-lo, e ele e engolido (23029 particulas dentro da colonia no C4).
        #
        # Como `p = B*excess^2`, ligar isso NAO pressuriza o dominio: o agar
        # nao-comprimido esta em `rho/rho0 = 1.00` e continua com `p ~ 0`. So o
        # comprimido responde — que e onde se quer resposta. Medido: no anel de contato
        # daria `p` max 7.0e-3 contra 1.5e-3 da colonia, mesma ordem.
        #
        # Tem de estar ligado desde t=0. Ligar num run ja formado liberaria de uma vez
        # a energia armazenada nos `rho/rho0` de 4.8 (C4) a 16.7 (A1).
        if rho_b < 0.1:
            fade_rep = self.agar_fade
            fade_att = 0.0
        elif rho_b < 0.5:
            t = (rho_b - 0.1) / 0.4
            s = t * t * (3.0 - 2.0 * t)
            # continuo na interface: sai de `agar_fade` e chega a 1 em rho_b=0.5
            fade_rep = self.agar_fade + (1.0 - self.agar_fade) * s
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
            excess = ratio - 1.0
            d_p[d_idx] = self.B * excess * excess * fade_rep
        else:
            deficit = 1.0 - ratio
            if deficit > 0.3:
                deficit = 0.3
            d_p[d_idx] = -self.B_tension * deficit * fade_att


class OsmoticForce(Equation):
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

        if rho_b >= self.rho_b_threshold:
            gate_raw = (rho_b - self.rho_b_threshold) / 0.4
            if gate_raw > 1.0:
                gate = 1.0
            else:
                gate = gate_raw * gate_raw * (3.0 - 2.0 * gate_raw)

            vol_j = s_m[s_idx] / s_rho[s_idx]
            rho_diff = s_rho[s_idx] - d_rho[d_idx]

            acc_x = self.k_osm * gate * vol_j * rho_diff * DWIJ[0]
            acc_y = self.k_osm * gate * vol_j * rho_diff * DWIJ[1]

            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_au_osm[d_idx] += (acc_x * acc_x + acc_y * acc_y) ** 0.5


class OsmolyteProduction(Equation):
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

        qs = rho_b * rho_b / (rho_b * rho_b + 0.01)
        production = self.k_o * qs * (1.0 - d_c_o[d_idx])

        decay = self.lambda_o * d_c_o[d_idx]

        d_a_c_o[d_idx] += production - decay

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
        s_is_filler,
        s_is_matrix,
        DWIJ,
    ):
        if s_is_filler[s_idx] < 0.5 or s_is_matrix[s_idx] > 0.5:
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
        d_is_filler,
    ):
        rho_b = d_rho_b_grown[d_idx]
        # Gate: swarmers na borda apenas, pico em rho_b=0.35
        if rho_b >= 0.2 and rho_b <= 0.6 and d_is_filler[d_idx] < 0.5:
            if rho_b < 0.4:
                t = (rho_b - 0.2) / 0.2
            else:
                t = (0.6 - rho_b) / 0.2
            gate = t * t * (3.0 - 2.0 * t)

            gx = d_grad_cs_x[d_idx]
            gy = d_grad_cs_y[d_idx]
            mag = (gx * gx + gy * gy) ** 0.5 + 1e-9

            acc_x = -self.f0 * gate * gx / mag
            acc_y = -self.f0 * gate * gy / mag
            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_au_flag[d_idx] = (acc_x * acc_x + acc_y * acc_y) ** 0.5


class OxigenConsumption(Equation):
    def __init__(
        self, dest, sources, D_n=0.02, D_n_int=1e-4, k_n=0.5, filler_transparent=0
    ):
        self.D_n = D_n
        self.D_n_int = D_n_int
        self.k_n = k_n
        self.filler_transparent = filler_transparent
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

    def post_loop(
        self, d_idx, d_a_c_n, d_c_n, d_rho_b_grown, d_is_filler, d_is_matrix, d_is_env
    ):
        # A banda [0, 0.01) e o PISO DE ENVELOPE, nao celula metabolizante: marca o agar
        # que a colonia ja ocupou para que o campo `rho_b` apareca continuo. Sem esta
        # isencao, 8712 particulas a 1e-3 somariam 66% ao consumo real de nutriente.
        if d_is_env[d_idx] > 0.5:
            consumption = 0.0
        elif self.filler_transparent == 1 and d_is_filler[d_idx] > 0.5:
            consumption = 0.0
        elif d_is_matrix[d_idx] > 0.5:
            consumption = 0.0
        else:
            consumption = self.k_n * d_rho_b_grown[d_idx] * d_c_n[d_idx]
        d_a_c_n[d_idx] -= consumption


class NutrientSource(Equation):
    """Reposicao do nutriente pelo agar — regime NUTRIENT-RICH.

    [T2] Srinivasan, Kaplan & Mahadevan 2019 (eLife 8, e42697): swarming de
    P. aeruginosa e regime nutrient-rich, com `c ~ c0` constante. Um `c_n` que
    esgota monotonicamente descreve BIOFILME — que a propria §3.0 identifica como
    o regime errado para este organismo. Fisicamente: a placa e um reservatorio 3D
    e a simulacao e um corte 2D, entao a reposicao vertical nao esta representada.

    Sem este termo a janela util termina em t ~ 55 s (licao #47) e o crescimento
    so consegue um fator `exp(r_growth*T) = 2.46x` — a semente e o resultado.

    Equilibrio local sob consumo:  c_n = k_src / (k_src + k_n * rho_b).

    Roda no post_loop APOS `OxigenConsumption` (ordem das equacoes no Group), que
    e quem zera `d_a_c_n` no initialize.
    """

    def __init__(self, dest, sources, k_src):
        self.k_src = k_src
        super(NutrientSource, self).__init__(dest, sources)

    def post_loop(self, d_idx, d_a_c_n, d_c_n):
        d_a_c_n[d_idx] += self.k_src * (1.0 - d_c_n[d_idx])
