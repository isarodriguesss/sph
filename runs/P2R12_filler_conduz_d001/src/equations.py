from pysph.sph.equation import Equation


class BiomassGrowth(Equation):
    def __init__(self, dest, sources, r_growth, rho_max):
        self.r_growth = r_growth
        self.rho_max = rho_max
        super(BiomassGrowth, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho_b_grown, d_a_rho_b_grown, d_m, d_am, d_c_n, d_is_filler):
        d_a_rho_b_grown[d_idx] = 0.0
        d_am[d_idx] = 0.0
        if (
            d_rho_b_grown[d_idx] > 1e-12
            and d_rho_b_grown[d_idx] < 0.8
            and d_is_filler[d_idx] < 0.5
        ):
            c_n = d_c_n[d_idx]
            if c_n < 0.4:
                c_n_factor = 0.0
            elif c_n > 0.8:
                c_n_factor = 1.0
            else:
                t = (c_n - 0.4) / 0.4
                c_n_factor = t * t * (3.0 - 2.0 * t)

            rate = self.r_growth * (1.0 - d_rho_b_grown[d_idx] / self.rho_max) * c_n_factor
            d_a_rho_b_grown[d_idx] = rate * d_rho_b_grown[d_idx]
            d_am[d_idx] = d_m[d_idx] * d_a_rho_b_grown[d_idx] / self.rho_max


class BiomassColonization(Equation):
    def __init__(self, dest, sources, k_col, rho_max, cs_min, filler_donor=0.0):
        self.k_col = k_col
        self.rho_max = rho_max
        self.cs_min = cs_min
        self.filler_donor = filler_donor
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
        WIJ,
        d_rho_b_w2,
        d_rho_b_smooth,
    ):
        if s_is_filler[s_idx] < 0.5 or self.filler_donor > 0.5:
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
            if den > 1e-9:  # sem doador na vizinhanca: evita 0/0
                deficit = d_rho_b_w2[d_idx] / den - d_rho_b_grown[d_idx]
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


class KernelSum(Equation):
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
        if det < self.det_min:  # suporte degenerado: sem correcao, evita amplificar
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
    def __init__(self, dest, sources, shift_coeff, shift_cap, rho_b_min=0.6, rho_b_pin=0.8):
        self.shift_coeff = shift_coeff
        self.shift_cap = shift_cap
        self.rho_b_min = rho_b_min
        self.rho_b_pin = rho_b_pin
        super().__init__(dest, sources)

    def initialize(self, d_idx, d_shift_dC_x, d_shift_dC_y):
        d_shift_dC_x[d_idx] = 0.0
        d_shift_dC_y[d_idx] = 0.0

    def loop(self, d_idx, s_idx, s_m, s_rho, DWIJ, d_shift_dC_x, d_shift_dC_y):
        vol_j = s_m[s_idx] / s_rho[s_idx]
        d_shift_dC_x[d_idx] += vol_j * DWIJ[0]
        d_shift_dC_y[d_idx] += vol_j * DWIJ[1]

    def post_loop(
        self, d_idx, d_shift_dC_x, d_shift_dC_y, d_shift_x, d_shift_y, d_rho_b_grown, d_h
    ):
        rho_b = d_rho_b_grown[d_idx]
        sx = 0.0
        sy = 0.0
        if rho_b >= self.rho_b_min and rho_b < self.rho_b_pin:
            h = d_h[d_idx]
            D = self.shift_coeff * h * h
            sx = -D * d_shift_dC_x[d_idx]
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
        cs_max=0.5,
        hill_k=0.1,
        lambda_bio_ratio=2.0,
        filler_conduz=0.0,
        filler_D=0.0,
    ):
        self.hill_k2 = hill_k * hill_k
        self.filler_conduz = filler_conduz
        self.filler_D = filler_D
        self.lambda_bio = lambda_ * lambda_bio_ratio
        self.D = D
        self.D_ext = D_ext
        self.sigma = sigma
        self.lambda_ = lambda_
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
    ):
        par = 0
        if d_is_filler[d_idx] < 0.5 and s_is_filler[s_idx] < 0.5:
            par = 1
        elif self.filler_conduz > 0.5:
            par = 1
        if par == 1:
            rb_d = d_rho_b_grown[d_idx]
            if d_is_filler[d_idx] > 0.5:
                rb_d = 0.0
            rb_s = s_rho_b_grown[s_idx]
            if s_is_filler[s_idx] > 0.5:
                rb_s = 0.0
            rho_b_avg = 0.5 * (rb_d + rb_s)
            if rho_b_avg < 0.1:
                D_eff = self.D_ext
            elif rho_b_avg < 0.5:
                t = (rho_b_avg - 0.1) / 0.4
                gate = t * t * (3.0 - 2.0 * t)
                D_eff = self.D_ext + (self.D - self.D_ext) * gate
            else:
                D_eff = self.D
            if self.filler_D > 0.0:
                if d_is_filler[d_idx] > 0.5 or s_is_filler[s_idx] > 0.5:
                    D_eff = self.filler_D

            cs_ij = d_cs[d_idx] - s_cs[s_idx]
            rij_sq = RIJ**2 + 0.01 * d_h[d_idx] ** 2
            dot_product = XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]
            term = (s_m[s_idx] / s_rho[s_idx]) * (cs_ij / rij_sq) * dot_product
            d_a_c_s[d_idx] += 2.0 * D_eff * term

    def post_loop(self, d_idx, d_rho_b_grown, d_a_c_s, d_cs, d_noise, d_c_n, d_is_filler):
        if d_is_filler[d_idx] > 0.5:
            if self.filler_conduz > 0.5:
                sink_f = self.lambda_ * 0.5 * d_cs[d_idx]
                d_a_c_s[d_idx] -= sink_f
            else:
                d_a_c_s[d_idx] = 0.0
        else:
            rho_b = d_rho_b_grown[d_idx]
            qs = rho_b * rho_b / (rho_b * rho_b + self.hill_k2)

            saturation = 1.0 - d_cs[d_idx] / self.cs_max
            if saturation < 0.0:
                saturation = 0.0

            c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + 0.1)

            production = self.sigma * qs * saturation * d_noise[d_idx] * c_n_factor

            if rho_b < 0.1:
                lambda_eff = self.lambda_ * 0.5
            else:
                lambda_eff = self.lambda_bio

            sink = lambda_eff * d_cs[d_idx]
            d_a_c_s[d_idx] += production - sink


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
        d_Lxx,
        d_Lxy,
        d_Lyx,
        d_Lyy,
        DWIJ,
    ):
        grad_mag = d_grad_rho_b_mag[d_idx]
        if grad_mag >= self.grad_low and d_is_filler[d_idx] < 0.5 and s_is_filler[s_idx] < 0.5:
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
    def __init__(self, dest, sources, gamma_base, gamma_mature=0.0):
        self.gamma_base = -gamma_base
        self.gamma_mature = -gamma_mature
        super(LinearDrag, self).__init__(dest, sources)

    def loop(self, d_idx, d_u, d_v, d_au, d_av, d_au_drag, d_rho_b_grown, d_ax_drag, d_ay_drag):
        rho_b = d_rho_b_grown[d_idx]
        gamma_eff = self.gamma_base + self.gamma_mature * rho_b * rho_b
        acc_x = gamma_eff * d_u[d_idx]
        acc_y = gamma_eff * d_v[d_idx]
        d_ax_drag[d_idx] = acc_x
        d_ay_drag[d_idx] = acc_y
        d_au[d_idx] += acc_x
        d_av[d_idx] += acc_y
        d_au_drag[d_idx] = (acc_x * acc_x + acc_y * acc_y) ** 0.5


class ViscousForce(Equation):
    def __init__(self, dest, sources, mu):
        self.mu = mu
        super(ViscousForce, self).__init__(dest, sources)

    def initialize(self, d_idx, d_ax_vis, d_ay_vis):
        d_ax_vis[d_idx] = 0.0
        d_ay_vis[d_idx] = 0.0

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
        d_ax_vis,
        d_ay_vis,
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
        d_ax_vis[d_idx] += acc_x
        d_ay_vis[d_idx] += acc_y


class BiomassEOS(Equation):
    def __init__(self, dest, sources, rho0, c0, gamma_eos=7.0, tension_ratio=0.02):
        self.rho0 = rho0
        self.c0 = c0
        self.B = rho0 * c0 * c0 / gamma_eos
        self.B_tension = self.B * tension_ratio
        super(BiomassEOS, self).__init__(dest, sources)

    def loop(self, d_idx, d_rho, d_p, d_rho_b_grown):
        ratio = d_rho[d_idx] / self.rho0
        rho_b = d_rho_b_grown[d_idx]

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


class FlagellarForce(Equation):
    def __init__(self, dest, sources, f0, gate_lo=0.2, gate_hi=0.6):
        self.f0 = f0
        self.gate_lo = gate_lo
        self.gate_hi = gate_hi
        self.gate_mid = 0.5 * (gate_lo + gate_hi)
        super().__init__(dest, sources)

    def initialize(self, d_idx, d_au_flag, d_ax_flag, d_ay_flag, d_grad_cs_x, d_grad_cs_y):
        d_au_flag[d_idx] = 0.0
        d_ax_flag[d_idx] = 0.0
        d_ay_flag[d_idx] = 0.0
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
        DWIJ,
    ):
        if s_is_filler[s_idx] < 0.5:
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
        d_ax_flag,
        d_ay_flag,
        d_is_filler,
    ):
        rho_b = d_rho_b_grown[d_idx]
        if rho_b >= self.gate_lo and rho_b <= self.gate_hi and d_is_filler[d_idx] < 0.5:
            if rho_b < self.gate_mid:
                t = (rho_b - self.gate_lo) / (self.gate_mid - self.gate_lo)
            else:
                t = (self.gate_hi - rho_b) / (self.gate_mid - self.gate_lo)
            gate = t * t * (3.0 - 2.0 * t)

            gx = d_grad_cs_x[d_idx]
            gy = d_grad_cs_y[d_idx]
            mag = (gx * gx + gy * gy) ** 0.5 + 1e-9  # evita divisao por zero

            acc_x = -self.f0 * gate * gx / mag
            acc_y = -self.f0 * gate * gy / mag
            d_au[d_idx] += acc_x
            d_av[d_idx] += acc_y
            d_ax_flag[d_idx] = acc_x
            d_ay_flag[d_idx] = acc_y
            d_au_flag[d_idx] = (acc_x * acc_x + acc_y * acc_y) ** 0.5


class OxigenConsumption(Equation):
    def __init__(self, dest, sources, D_n=0.02, D_n_int=1e-4, k_n=0.5):
        self.D_n = D_n
        self.D_n_int = D_n_int
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
        consumption = self.k_n * d_rho_b_grown[d_idx] * d_c_n[d_idx]
        d_a_c_n[d_idx] -= consumption


class NutrientSource(Equation):
    def __init__(self, dest, sources, k_src):
        self.k_src = k_src
        super(NutrientSource, self).__init__(dest, sources)

    def post_loop(self, d_idx, d_a_c_n, d_c_n):
        d_a_c_n[d_idx] += self.k_src * (1.0 - d_c_n[d_idx])
