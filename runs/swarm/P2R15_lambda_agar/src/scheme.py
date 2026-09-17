from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
from pysph.sph.wc.basic import MomentumEquation
from .equations import (
    BiomassColonization,
    BiomassEOS,
    BiomassGradient,
    BiomassGrowth,
    FlagellarForce,
    KernelGradientCorrection,
    KernelSum,
    LinearDrag,
    MarangoniForce,
    NutrientSource,
    OxigenConsumption,
    ParticleShift,
    SurfactantEquation,
    ViscousForce,
)


class CustomEulerStep(EulerStep):
    def stage1(
        self,
        d_idx,
        d_m,
        d_am,
        d_u,
        d_v,
        d_au,
        d_av,
        d_x,
        d_y,
        d_rho_b_grown,
        d_a_rho_b_grown,
        d_cs,
        d_a_c_s,
        d_c_n,
        d_a_c_n,
        d_shift_x,
        d_shift_y,
        d_is_filler,
        d_is_wake,
        dt,
    ):
        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_cs[d_idx] += dt * d_a_c_s[d_idx]
        d_c_n[d_idx] += dt * d_a_c_n[d_idx]

        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0))
        d_cs[d_idx] = max(1e-9, d_cs[d_idx])
        d_c_n[d_idx] = max(1e-9, min(d_c_n[d_idx], 1.0))

        if (
            d_rho_b_grown[d_idx] >= 0.8
            or d_c_n[d_idx] < 0.6
            or d_is_filler[d_idx] > 0.5
        ) and d_is_wake[d_idx] < 0.5:
            d_u[d_idx] = 0.0
            d_v[d_idx] = 0.0
        else:
            d_u[d_idx] += dt * d_au[d_idx]
            d_v[d_idx] += dt * d_av[d_idx]

        d_x[d_idx] += dt * d_u[d_idx]
        d_y[d_idx] += dt * d_v[d_idx]
        d_m[d_idx] += dt * d_am[d_idx]

        d_x[d_idx] += d_shift_x[d_idx]
        d_y[d_idx] += d_shift_y[d_idx]


class MyBiomassScheme(Scheme):
    def __init__(
        self,
        fluids,
        solids,
        dim,
        mu,
        gamma,
        beta,
        sigma,
        D,
        D_ext,
        lambda_,
        r_growth,
        rho_max,
        c0=10.0,
        alpha_mon=0.5,
        D_n=0.02,
        D_n_int=1e-4,
        k_n=0.5,
        k_src=0.0,
        k_col=0.0,
        cs_max=0.5,
        col_cs_min=0.3,
        col_filler_donor=0.0,
        hill_k=0.1,
        lambda_bio_ratio=2.0,
        flag_gate_lo=0.2,
        flag_gate_hi=0.6,
        flag_f0=3.0,
        use_shift=False,
        shift_coeff=0.5,
        shift_cap=0.05,
        shift_rho_b_min=0.6,
        use_kgc=False,
        kgc_det_min=0.25,
        filler_cs_conduz=0.0,
        filler_cs_D=0.0,
        filler_cs_D_interno=0.0,
        filler_cs_lambda=0.5,
        agar_cs_lambda=0.5,
    ):
        self.mu = mu
        self.gamma = gamma
        self.beta = beta
        self.sigma = sigma
        self.D = D
        self.D_ext = D_ext
        self.lambda_ = lambda_
        self.r_growth = r_growth
        self.rho_max = rho_max
        self.c0 = c0
        self.alpha_mon = alpha_mon
        self.D_n = D_n
        self.D_n_int = D_n_int
        self.k_n = k_n
        self.k_src = k_src
        self.k_col = k_col
        self.cs_max = cs_max
        self.col_cs_min = col_cs_min
        self.col_filler_donor = col_filler_donor
        self.hill_k = hill_k
        self.lambda_bio_ratio = lambda_bio_ratio
        self.flag_gate_lo = flag_gate_lo
        self.flag_gate_hi = flag_gate_hi
        self.flag_f0 = flag_f0
        self.use_shift = use_shift
        self.shift_coeff = shift_coeff
        self.shift_cap = shift_cap
        self.shift_rho_b_min = shift_rho_b_min
        self.use_kgc = use_kgc
        self.kgc_det_min = kgc_det_min
        self.filler_cs_conduz = filler_cs_conduz
        self.filler_cs_D = filler_cs_D
        self.filler_cs_D_interno = filler_cs_D_interno
        self.filler_cs_lambda = filler_cs_lambda
        self.agar_cs_lambda = agar_cs_lambda
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[SummationDensity(dest="fluid", sources=["fluid"])],
            real=False,
        )

        # EOS em Group proprio: precisa de `rho` ja finalizado pelo SummationDensity
        equations_eos = Group(
            equations=[
                BiomassEOS(
                    dest="fluid", sources=None, rho0=1.0, c0=self.c0, tension_ratio=0.30
                ),
            ],
            real=False,
        )

        equations_kernel_sum = Group(
            equations=[KernelSum(dest="fluid", sources=["fluid"])],
            real=False,
        )

        # ordem: BiomassGradient antes de MarangoniForce (gate de interface);
        # SurfactantEquation antes de FlagellarForce (cs atualizado)
        equations_main = Group(
            equations=[
                MomentumEquation(
                    dest="fluid",
                    sources=["fluid"],
                    c0=self.c0,
                    alpha=self.alpha_mon,
                    beta=0.0,
                ),
                BiomassGrowth(
                    dest="fluid",
                    sources=None,
                    r_growth=self.r_growth,
                    rho_max=self.rho_max,
                ),
                BiomassColonization(
                    dest="fluid",
                    sources=["fluid"],
                    k_col=self.k_col,
                    rho_max=self.rho_max,
                    cs_min=self.col_cs_min,
                    filler_donor=self.col_filler_donor,
                ),
                BiomassGradient(dest="fluid", sources=["fluid"]),
                MarangoniForce(dest="fluid", sources=["fluid"], beta=self.beta),
                ViscousForce(dest="fluid", sources=["fluid"], mu=self.mu),
                LinearDrag(
                    dest="fluid",
                    sources=None,
                    gamma_base=self.gamma,
                    gamma_mature=self.gamma * 1.5,
                ),
                SurfactantEquation(
                    dest="fluid",
                    sources=["fluid"],
                    D=self.D,
                    D_ext=self.D_ext,
                    sigma=self.sigma,
                    lambda_=self.lambda_,
                    cs_max=self.cs_max,
                    hill_k=self.hill_k,
                    lambda_bio_ratio=self.lambda_bio_ratio,
                    filler_conduz=self.filler_cs_conduz,
                    filler_D=self.filler_cs_D,
                    filler_D_interno=self.filler_cs_D_interno,
                    filler_lambda_ratio=self.filler_cs_lambda,
                    agar_lambda_ratio=self.agar_cs_lambda,
                ),
                OxigenConsumption(
                    dest="fluid",
                    sources=["fluid"],
                    D_n=self.D_n,
                    D_n_int=self.D_n_int,
                    k_n=self.k_n,
                ),
                NutrientSource(dest="fluid", sources=None, k_src=self.k_src),
                FlagellarForce(
                    dest="fluid",
                    sources=["fluid"],
                    f0=self.flag_f0,
                    gate_lo=self.flag_gate_lo,
                    gate_hi=self.flag_gate_hi,
                ),
            ],
        )

        groups = [equations_pre, equations_eos, equations_kernel_sum]

        if self.use_kgc:
            groups.append(
                Group(
                    equations=[
                        KernelGradientCorrection(
                            dest="fluid", sources=["fluid"], det_min=self.kgc_det_min
                        ),
                    ],
                    real=False,
                )
            )

        groups.append(equations_main)

        if self.use_shift:
            groups.append(
                Group(
                    equations=[
                        ParticleShift(
                            dest="fluid",
                            sources=["fluid"],
                            shift_coeff=self.shift_coeff,
                            shift_cap=self.shift_cap,
                            rho_b_min=self.shift_rho_b_min,
                        ),
                    ],
                )
            )

        return groups

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
