from pysph.sph.scheme import Scheme
from pysph.sph.integrator import EulerIntegrator
from pysph.sph.integrator_step import EulerStep
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
from pysph.sph.wc.basic import MomentumEquation
from .equations import (
    BiomassEOS,
    BiomassGrowth,
    FlagellarForce,
    KernelGradientCorrection,
    KernelSum,
    MarangoniForce,
    OxigenConsumption,
    ParticleShift,
    SurfactantEquation,
    LinearDrag,
    ViscousForce,
    BiomassGradient,
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
        d_c_o,
        d_a_c_o,
        d_c_n,
        d_a_c_n,
        d_shift_x,
        d_shift_y,
        d_is_filler,
        dt,
    ):
        d_rho_b_grown[d_idx] += dt * d_a_rho_b_grown[d_idx]
        d_cs[d_idx] += dt * d_a_c_s[d_idx]
        d_c_o[d_idx] += dt * d_a_c_o[d_idx]
        d_c_n[d_idx] += dt * d_a_c_n[d_idx]

        d_rho_b_grown[d_idx] = max(0.0, min(d_rho_b_grown[d_idx], 1.0))
        d_cs[d_idx] = max(1e-9, d_cs[d_idx])
        d_c_o[d_idx] = max(0.0, min(d_c_o[d_idx], 1.0))
        d_c_n[d_idx] = max(1e-9, min(d_c_n[d_idx], 1.0))

        if (
            d_rho_b_grown[d_idx] >= 0.8
            or d_c_n[d_idx] < 0.6
            or d_is_filler[d_idx] > 0.5
        ):
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
        p0=0.0,
        D_o=1e-3,
        k_o=0.5,
        lambda_o=0.05,
        Q0=5e-4,
        D_n=0.02,
        D_n_int=1e-4,
        k_n=0.5,
        use_shift=False,
        shift_coeff=0.5,
        shift_cap=0.05,
        shift_rho_b_min=0.6,
        use_kgc=False,
        kgc_det_min=0.25,
    ):
        self.use_shift = use_shift
        self.shift_coeff = shift_coeff
        self.shift_cap = shift_cap
        self.shift_rho_b_min = shift_rho_b_min
        self.use_kgc = use_kgc
        self.kgc_det_min = kgc_det_min
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
        self.p0 = p0
        self.D_o = D_o
        self.k_o = k_o
        self.lambda_o = lambda_o
        self.Q0 = Q0
        self.D_n = D_n
        self.D_n_int = D_n_int
        self.k_n = k_n
        super(MyBiomassScheme, self).__init__(fluids, solids, dim=dim)

    def get_equations(self):
        equations_pre = Group(
            equations=[
                SummationDensity(dest="fluid", sources=["fluid"]),
                BiomassEOS(
                    dest="fluid",
                    sources=None,
                    rho0=1.0,
                    c0=self.c0,
                    tension_ratio=0.30,
                ),
            ],
            real=False,
        )

        equations_kernel_sum = Group(
            equations=[
                KernelSum(dest="fluid", sources=["fluid"]),
            ],
            real=False,
        )

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
                # BiomassGradient DEVE vir antes de MarangoniForce
                # (Marangoni usa grad_rho_b_mag como gate de interface)
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
                    lambda_ext_ratio=5.0,
                    k_consume=0.0,
                    cs_max=0.5,
                ),
                OxigenConsumption(
                    dest="fluid",
                    sources=["fluid"],
                    D_n=self.D_n,
                    D_n_int=self.D_n_int,
                    k_n=self.k_n,
                ),
                # OsmolyteProduction(
                #     dest="fluid",
                #     sources=["fluid"],
                #     D_o=self.D_o,
                #     k_o=self.k_o,
                #     lambda_o=self.lambda_o,
                #     Q0=self.Q0,
                # ),
                FlagellarForce(
                    dest="fluid",
                    sources=["fluid"],
                    f0=3.0,
                ),
            ],
        )

        groups = [equations_pre, equations_kernel_sum]

        if self.use_kgc:
            equations_kgc = Group(
                equations=[
                    KernelGradientCorrection(
                        dest="fluid",
                        sources=["fluid"],
                        det_min=self.kgc_det_min,
                    ),
                ],
                real=False,
            )
            groups.append(equations_kgc)

        groups.append(equations_main)

        if self.use_shift:
            equations_shift = Group(
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
            groups.append(equations_shift)

        return groups

    def get_integrator(self):
        return EulerIntegrator(fluid=CustomEulerStep())
