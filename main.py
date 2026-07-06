import csv
import numpy as np
from scipy.spatial import cKDTree
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme

LOG_FILE = "log.csv"
LOG_HEADER = [
    "t",
    "iteration",
    "max_v",
    "mean_v",
    "n_fast",
    "a_marangoni",
    "a_drag",
    "a_pressure",
    "a_flag",
    "a_total",
    "min_cs",
    "max_cs",
    "mean_cs",
    "constrast_cs",
    "min_c_n",
    "max_c_n",
    "mean_c_n",
    "contrast_c_n",
    "mass_total",
    "pass_n_spawned",
    "min_sig_arms",
    "mean_sig_arms",
    "frac_lowsig_arms",
]

x_dim, y_dim = 187, 187

x_min_domain, x_max_domain = -5.0, 5.0
y_min_domain, y_max_domain = -5.0, 5.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)


mu = 0.020
gamma = 60.0
beta = 5.0
sigma = 10.0
D = 1.5e-3
D_ext = 0.08
lambda_ = 0.15
r_growth = 0.02
rho_max = 1.0
alpha_mon = 0.12

D_o = 0.04
k_o = 0.5
lambda_o = 0.05
Q0 = 5.0

D_n = 0.02
D_n_int = 1e-4
k_n = 0.5

dt_global = 0.001
total_sim_time = 100.0
print_freq = 200

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35

use_splitting = False

use_shift = False
SHIFT_COEFF = 0.5
SHIFT_CAP = 0.05
SHIFT_RHO_B_MIN = 0.6

use_kgc = False
KGC_DET_MIN = 0.25

use_insert = True
INSERT_FREQ = 200
INSERT_SIGMA_TRIG = 0.85 # 15% de erro
INSERT_RHO_B_MIN = 0.1  # só vácuo estrutural no núcleo. 0.1 para tudo
INSERT_RHO_B_CORE = 0.5  # rho_b > core: núcleo (frozen isotrópico); abaixo: braço
INSERT_TANGENTIAL_EPS = 0.5  # |offset·n̂| < eps*dx → candidato tangencial ao braço
INSERT_PROX = 0.7
INSERT_MAX = 100

use_pass_n = False
PASS_N_FREQ = 100
PASS_N_MAX_PARENTS = 100
PASS_N_SIGMA_TRIG = 0.95
PASS_N_ALPHA = 0.75
PASS_N_EPSILON = 0.5
PASS_N_MAX_GEN = 1
PASS_N_RHO_B_MIN = 0.15
PASS_N_RHO_B_MAX = 0.95
PASS_N_PROXIMITY_MIN = 0.4


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = (
            None  # snapshot da massa total em t=0 (preenchido em post_step)
        )
        self._pass_n_spawned_since_log = 0  # acumula spawns entre linhas de log

    def create_particles(self):
        fluid_solid = create_initial_state(
            x_dim,
            y_dim,
            rho_max,
            dt_global,
            x_min=x_min_domain,
            x_max=x_max_domain,
            y_min=y_min_domain,
            y_max=y_max_domain,
        )

        for pa in fluid_solid:
            if pa.name == "fluid":
                pa.add_property("noise")
                pa.noise[:] = (
                    1.0
                    + 0.6 * np.sin(8 * np.arctan2(pa.y, pa.x))
                    + 0.01 * np.random.rand(len(pa.x))
                )
                pa.add_property("dt_force")
                pa.add_property("dt_cfl")
                pa.add_property("au_mar")
                pa.add_property("ax_mar")
                pa.add_property("ay_mar")
                pa.add_property("au_drag")
                pa.add_property("grad_rho_b_x")
                pa.add_property("grad_rho_b_y")
                pa.add_property("grad_rho_b_mag")
                pa.add_property("grad_cs_x")
                pa.add_property("grad_cs_y")
                pa.add_property("grad_co_x")
                pa.add_property("grad_co_y")
                pa.add_property("sigma_a")
                pa.sigma_a[:] = 1.0
                pa.add_property("shift_dC_x")
                pa.add_property("shift_dC_y")
                pa.add_property("shift_x")
                pa.add_property("shift_y")
                pa.shift_x[:] = 0.0
                pa.shift_y[:] = 0.0
                pa.add_property("Mxx")
                pa.add_property("Mxy")
                pa.add_property("Myx")
                pa.add_property("Myy")
                pa.add_property("Lxx")
                pa.add_property("Lxy")
                pa.add_property("Lyx")
                pa.add_property("Lyy")
                pa.Lxx[:] = 1.0
                pa.Lxy[:] = 0.0
                pa.Lyx[:] = 0.0
                pa.Lyy[:] = 1.0
                pa.add_property("is_filler")
                pa.is_filler[:] = 0.0
                pa.add_property("gen")
                pa.gen[:] = 0.0
                pa.add_property("ax_drag")
                pa.add_property("ay_drag")
                pa.add_property("au_flag")
                pa.add_property("x_spawn_ref")
                pa.add_property("y_spawn_ref")
                pa.x_spawn_ref[:] = pa.x[:]
                pa.y_spawn_ref[:] = pa.y[:]
                pa.add_output_arrays(
                    [
                        "rho_b_grown",
                        "cs",
                        "c_o",
                        "u",
                        "v",
                        "p",
                        "noise",
                        "au_flag",
                        "au_mar",
                        "sigma_a",
                    ]
                )
            elif pa.name == "solid":
                pa.add_property("p")

        return fluid_solid

    def create_scheme(self):
        return MyBiomassScheme(
            fluids=["fluid"],
            solids=["solid"],
            dim=2,
            mu=mu,
            gamma=gamma,
            beta=beta,
            sigma=sigma,
            D=D,
            D_ext=D_ext,
            lambda_=lambda_,
            r_growth=r_growth,
            rho_max=rho_max,
            c0=c0,
            alpha_mon=alpha_mon,
            D_o=D_o,
            k_o=k_o,
            lambda_o=lambda_o,
            Q0=Q0,
            D_n=D_n,
            D_n_int=D_n_int,
            k_n=k_n,
            use_shift=use_shift,
            shift_coeff=SHIFT_COEFF,
            shift_cap=SHIFT_CAP,
            shift_rho_b_min=SHIFT_RHO_B_MIN,
            use_kgc=use_kgc,
            kgc_det_min=KGC_DET_MIN,
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        solver = Solver(
            dim=2,
            integrator=self.scheme.get_integrator(),
            kernel=kernel,
            dt=5e-5,
            adaptive_timestep=True,
            cfl=0.4,
        )
        solver.tf = total_sim_time
        solver.set_print_freq(print_freq)
        return solver

    def post_step(self, solver):
        # Print stats
        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            v_mag = np.sqrt(fluid.u**2 + fluid.v**2)
            max_v = np.max(v_mag)
            mean_v = np.mean(v_mag)
            n_fast = int(np.sum(v_mag > 0.1))

            a_mar = np.max(np.abs(fluid.au_mar))
            a_drag = np.max(np.abs(fluid.au_drag))
            a_total = np.max(np.sqrt(fluid.au**2 + fluid.av**2))
            ax_p = fluid.au - fluid.ax_mar - fluid.ax_drag
            ay_p = fluid.av - fluid.ay_mar - fluid.ay_drag
            a_pressure = np.max(np.sqrt(ax_p**2 + ay_p**2))
            a_flag = np.max(np.abs(fluid.au_flag))

            min_cs = np.min(fluid.cs)
            max_cs = np.max(fluid.cs)
            mean_cs = np.mean(fluid.cs)
            contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

            min_c_n = np.min(fluid.c_n)
            max_c_n = np.max(fluid.c_n)
            mean_c_n = np.mean(fluid.c_n)
            contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

            mass_total = float(np.sum(fluid.m))

            arms_mask = (fluid.rho_b_grown >= 0.1) & (fluid.rho_b_grown < 0.5)
            n_arms = int(np.sum(arms_mask))
            if n_arms > 0:
                sig_arms = fluid.sigma_a[arms_mask]
                min_sig_arms = float(np.min(sig_arms))
                mean_sig_arms = float(np.mean(sig_arms))
                frac_lowsig_arms = float(np.mean(sig_arms < 0.85))
            else:
                min_sig_arms = mean_sig_arms = 1.0
                frac_lowsig_arms = 0.0

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
            )
            print(
                f"sigma_a braços (rho_b∈[0.1,0.5]): min={min_sig_arms:.3f} "
                f"mean={mean_sig_arms:.3f} frac<0.85={frac_lowsig_arms:.2%}"
            )
            print("Acelerações:")
            print(f"  > Marangoni (líq): {a_mar:.2f}")
            print(f"  > Drag:            {a_drag:.2f} (Freio)")
            print(f"  > Pressão (est):   {a_pressure:.2f}")
            print(f"  > Total |a|:       {a_total:.2f}")

            with open(LOG_FILE, "a", newline="") as f:
                csv.writer(f).writerow(
                    [
                        f"{solver.t:.4f}",
                        solver.count,
                        f"{max_v:.6f}",
                        f"{mean_v:.6f}",
                        n_fast,
                        f"{a_mar:.4f}",
                        f"{a_drag:.4f}",
                        f"{a_pressure:.4f}",
                        f"{a_flag:.4f}",
                        f"{a_total:.4f}",
                        f"{min_cs:.4f}",
                        f"{max_cs:.4f}",
                        f"{mean_cs:.4f}",
                        f"{contrast_cs:.4f}",
                        f"{min_c_n:.4f}",
                        f"{max_c_n:.4f}",
                        f"{mean_c_n:.4f}",
                        f"{contrast_c_n:.4f}",
                        f"{mass_total:.6e}",
                        self._pass_n_spawned_since_log,
                        f"{min_sig_arms:.4f}",
                        f"{mean_sig_arms:.4f}",
                        f"{frac_lowsig_arms:.4f}",
                    ]
                )
                self._pass_n_spawned_since_log = 0  # reseta após registrar

        if use_splitting:
            if solver.count > 0 and solver.count % 500 == 0:
                fluid = self.particles[0]
                min_rho_b_for_division = 0.05
                mature_by_mass_indices = np.where(fluid.m > 1.99 * fluid.m0)[0]
                mature_by_rho_b_indices = np.where(
                    fluid.rho_b_grown[mature_by_mass_indices] > min_rho_b_for_division
                )[0]
                indices_to_split = []
                for idx in mature_by_rho_b_indices:
                    original_idx = mature_by_mass_indices[idx]
                    if np.random.rand() < prob_of_splitting:
                        indices_to_split.append(original_idx)

                if len(indices_to_split) > 0:
                    print(
                        f"\n--- Divisão Celular em t={solver.t:.2f}: {len(indices_to_split)} partículas se dividindo. ---"
                    )

                    daughters = fluid.empty_clone()

                    props_to_copy = [
                        "x",
                        "y",
                        "m",
                        "h",
                        "rho",
                        "rho_b_grown",
                        "cs",
                        "u",
                        "v",
                        "au",
                        "av",
                        "a_rho_b_grown",
                        "a_c_s",
                        "m0",
                    ]

                    for parent_idx in indices_to_split:
                        parent_props = {
                            prop: getattr(fluid, prop)[parent_idx]
                            for prop in props_to_copy
                        }

                        for i in range(2):
                            daughter_data = parent_props.copy()
                            daughter_data["m"] = parent_props["m"] / 2.0
                            daughter_data["m0"] = parent_props["m0"]
                            # daughter_data["rho_b_grown"] = parent_props["rho_b_grown"] / 2.0
                            dx_local = parent_props["h"] * 0.5
                            offset_x = dx_local * (np.random.rand() - 0.5)
                            offset_y = dx_local * (np.random.rand() - 0.5)

                            daughter_data = parent_props.copy()
                            daughter_data["x"] = parent_props["x"] + offset_x
                            daughter_data["y"] = parent_props["y"] + offset_y

                            data_to_add = {
                                key: [value] for key, value in daughter_data.items()
                            }

                            daughters.add_particles(**data_to_add)

                    fluid.append_parray(daughters)

                    fluid.remove_particles(indices_to_split)

                    solver.nnps.update()

        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            splittable_gen_mask = fluid.gen < PASS_N_MAX_GEN

            sigma_a = fluid.sigma_a
            split_mask = (
                colony_mask & (sigma_a < PASS_N_SIGMA_TRIG) & splittable_gen_mask
            )
            split_idx = np.where(split_mask)[0]

            if len(split_idx) > 0:
                if len(split_idx) > PASS_N_MAX_PARENTS:
                    order = np.argsort(sigma_a[split_idx])
                    split_idx = split_idx[order][:PASS_N_MAX_PARENTS]

                eps = PASS_N_EPSILON
                alpha = PASS_N_ALPHA
                prox_min = PASS_N_PROXIMITY_MIN * dx
                prox_min_sq = prox_min * prox_min

                # Padrão hexagonal: 6 vértices a 60°
                angles_hex = np.arange(6) * (np.pi / 3.0)
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                tree_mask = np.ones(len(fluid.x), dtype=bool)
                tree_mask[split_idx] = False
                positions = np.column_stack([fluid.x[tree_mask], fluid.y[tree_mask]])
                tree = cKDTree(positions)

                daughters = fluid.empty_clone()
                new_positions = []  # filhas já adicionadas nesta call
                mothers_used = []  # índices de mães que efetivamente splitaram
                n_daughters_total = 0

                # Snapshot dos campos (índices estáveis até o remove)
                x_arr = fluid.x[split_idx]
                y_arr = fluid.y[split_idx]
                h_arr = fluid.h[split_idx]
                m_arr = fluid.m[split_idx]
                u_arr = fluid.u[split_idx]
                v_arr = fluid.v[split_idx]
                rho_arr = fluid.rho[split_idx]
                rhob_arr = fluid.rho_b_grown[split_idx]
                cs_arr = fluid.cs[split_idx]
                co_arr = fluid.c_o[split_idx]
                cn_arr = fluid.c_n[split_idx]
                noise_arr = fluid.noise[split_idx]
                sigma_arr = fluid.sigma_a[split_idx]
                gen_arr = fluid.gen[split_idx]

                for k in range(len(split_idx)):
                    x_m, y_m = float(x_arr[k]), float(y_arr[k])
                    h_m = float(h_arr[k])
                    m_m = float(m_arr[k])
                    r_off = eps * h_m

                    valid_vertices = []
                    for j in range(6):
                        vx = x_m + r_off * cos_a[j]
                        vy = y_m + r_off * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox_min:
                            continue
                        too_close = False
                        for nx, ny in new_positions:
                            if (vx - nx) ** 2 + (vy - ny) ** 2 < prox_min_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        valid_vertices.append((vx, vy))

                    n_d = 1 + len(valid_vertices)
                    if n_d != 7:
                        continue

                    m_d = m_m / 7.0
                    h_d = alpha * h_m
                    u_d = float(u_arr[k])
                    v_d = float(v_arr[k])

                    xs = [x_m] + [p[0] for p in valid_vertices]
                    ys = [y_m] + [p[1] for p in valid_vertices]

                    new_positions.append((x_m, y_m))
                    new_positions.extend(valid_vertices)

                    data = {
                        "x": xs,
                        "y": ys,
                        "m": [m_d] * 7,
                        "h": [h_d] * 7,
                        "rho": [float(rho_arr[k])] * 7,
                        "rho_b_grown": [float(rhob_arr[k])] * 7,
                        "cs": [float(cs_arr[k])] * 7,
                        "c_o": [float(co_arr[k])] * 7,
                        "c_n": [float(cn_arr[k])] * 7,
                        "u": [u_d] * 7,
                        "v": [v_d] * 7,
                        "noise": [float(noise_arr[k])] * 7,
                        "sigma_a": [float(sigma_arr[k])] * 7,
                        "is_filler": [0.0] * 7,
                        "gen": [float(gen_arr[k]) + 1.0] * 7,
                    }
                    daughters.add_particles(**data)
                    mothers_used.append(int(split_idx[k]))
                    n_daughters_total += 7

                if mothers_used:
                    fluid.append_parray(daughters)
                    fluid.remove_particles(np.asarray(mothers_used, dtype=np.uint32))
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_daughters_total
                    print(
                        f"Pass N v2.4 t={solver.t:.1f}s: "
                        f"{len(mothers_used)} mães → {n_daughters_total} filhas "
                        f"(n_d=7 estrito, ε={eps}, α={alpha}, σ_trig={PASS_N_SIGMA_TRIG}, gen≤1)"
                    )

        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx

            void_mask = (fluid.rho_b_grown > INSERT_RHO_B_MIN) & (
                fluid.sigma_a < INSERT_SIGMA_TRIG
            )
            void_idx = np.where(void_mask)[0]

            if len(void_idx) > 0:
                positions = np.column_stack([fluid.x, fluid.y])
                tree = cKDTree(positions)
                prox = INSERT_PROX * dx  # distância de segurança entre duas partículas
                prox_sq = prox * prox  # o quadrado da distância de segurança

                tang_eps = INSERT_TANGENTIAL_EPS * dx

                angles_hex = np.arange(6) * (
                    np.pi / 3.0
                )  # angulos hexagonais 0°, 60°, 120°, 180°, 240°, 300°
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                new_x = []
                new_y = []
                parent_idx = []
                new_u = []
                new_v = []
                new_filler = []
                added_pts = []  # dedupe entre candidatos da mesma call

                for k in void_idx:
                    xk = float(fluid.x[k])
                    yk = float(fluid.y[k])

                    is_arm = float(fluid.rho_b_grown[k]) < INSERT_RHO_B_CORE
                    nx_hat = 0.0
                    ny_hat = 0.0
                    if is_arm:
                        gx = float(fluid.grad_rho_b_x[k])
                        gy = float(fluid.grad_rho_b_y[k])
                        gmag = (gx * gx + gy * gy) ** 0.5
                        if gmag > 1e-9:
                            nx_hat = gx / gmag
                            ny_hat = gy / gmag
                        else:
                            is_arm = False

                    for j in range(6):
                        ox = dx * cos_a[j]
                        oy = dx * sin_a[j]
                        if is_arm:
                            proj = ox * nx_hat + oy * ny_hat
                            if proj < 0.0:
                                proj = -proj
                            if proj > tang_eps:
                                continue
                        vx = xk + ox
                        vy = yk + oy
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox:
                            continue
                        too_close = False
                        for ax_, ay_ in added_pts:
                            if (vx - ax_) ** 2 + (vy - ay_) ** 2 < prox_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        new_x.append(vx)
                        new_y.append(vy)
                        parent_idx.append(int(k))
                        if is_arm:
                            new_u.append(float(fluid.u[k]))
                            new_v.append(float(fluid.v[k]))
                            new_filler.append(2.0)
                        else:
                            new_u.append(0.0)
                            new_v.append(0.0)
                            new_filler.append(1.0)
                        added_pts.append((vx, vy))
                    if len(new_x) >= INSERT_MAX:
                        break

                if len(new_x) > 0:
                    n_ins = len(new_x)
                    p = np.asarray(parent_idx, dtype=int)
                    inserted = fluid.empty_clone()
                    data = {
                        "x": new_x,
                        "y": new_y,
                        "m": [m_target] * n_ins,
                        "h": list(fluid.h[p]),
                        "rho": list(fluid.rho[p]),
                        "rho_b_grown": list(fluid.rho_b_grown[p]),
                        "cs": list(fluid.cs[p]),
                        "c_o": list(fluid.c_o[p]),
                        "c_n": list(fluid.c_n[p]),
                        "u": new_u,
                        "v": new_v,
                        "noise": list(fluid.noise[p]),
                        "is_filler": new_filler,
                    }
                    inserted.add_particles(**data)
                    fluid.append_parray(inserted)
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_ins
                    print(
                        f"C3 inserção t={solver.t:.1f}s: {n_ins} partículas inseridas"
                    )


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
