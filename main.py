import csv
import numpy as np
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
    # Pass M-B — campo de nutriente consumivel c_n (Frente 6)
    "min_c_n",
    "max_c_n",
    "mean_c_n",
    "contrast_c_n",
    "mass_total",  # soma de m — cresce por BiomassGrowth (logistico)
]

x_dim, y_dim = 150, 150  # Pass I.8: resolucao aumentada (era 100x100, dx 0.06→0.04)

# Domínio 6×6 centrado na origem
x_min_domain, x_max_domain = -3.0, 3.0
y_min_domain, y_max_domain = -3.0, 3.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)


# ── Parâmetros Calibrados (Coesão Viscosa + Interface Gateada) ──────────
# Meta: v_term ≈ 0.1, acelerações totais 10–50, Marangoni só na interface.
# Coesão vem de: (a) EOS com ramo atrativo (p<0 se rho<rho0),
#                (b) viscosidade maior, (c) Monaghan artificial viscosity forte,
#                (d) kernel com ~35 vizinhos (h_factor=1.8).
mu = 0.020  # K.23: I.2 revert parcial — fortalecer coesão viscosa (I.2 era 0.012, pre-I.2 era 0.025)
gamma = 60.0  # Drag: a_drag = gamma * v_term = 60 * 0.1 = 6
beta = 1.0  # Marangoni (com gate de interface, só ~60% ativo em média)
sigma = 1.5  # Pass I.7: boost +67% compensa drenagem por D_ext (era 1.2)
D = 1.5e-3  # Pass I.3: D_int dentro do biofilme — gradiente afiado na interface
D_ext = 0.08  # K.18: 4x — L_D_ext=0.298 (~2x maior); habilita focalizacao Mullins-Sekerka pos-K.17
lambda_ = 0.15  # Decaimento: confina cs mas permite penetracao de ~L_D_ext no exterior
r_growth = 0.05  # K.21: 0.4→0.15 — estende vida do swarmer ring (~60s→~150s) evitando trap K.17 quando rho_b satura globalmente
rho_max = 1.0
alpha_mon = 0.12  # K.23: I.2 revert parcial — previne instabilidade de tração SPH (I.2 era 0.06, pre-I.2 era 0.15)

# Pass M-A (Frente 6): osmolito c_o secretado pelas bacterias.
# Mecanismo: c_o satura nas baias (agar confinado entre dendritos), fica ~0
# nas pontas (agar virgem). |∇c_o| dispara influxo de massa via van't Hoff.
D_o = 0.04  # Pass M-A.2: L_D_o = sqrt(0.04/0.05) = 0.894 ~ d_tip-tip (Mullins-Sekerka)
k_o = 0.5  # taxa de producao por bacteria
lambda_o = 0.05  # decaimento lento (osmolitos persistem mais que cs)
Q0 = 5.0  # Pass M-A.2: forca osmotica Darcy — a_osm_tip ~ Q0*gate*|grad_c_o| ~ 2.8

# Pass M-B (Frente 6): nutriente consumivel c_n.
# c_n inicia em 1.0 (agar virgem). Bacterias consomem na taxa k_n*rho_b.
# Difusao tem que dominar consumo (tau_cons/tau_dif > 4) para evitar morte
# quimica global. Lição da rodada inicial M-B (k_n=1.0, D_n=1e-3): consumo
# dominava → motor cs morria em t<2s. Calibracao corrigida: razao = 25.
D_n = 0.02  # difusao do nutriente no agar (~D_ext/4)
k_n = 0.5  # taxa de consumo por unidade de biomassa

dt_global = 0.001
total_sim_time = 100.0
print_freq = 200

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 0.35  # EOS: B = 1.5²/7 ≈ 0.32 (repulsão suave, atração ~0.1)

use_splitting = False


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = (
            None  # snapshot da massa total em t=0 (preenchido em post_step)
        )

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
                # debug das acelerações (componentes vetoriais + magnitude)
                pa.add_property("au_mar")  # |aceleração Marangoni| (líquida)
                pa.add_property("ax_mar")  # aceleração Marangoni componente x
                pa.add_property("ay_mar")  # aceleração Marangoni componente y
                pa.add_property("au_drag")  # |aceleração Drag|
                # gradiente da densidade
                pa.add_property("grad_rho_b_x")
                pa.add_property("grad_rho_b_y")
                pa.add_property("grad_rho_b_mag")
                # gradiente do surfactante (usado por FlagellarForce)
                pa.add_property("grad_cs_x")
                pa.add_property("grad_cs_y")
                # gradiente do osmolito c_o (Pass M-A — usado por OsmolyteProduction)
                pa.add_property("grad_co_x")
                pa.add_property("grad_co_y")
                pa.add_property("ax_drag")
                pa.add_property("ay_drag")
                # flag
                pa.add_property("au_flag")
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
            k_n=k_n,
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
        # dt = solver.dt
        # bact = self.particles[2]

        # bact.x += dt * bact.u
        # bact.y += dt * bact.v

        # domain_width = x_max_domain - x_min_domain
        # domain_height = y_max_domain - y_min_domain

        # bact.x[:] = (bact.x - x_min_domain) % domain_width + x_min_domain
        # bact.y[:] = (bact.y - y_min_domain) % domain_height + y_min_domain

        # if solver.count % trajectory_store_interval == 0:
        #     for i in range(self.n_bact):
        #         self.all_bact_trajectories[i].append([bact.y[i], bact.x[i]])

        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            v_mag = np.sqrt(fluid.u**2 + fluid.v**2)
            max_v = np.max(v_mag)
            mean_v = np.mean(v_mag)
            n_fast = int(np.sum(v_mag > 0.1))

            # Marangoni líquido (magnitude do vetor, não soma de normas)
            a_mar = np.max(np.abs(fluid.au_mar))
            a_drag = np.max(np.abs(fluid.au_drag))
            # Aceleração total (inclui pressão via MomentumEquation + tudo)
            a_total = np.max(np.sqrt(fluid.au**2 + fluid.av**2))
            # ax_pressure = au_total - ax_mar - ax_drag
            ax_p = fluid.au - fluid.ax_mar - fluid.ax_drag
            ay_p = fluid.av - fluid.ay_mar - fluid.ay_drag
            a_pressure = np.max(np.sqrt(ax_p**2 + ay_p**2))
            a_flag = np.max(np.abs(fluid.au_flag))

            # 3. Estatísticas do Surfactante (cs)
            min_cs = np.min(fluid.cs)
            max_cs = np.max(fluid.cs)
            mean_cs = np.mean(fluid.cs)
            contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

            # 4. Pass M-B — Estatísticas do nutriente c_n
            # c_n inicia em 1.0 (agar virgem). Bacterias consomem na taxa k_n*rho_b.
            # mean_c_n esperado: cair de 1.0 mas estabilizar em ~0.3-0.5 (difusao
            # repoe nutriente do exterior). Se cai para <0.1, motor cs vai morrer.
            # contrast_c_n alto = baias esgotadas vs pontas em agar virgem (selecao).
            min_c_n = np.min(fluid.c_n)
            max_c_n = np.max(fluid.c_n)
            mean_c_n = np.mean(fluid.c_n)
            contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

            # 5. Massa total
            # Cresce por BiomassGrowth (logistico). Sem bug osmotico desde M-A.2.
            mass_total = float(np.sum(fluid.m))

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
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
                    ]
                )

        if use_splitting:
            print("Iniciando processo de divisão celular...")
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

        else:
            pass


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
