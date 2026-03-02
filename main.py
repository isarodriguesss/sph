import numpy as np
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme

x_dim, y_dim = 100, 100
# x_dim, y_dim = 128, 128  # Dimensões originais

x_min_domain, x_max_domain = -1.0, 5.0
y_min_domain, y_max_domain = -1.0, 5.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)

# expansões maiores
# x_min_domain, x_max_domain = -6.0, 6.0
# y_min_domain, y_max_domain = -6.0, 6.0

mu = 0.001
# mu = 0.02  # Valor anterior
gamma = 7.0
# gamma = 40.0  # Valor anterior
beta = 10.0
# beta = 0.5  # Valor para evitar instabilidade
sigma = 4.0
# sigma = 2.0  # Valor anterior
D = 1e-5
lambda_ = 0.01
# lambda_ = 0.05 # Valor anterior
r_growth = 2.0
# r_growth = 1.5  # Valor anterior
rho_max = 1.0

dt_global = 0.001
total_sim_time = 20.0
# total_sim_time = 100.0 # expensões maiores
print_freq = 200
# print_freq = 2500 # expensões maiores

trajectory_store_interval = 20

prob_of_splitting = 0.03
c0 = 10.0
# c0 = 1.0  # Valor anterior

use_splitting = False


class SwarmApp(Application):
    def initialize(self):
        pass

    def create_particles(self):
        fluid_solid = create_initial_state(x_dim, y_dim, rho_max, dt_global)

        for pa in fluid_solid:
            if pa.name == "fluid":
                pa.add_property("noise")
                pa.noise[:] = 1.0 + 0.1 * (np.random.rand(len(pa.x)))
                pa.add_property("dt_force")
                pa.add_property("dt_cfl")
                pa.add_output_arrays(["rho_b_grown", "cs", "u", "v", "p", "noise"])
                # debug da aceleração
                pa.add_property("au_mar")  # aceleração x Marangoni
                pa.add_property("au_osm")  # aceleração x Osmótica
                pa.add_property("au_pres")  # aceleração x Pressão (estimada)
                pa.add_property("au_drag")  # aceleração x Drag
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
            lambda_=lambda_,
            r_growth=r_growth,
            rho_max=rho_max,
            c0=c0,
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        solver = Solver(
            dim=2,
            integrator=self.scheme.get_integrator(),
            kernel=kernel,
            dt=1e-4,
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
            max_v = np.max(np.sqrt(fluid.u**2 + fluid.v**2))

            # Pegamos o valor absoluto máximo de cada componente de aceleração
            a_mar = np.max(np.abs(fluid.au_mar))
            a_osm = np.max(np.abs(fluid.au_osm))
            a_drag = np.max(np.abs(fluid.au_drag))
            # O resto da aceleração au vem da MomentumEquation (pressão)
            a_total = np.max(np.abs(fluid.au))

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print("Acelerações (au):")
            print(f"  > Marangoni: {a_mar:.2f}")
            print(f"  > Osmótica:  {a_osm:.2f}")
            print(f"  > Drag:      {a_drag:.2f} (Freio)")
            print(f"  > Total au:  {a_total:.2f}")

            # Resetamos os acumuladores de print para o próximo passo
            fluid.au_mar[:] = 0.0
            fluid.au_osm[:] = 0.0
            fluid.au_drag[:] = 0.0

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
