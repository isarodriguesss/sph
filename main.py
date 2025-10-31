import numpy as np
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver
from pysph.base.utils import get_particle_array

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme

x_dim, y_dim = 64, 64
# x_dim, y_dim = 128, 128  # Dimensões originais

x_min_domain, x_max_domain = -1.0, 5.0
y_min_domain, y_max_domain = -1.0, 5.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)

# expansões maiores
# x_min_domain, x_max_domain = -6.0, 6.0
# y_min_domain, y_max_domain = -6.0, 6.0

mu = 1.0
# mu = 0.07  # Valor para evitar instabilidade
gamma = 50.0
# gamma = 13.0  # Valor para evitar instabilidade
beta = 0.1
# beta = 0.5  # Valor para evitar instabilidade
sigma = 1.0
D = 0.001
lambda_ = 0.1
r_growth = 1.0
rho_max = 1.0

dt_global = 0.001
total_sim_time = 20.0
# total_sim_time = 100.0 # expensões maiores
print_freq = 500
# print_freq = 2500 # expensões maiores

trajectory_store_interval = 20

prob_of_splitting = 0.007


class SwarmApp(Application):
    def initialize(self):
        self.n_bact = 300
        self.all_bact_trajectories = [[] for _ in range(self.n_bact)]

    def create_particles(self):
        fluid_solid = create_initial_state(x_dim=x_dim, y_dim=y_dim, rho_max=rho_max)
        x = np.linspace(x_min_domain, x_max_domain, x_dim)
        y = np.linspace(y_min_domain, y_max_domain, y_dim)
        center_x = (x.min() + x.max()) / 2.0
        center_y = (y.min() + y.max()) / 2.0
        r0 = (x[1] - x[0]) * 2
        angles = np.linspace(0, 2 * np.pi, self.n_bact, endpoint=False)

        xb = center_x + r0 * np.cos(angles)
        yb = center_y + r0 * np.sin(angles)

        bact = get_particle_array(
            name="bact",
            x=xb,
            y=yb,
            m=np.ones_like(xb),
            h=np.ones_like(xb) * (x[1] - x[0]) * 1.2,
        )

        bact.add_property("u", default=0.0)
        bact.add_property("v", default=0.0)

        bact.add_property("rho", default=1.0)

        self.particles = fluid_solid + [bact]

        for pa in self.particles:
            if pa.name == "fluid":
                pa.add_output_arrays(["rho_b_grown", "cs", "u", "v", "p"])
            elif pa.name == "bact":
                pa.add_output_arrays(["u", "v"])

        return self.particles

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
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        solver = Solver(
            dim=2,
            integrator=self.scheme.get_integrator(),
            kernel=kernel,
            dt=1e-3,
            adaptive_timestep=False,
        )
        solver.tf = total_sim_time
        solver.set_print_freq(print_freq)
        return solver

    def post_step(self, solver):
        dt = solver.dt
        bact = self.particles[2]

        bact.x += dt * bact.u
        bact.y += dt * bact.v

        domain_width = x_max_domain - x_min_domain
        domain_height = y_max_domain - y_min_domain

        bact.x[:] = (bact.x - x_min_domain) % domain_width + x_min_domain
        bact.y[:] = (bact.y - y_min_domain) % domain_height + y_min_domain

        if solver.count % trajectory_store_interval == 0:
            for i in range(self.n_bact):
                self.all_bact_trajectories[i].append([bact.y[i], bact.x[i]])

        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            max_vel = np.max(np.sqrt(fluid.u**2 + fluid.v**2))

            print(
                f"  t={solver.t:.2e}s ({((solver.t) / total_sim_time) * 100:.1f}%), "
                f"dt={solver.dt:.2e}s, Max c_s: {np.max(fluid.cs):.2e}, "
                f"Max rho_b: {np.max(fluid.rho_b_grown):.2f}, Max |v|: {max_vel:.2e}"
            )

        if solver.count > 0 and solver.count % 100 == 0:
            fluid = self.particles[0]
            min_rho_b_for_division = 0.05
            mature_by_mass_indices = np.where(fluid.m > 1.99 * fluid.m0)[0]
            mature_by_rho_b_indices = np.where(fluid.rho_b_grown[mature_by_mass_indices] > min_rho_b_for_division)[0]
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
                        prop: getattr(fluid, prop)[parent_idx] for prop in props_to_copy
                    }

                    for i in range(2):
                        daughter_data = parent_props.copy()
                        daughter_data["m"] = parent_props["m"] / 2.0
                        daughter_data["m0"] = parent_props["m0"]
                        daughter_data["rho_b_grown"] = parent_props["rho_b_grown"] / 2.0
                        dx_local = parent_props["h"] / 4.0
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


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
