from scipy.interpolate import griddata
import numpy as np
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import create_initial_state
from src.scheme import MyBiomassScheme

x_dim, y_dim = 48, 48
# x_dim, y_dim = 128, 128  # Dimensões originais

x_min_domain, x_max_domain = -1.0, 5.0
y_min_domain, y_max_domain = -1.0, 5.0

mu = 0.05
gamma = 10.0
beta = 0.5
sigma = 1.0
D = 0.01
lambda_ = 0.1
r_growth = 0.2
rho_max = 1.0

dt_global = 0.001
total_sim_time = 10.0
print_freq = 100


class SwarmApp(Application):
    def initialize(self):
        self.n_bact = 300
        x = np.linspace(x_min_domain, x_max_domain, x_dim)
        y = np.linspace(y_min_domain, y_max_domain, y_dim)
        center_x = (x.min() + x.max()) / 2.0
        center_y = (y.min() + y.max()) / 2.0
        r0 = (x[1] - x[0]) * 2
        angles = np.linspace(0, 2 * np.pi, self.n_bact, endpoint=False)
        self.bact_positions = np.array(
            [[center_y + r0 * np.sin(a), center_x + r0 * np.cos(a)] for a in angles]
        )
        self.all_bact_trajectories = [[] for _ in range(self.n_bact)]

    def create_particles(self):
        self.particles = create_initial_state(rho_max=rho_max, dt=dt_global)
        for pa in self.particles:
            if pa.name == "fluid":
                pa.add_output_arrays(["rho_b_grown", "cs", "u", "v", "p"])
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
        solver = Solver(dim=2, integrator=self.scheme.get_integrator(), kernel=kernel)
        solver.set_adaptive_timestep(True)
        solver.cfl = 0.3
        solver.tf = total_sim_time
        solver.set_print_freq(print_freq)
        return solver

    def post_step(self, solver):
        fluid = self.particles[0]
        x_min, x_max = -1, 5
        y_min, y_max = -1, 5

        points = np.vstack((fluid.x, fluid.y)).T

        v_bact_x = griddata(
            points,
            fluid.u,
            (self.bact_positions[:, 1], self.bact_positions[:, 0]),
            method="nearest",
        )
        v_bact_y = griddata(
            points,
            fluid.v,
            (self.bact_positions[:, 1], self.bact_positions[:, 0]),
            method="nearest",
        )

        dt = solver.dt
        self.bact_positions[:, 1] += v_bact_x * dt
        self.bact_positions[:, 0] += v_bact_y * dt

        self.bact_positions[:, 1] = (self.bact_positions[:, 1] - x_min) % (
            x_max - x_min
        ) + x_min
        self.bact_positions[:, 0] = (self.bact_positions[:, 0] - y_min) % (
            y_max - y_min
        ) + y_min

        if solver.count % (print_freq * 2) == 0:
            print(
                f"  Progresso: {((solver.t) / total_sim_time) * 100:.1f}% completo. "
                f"Max cs: {np.max(fluid.cs):.2e}, Max rho_b: {np.max(fluid.rho_b_grown):.2f}, "
                f"Max |v|: {np.max(np.sqrt(fluid.u**2 + fluid.v**2)):.2e}"
            )


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
