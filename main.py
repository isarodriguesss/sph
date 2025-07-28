import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import os

from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from scipy.interpolate import RegularGridInterpolator

from src.particles import create_biomass_surfactant_particles
from src.scheme import MyBiomassScheme

x_dim, y_dim = 128, 128
domain_min_x, domain_max_x = -1.0, 5.0
domain_min_y, domain_max_y = -1.0, 5.0
dx = (domain_max_x - domain_min_x) / x_dim
dy = (domain_max_y - domain_min_y) / y_dim
h = 1.2 * dx

mu = 0.05
gamma = 10.0
beta = 0.5
sigma = 1.0
D = 0.01
lambda_ = 0.1
r_growth = 0.2
rho_max = 1.0

dt_global = 0.001
total_sim_time = 40.0
trajectory_store_interval = max(1, int(0.02 / dt_global))

n_bact = 300
r0_bact_initial = dx * 2

class BiomassSurfactantApp(Application):
    def initialize(self):
        self.particles = create_biomass_surfactant_particles(
            domain_min_x, domain_max_x, domain_min_y, domain_max_y,
            h, (dx, dy), rho_max, x_dim, y_dim,
        )
        center_x = (domain_min_x + domain_max_x) / 2.0
        center_y = (domain_min_y + domain_max_y) / 2.0
        angles = np.linspace(0, 2 * np.pi, n_bact, endpoint=False)
        self.bact_positions = np.array([[center_y + r0_bact_initial * np.sin(a), center_x + r0_bact_initial * np.cos(a)] for a in angles])
        self.all_bact_trajectories = [[] for _ in range(n_bact)]

        self.periodic_domain_min = np.array([domain_min_x, domain_min_y, -h])
        self.periodic_domain_max = np.array([domain_max_x, domain_max_y, h])

        self.x_grid_imp, self.y_grid_imp = np.meshgrid(
            np.linspace(domain_min_x, domain_max_x, x_dim),
            np.linspace(domain_min_y, domain_max_y, y_dim)
        )

        self.imposed_vx_grid = -np.sin(self.y_grid_imp * np.pi / (domain_max_y - domain_min_y)) * np.cos(self.x_grid_imp * np.pi / (domain_max_x - domain_min_x)) * 0.1
        self.imposed_vy_grid = np.cos(self.y_grid_imp * np.pi / (domain_max_y - domain_min_y)) * np.sin(self.x_grid_imp * np.pi / (domain_max_x - domain_min_x)) * 0.1

    def create_particles(self):
        particles = create_biomass_surfactant_particles(
            domain_min_x, domain_max_x, domain_min_y, domain_max_y,
            h, (dx, dy), rho_max, x_dim, y_dim
        )

        for pa in particles:
            if pa.name == 'fluid': 
                pa.add_property('rho_b_grown', type='double', default=0.0)
                pa.add_property('c_s', type='double', default=0.0)
                pa.add_output_arrays(['rho', 'm', 'x', 'y', 'u', 'v', 'h', 'rho_b_grown', 'c_s', 'ax', 'ay'])

        return particles

    def create_scheme(self):
        return MyBiomassScheme(fluids=['fluid'], solids=[], dim=2,
                                rho_max=rho_max, r_growth=r_growth,
                                sigma=sigma, lambda_=lambda_, beta=beta, gamma=gamma, D=D, mu=mu,
                                periodic_domain=(self.periodic_domain_min, self.periodic_domain_max))
    
    def create_solver(self):
        kernel = CubicSpline(dim=2)
        scheme = self.create_scheme()
        solver = Solver(dim=2,
                        integrator=scheme.get_integrator(),
                        kernel=kernel,
                        tf=total_sim_time,
                        dt=dt_global,
                        set_adaptive_timestep=False,
                        set_print_freq=trajectory_store_interval,
                        )
        return solver

    def create_tools(self):
        return []
    
    def post_step(self, solver):
        fluid_array = self.particles[0]

        vx_interp_func_imposed = RegularGridInterpolator(
            (np.linspace(domain_min_y, domain_max_y, y_dim), np.linspace(domain_min_x, domain_max_x, x_dim)),
            self.imposed_vx_grid, bounds_error=False, fill_value=0
        )
        vy_interp_func_imposed = RegularGridInterpolator(
            (np.linspace(domain_min_y, domain_max_y, y_dim), np.linspace(domain_min_x, domain_max_x, x_dim)),
            self.imposed_vy_grid, bounds_error=False, fill_value=0
        )

        particle_coords = np.column_stack((fluid_array.y, fluid_array.x))
        fluid_array.u[:] = vx_interp_func_imposed(particle_coords)
        fluid_array.v[:] = vy_interp_func_imposed(particle_coords)

        vx_interp_func_bact = RegularGridInterpolator(
            (np.linspace(domain_min_y, domain_max_y, y_dim), np.linspace(domain_min_x, domain_max_x, x_dim)),
            fluid_array.u.reshape(y_dim, x_dim),
            bounds_error=False, fill_value=0
        )
        vy_interp_func_bact = RegularGridInterpolator(
            (np.linspace(domain_min_y, domain_max_y, y_dim), np.linspace(domain_min_x, domain_max_x, x_dim)),
            fluid_array.v.reshape(y_dim, x_dim),
            bounds_error=False, fill_value=0
        )

        v_particles_y = vy_interp_func_bact(self.bact_positions)
        v_particles_x = vx_interp_func_bact(self.bact_positions)

        dt = solver.dt
        self.bact_positions[:, 0] += v_particles_y * dt
        self.bact_positions[:, 1] += v_particles_x * dt

        self.bact_positions[:, 0] = (self.bact_positions[:, 0] - domain_min_y) % (domain_max_y - domain_min_y) + domain_min_y
        self.bact_positions[:, 1] = (self.bact_positions[:, 1] - domain_min_x) % (domain_max_x - domain_min_x) + domain_min_x

        if solver.count % trajectory_store_interval == 0:
            for i_bact in range(self.bact_positions.shape[0]):
                self.all_bact_trajectories[i_bact].append(self.bact_positions[i_bact, :].copy())

        rho_b_from_particles_grid = np.zeros((y_dim, x_dim), dtype=float)
        bact_y_indices = np.clip(np.floor((self.bact_positions[:, 0] - domain_min_y) / dy).astype(int), 0, y_dim - 1)
        bact_x_indices = np.clip(np.floor((self.bact_positions[:, 1] - domain_min_x) / dx).astype(int), 0, x_dim - 1)

        particle_mass_contribution_factor = 0.6
        particle_mass_contribution = rho_max * particle_mass_contribution_factor / (dx*dy)

        for i in range(len(self.bact_positions)):
            rho_b_from_particles_grid[bact_y_indices[i], bact_x_indices[i]] += particle_mass_contribution

        fluid_array.rho_b_grown[:] = np.clip(rho_b_from_particles_grid.ravel()[
            np.clip(np.floor((fluid_array.y - domain_min_y) / dy).astype(int), 0, y_dim - 1) * x_dim +
            np.clip(np.floor((fluid_array.x - domain_min_x) / dx).astype(int), 0, x_dim - 1)
        ], 0, rho_max)

        if solver.count % solver.set_print_freq == 0 or solver.count == int(total_sim_time / dt_global) - 1:
            print(f"  Progresso: {((solver.count + 1) / (total_sim_time / dt_global)) * 100:.1f}% completo. "
                    f"Max c_s: {np.max(fluid_array.c_s):.2e}, Max rho_b (partículas): {np.max(fluid_array.rho_b_grown):.2f}, "
                    f"Max u: {np.max(np.abs(fluid_array.u)):.2e}, Max v: {np.max(np.abs(fluid_array.v)):.2e}")

    def post_process(self, info):
        print("Simulação concluída. Gerando plot final...")

        final_trajectories_for_plot = [np.array(traj_list) for traj_list in self.all_bact_trajectories if traj_list]

        fig, ax = plt.subplots(figsize=(8, 8))

        fluid_array = self.particles[0]
        x_fluid = fluid_array.x
        y_fluid = fluid_array.y
        c_s_fluid = fluid_array.c_s

        x_grid_plot = np.linspace(domain_min_x, domain_max_x, x_dim)
        y_grid_plot = np.linspace(domain_min_y, domain_max_y, y_dim)
        _X_plot, _Y_plot = np.meshgrid(x_grid_plot, y_grid_plot)

        cs_grid = np.zeros((y_dim, x_dim))
        count_cs_grid = np.zeros((y_dim, x_dim))

        for i in range(len(x_fluid)):
            gx = int((x_fluid[i] - domain_min_x) / dx)
            gy = int((y_fluid[i] - domain_min_y) / dy)
            if 0 <= gx < x_dim and 0 <= gy < y_dim:
                cs_grid[gy, gx] += c_s_fluid[i]
                count_cs_grid[gy, gx] += 1

        cs_grid[count_cs_grid > 0] /= count_cs_grid[count_cs_grid > 0]
        cs_grid = np.ma.masked_where(count_cs_grid == 0, cs_grid)

        cs_min_plot = np.percentile(c_s_fluid[c_s_fluid > 1e-9], 1) if np.any(c_s_fluid > 1e-9) else 1e-7
        cs_max_plot = np.max(c_s_fluid)
        if cs_max_plot <= cs_min_plot: cs_max_plot = cs_min_plot * 1.1 + 1e-9

        levels_cs_contourf = np.logspace(np.log10(cs_min_plot), np.log10(cs_max_plot), 30)
        norm_cs = LogNorm(vmin=cs_min_plot, vmax=cs_max_plot)

        cont_f = ax.contourf(_X_plot, _Y_plot, cs_grid, levels=levels_cs_contourf, cmap='viridis', alpha=0.6, norm=norm_cs, zorder=1)
        cbar = fig.colorbar(cont_f, ax=ax, label='Concentração de Surfactante (c_s)', fraction=0.046, pad=0.04)

        for traj_array in final_trajectories_for_plot:
            if traj_array.shape[0] > 1:
                ax.plot(traj_array[:, 1], traj_array[:, 0], '-', color='darkorange', alpha=0.5, linewidth=0.7, zorder=3)

        ax.set_title(f"Movimento Bacteriano Dinâmico (T={total_sim_time:.2f}, dt={dt_global:.1e})")
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(True, linestyle=':', alpha=0.4, zorder=0)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(domain_min_x, domain_max_x)
        ax.set_ylim(domain_min_y, domain_max_y)
        plt.tight_layout()

        output_plot_dir = os.path.join(info.output_directory, 'plots')
        os.makedirs(output_plot_dir, exist_ok=True)
        output_path = os.path.join(output_plot_dir, 'simulacao_final.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')

        plt.show()

if __name__ == '__main__':
    app = BiomassSurfactantApp()
    app.run()