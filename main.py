import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
import os

from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from scipy.interpolate import RegularGridInterpolator, griddata

from src.particles import create_biomass_surfactant_particles
from src.scheme import MyBiomassScheme

# Parâmetros da Simulação
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
# Ajuste a frequência de print e armazenamento para algo razoável
print_freq = 500
trajectory_store_interval = 10

n_bact = 300
r0_bact_initial = dx * 2


class BiomassSurfactantApp(Application):

    def initialize(self):
        # 1. Criação das partículas (inalterado)
        self.particles = create_biomass_surfactant_particles(
            domain_min_x, domain_max_x, domain_min_y, domain_max_y,
            h, (dx, dy), rho_max, x_dim, y_dim,
        )
        
        fluid_pa = self.particles[0]
        
        # 2. Adicionar propriedades de TAXA com a convenção 'a_propriedade'
        if not hasattr(fluid_pa, 'a_rho_b_grown'):
            fluid_pa.add_property('a_rho_b_grown', default=0.0)
        if not hasattr(fluid_pa, 'a_c_s'):
            fluid_pa.add_property('a_c_s', default=0.0)

        # 3. Setup de bactérias e domínio (inalterado)
        center_x = (domain_min_x + domain_max_x) / 2.0
        center_y = (domain_min_y + domain_max_y) / 2.0
        angles = np.linspace(0, 2 * np.pi, n_bact, endpoint=False)
        self.bact_positions = np.array([[center_y + r0_bact_initial * np.sin(a), center_x + r0_bact_initial * np.cos(a)] for a in angles])
        self.all_bact_trajectories = [[] for _ in range(n_bact)]
        self.periodic_domain_min = np.array([domain_min_x, domain_min_y, -h])
        self.periodic_domain_max = np.array([domain_max_x, domain_max_y, h])

        # 4. Setup do campo de fluxo imposto (inalterado)
        self.x_grid_imp, self.y_grid_imp = np.meshgrid(
            np.linspace(domain_min_x, domain_max_x, x_dim),
            np.linspace(domain_min_y, domain_max_y, y_dim)
        )
        self.imposed_vx_grid = -np.sin(self.y_grid_imp * np.pi / (domain_max_y - domain_min_y)) * np.cos(self.x_grid_imp * np.pi / (domain_max_x - domain_min_x)) * 0.1
        self.imposed_vy_grid = np.cos(self.y_grid_imp * np.pi / (domain_max_y - domain_min_y)) * np.sin(self.x_grid_imp * np.pi / (domain_max_x - domain_min_x)) * 0.1

        # 5. **** MUDANÇA CRÍTICA PARA PERFORMANCE ****
        # Crie os interpoladores AQUI, UMA SÓ VEZ, e armazene-os em 'self'.
        print("Criando interpoladores de grade uma única vez...")
        y_coords = np.linspace(domain_min_y, domain_max_y, y_dim)
        x_coords = np.linspace(domain_min_x, domain_max_x, x_dim)
        
        self.vx_interp_func_imposed = RegularGridInterpolator(
            (y_coords, x_coords), self.imposed_vx_grid, 
            bounds_error=False, fill_value=0
        )
        self.vy_interp_func_imposed = RegularGridInterpolator(
            (y_coords, x_coords), self.imposed_vy_grid, 
            bounds_error=False, fill_value=0
        )
        print("Interpoladores criados com sucesso. Iniciando a simulação...")


    def create_particles(self):
        for pa in self.particles:
            if pa.name == 'fluid': 
                pa.add_output_arrays(['rho', 'm', 'x', 'y', 'u', 'v', 'h', 'rho_b_grown', 'c_s', 'ax', 'ay', 'a_rho_b_grown', 'a_c_s'])
        return self.particles


    def create_scheme(self):
        return MyBiomassScheme(fluids=['fluid'], solids=[], dim=2,
                                 rho_max=rho_max, r_growth=r_growth,
                                 sigma=sigma, lambda_=lambda_, beta=beta, gamma=gamma, D=D, mu=mu)
    

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        scheme = self.create_scheme()
        solver = Solver(dim=2, integrator=scheme.get_integrator(), kernel=kernel)
        
        solver.tf = total_sim_time
        solver.dt = dt_global
        solver.set_adaptive_timestep(False)
        solver.set_print_freq(print_freq)
        return solver


    def create_tools(self):
        return []
    

    def post_step(self, solver):
        fluid_array = self.particles[0]

        # **** MUDANÇA CRÍTICA PARA PERFORMANCE ****
        # Obtenha as coordenadas e REUTILIZE os interpoladores pré-construídos.
        # Esta operação agora é extremamente rápida.
        particle_coords = np.column_stack((fluid_array.y, fluid_array.x))
        fluid_array.u[:] = self.vx_interp_func_imposed(particle_coords)
        fluid_array.v[:] = self.vy_interp_func_imposed(particle_coords)
        
        # Mova as bactérias usando o mesmo campo de fluxo imposto de forma eficiente.
        v_particles_y = self.vy_interp_func_imposed(self.bact_positions)
        v_particles_x = self.vx_interp_func_imposed(self.bact_positions)

        dt = solver.dt
        self.bact_positions[:, 0] += v_particles_y * dt
        self.bact_positions[:, 1] += v_particles_x * dt

        # Aplica condições de contorno periódicas às bactérias.
        self.bact_positions[:, 0] = (self.bact_positions[:, 0] - domain_min_y) % (domain_max_y - domain_min_y) + domain_min_y
        self.bact_positions[:, 1] = (self.bact_positions[:, 1] - domain_min_x) % (domain_max_x - domain_min_x) + domain_min_x

        if solver.count % trajectory_store_interval == 0:
            for i_bact in range(self.bact_positions.shape[0]):
                self.all_bact_trajectories[i_bact].append(self.bact_positions[i_bact, :].copy())

        # O print de progresso agora aparecerá em intervalos regulares.
        if solver.count % print_freq == 0:
            print(f"  Progresso: {((solver.count) / (total_sim_time / dt_global)) * 100:.1f}% completo. "
                  f"Max c_s: {np.max(fluid_array.c_s):.2e}, Max rho_b: {np.max(fluid_array.rho_b_grown):.2f}")


    def post_process(self, info):
        print("Simulação concluída. Gerando plot final...")

        final_trajectories_for_plot = [np.array(traj_list) for traj_list in self.all_bact_trajectories if traj_list]
        fig, ax = plt.subplots(figsize=(8, 8))
        fluid_array = self.particles[0]
        x_fluid, y_fluid, c_s_fluid = fluid_array.x, fluid_array.y, fluid_array.c_s

        x_grid_plot = np.linspace(domain_min_x, domain_max_x, x_dim)
        y_grid_plot = np.linspace(domain_min_y, domain_max_y, y_dim)
        _X_plot, _Y_plot = np.meshgrid(x_grid_plot, y_grid_plot)

        # Interpolação para o plot (pode ser simplificada, mas funcional)
        cs_grid = griddata((x_fluid, y_fluid), c_s_fluid, (_X_plot, _Y_plot), method='cubic', fill_value=1e-9)

        cs_min_plot = np.percentile(c_s_fluid[c_s_fluid > 1e-9], 1) if np.any(c_s_fluid > 1e-9) else 1e-7
        cs_max_plot = np.max(c_s_fluid)
        if cs_max_plot <= cs_min_plot: cs_max_plot = cs_min_plot * 1.1 + 1e-9

        levels_cs_contourf = np.logspace(np.log10(cs_min_plot), np.log10(cs_max_plot), 30)
        norm_cs = LogNorm(vmin=cs_min_plot, vmax=cs_max_plot)

        cont_f = ax.contourf(_X_plot, _Y_plot, cs_grid, levels=levels_cs_contourf, cmap='viridis', alpha=0.6, norm=norm_cs, zorder=1)
        fig.colorbar(cont_f, ax=ax, label='Concentração de Surfactante (c_s)', fraction=0.046, pad=0.04)

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

        output_plot_dir = os.path.join(self.output_dir, 'plots')
        os.makedirs(output_plot_dir, exist_ok=True)
        output_path = os.path.join(output_plot_dir, 'simulacao_final.png')
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot final salvo em: {output_path}")
        plt.show()


if __name__ == '__main__':
    app = BiomassSurfactantApp()
    app.run()