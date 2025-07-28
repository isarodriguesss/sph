import numpy as np
from pysph.base.utils import get_particle_array_wcsph

def create_biomass_surfactant_particles(domain_min_x, domain_max_x, domain_min_y, domain_max_y, h, initial_dx_dy, rho_max, x_dim, y_dim):

    _x = np.linspace(domain_min_x + h/2, domain_max_x - h/2, int((domain_max_x - domain_min_x)/initial_dx_dy[0]))
    _y = np.linspace(domain_min_y + h/2, domain_max_y - h/2, int((domain_max_y - domain_min_y)/initial_dx_dy[1]))
    _X, _Y = np.meshgrid(_x, _y)

    x = _X.ravel()
    y = _Y.ravel()
    z = np.zeros_like(x)
    u, v, w = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)

    rho = np.ones_like(x) * 1.0
    m = rho * (initial_dx_dy[0] * initial_dx_dy[1])

    cs = np.zeros_like(x) + 1e-9

    """ external_band_start_coord = 1.0
    external_band_width_coord = 0.6
    external_mask = (x >= external_band_start_coord) & (x <= external_band_start_coord + external_band_width_coord)
    cs[external_mask] = 0.3 """

    center_x = (domain_min_x + domain_max_x) / 2.0
    center_y = (domain_min_y + domain_max_y) / 2.0

    _X_grid, _Y_grid = np.meshgrid(np.linspace(domain_min_x, domain_max_x, x_dim), np.linspace(domain_min_y, domain_max_y, y_dim))

    rho_b_grid = np.exp(-((_X_grid - center_x)**2 + (_Y_grid - center_y)**2) / 0.05)
    perturb_offset_x = initial_dx_dy[0] * 5
    perturb_offset_y = initial_dx_dy[1] * 5
    perturb = 0.5 * np.exp(-((_X_grid - (center_x + perturb_offset_x))**2 + (_Y_grid - (center_y + perturb_offset_y))**2) / 0.01)
    rho_b_grid += perturb
    rho_b_grid = np.clip(rho_b_grid, 0, rho_max)

    rho_b_grown = rho_b_grid.ravel()

    fluid_properties = ['cs', 'rho', 'm', 'h', 'u', 'v', 'w', 'x', 'y', 'z', 'rho_b_grown']
    fluid = get_particle_array_wcsph(name='fluid', x=x, y=y, z=z, u=u, v=v, w=w, rho=rho, m=m, h=h, c_s=cs, rho_b_grown=rho_b_grown)

    return [fluid]