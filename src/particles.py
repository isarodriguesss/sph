import numpy as np
from pysph.base.utils import get_particle_array


def create_initial_state(x_dim=128, y_dim=128, rho_max=1.0, dt=0.001):
    x = np.linspace(-1, 5, x_dim)
    y = np.linspace(-1, 5, y_dim)
    X_grid, Y_grid = np.meshgrid(x, y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    x_part = X_grid.ravel()
    y_part = Y_grid.ravel()
    z_part = np.zeros_like(x_part)
    m_part = np.ones_like(x_part) * dx * dx
    h_part = np.ones_like(x_part) * 1.2 * dx
    rho_part = np.ones_like(x_part) * 1.0

    center_x = (x.min() + x.max()) / 2.0
    center_y = (y.min() + y.max()) / 2.0
    """ dist_from_center = np.sqrt((X_grid - center_x) ** 2 + (Y_grid - center_y) ** 2)
    initial_rho_b_radius = 0.25
    seed_mask = dist_from_center < initial_rho_b_radius """
    rho_b = np.exp(-((X_grid - center_x) ** 2 + (Y_grid - center_y) ** 2) / 0.05)
    perturb_offset_x = dx * 5
    perturb_offset_y = dy * 5
    perturb = 0.5 * np.exp(
        -(
            (X_grid - (center_x + perturb_offset_x)) ** 2
            + (Y_grid - (center_y + perturb_offset_y)) ** 2
        )
        / 0.01
    )
    rho_b += perturb
    noise = 0.15 * np.random.randn(*X_grid.shape)
    rho_b += rho_b * noise

    rho_b = np.clip(rho_b, 0, None)

    """ max_initial_rho_b = np.max(rho_b)
    if max_initial_rho_b > 0:
        rho_b = (rho_b / max_initial_rho_b) * 0.2 * rho_max """

    rho_b_grown_part = np.clip(rho_b.ravel(), 0, rho_max)

    cs_part = np.ones_like(x_part) * 1e-9

    fluid = get_particle_array(
        name="fluid",
        x=x_part,
        y=y_part,
        z=z_part,
        m=m_part,
        h=h_part,
        ah=np.zeros_like(x_part),
        rho=rho_part,
        rho_b_grown=rho_b_grown_part,
        cs=cs_part,
        a_rho_b_grown=np.zeros_like(x_part),
        a_c_s=np.zeros_like(x_part),
        m0=m_part.copy(),
        am=np.zeros_like(x_part),
    )

    num_layers = 2
    spacing = dx

    x_left = np.arange(x.min() - num_layers * spacing, x.min(), spacing)
    x_right = np.arange(
        x.max() + spacing, x.max() + (num_layers + 1) * spacing, spacing
    )
    y_walls = np.arange(
        y.min() - num_layers * spacing, y.max() + (num_layers + 1) * spacing, spacing
    )

    x_lr, y_lr = np.meshgrid(np.concatenate([x_left, x_right]), y_walls)

    y_top = np.arange(y.max() + spacing, y.max() + (num_layers + 1) * spacing, spacing)
    y_bottom = np.arange(y.min() - num_layers * spacing, y.min(), spacing)
    x_walls = x

    x_tb, y_tb = np.meshgrid(x_walls, np.concatenate([y_bottom, y_top]))

    x_solid = np.concatenate([x_lr.ravel(), x_tb.ravel()])
    y_solid = np.concatenate([y_lr.ravel(), y_tb.ravel()])

    m_solid = np.ones_like(x_solid) * dx * dy
    h_solid = np.ones_like(x_solid) * 1.2 * dx
    rho_solid = np.ones_like(x_solid) * 1.0

    solid = get_particle_array(
        name="solid",
        x=x_solid,
        y=y_solid,
        m=m_solid,
        h=h_solid,
        rho=rho_solid,
        cs=np.zeros_like(x_solid),
    )

    return [fluid, solid]
