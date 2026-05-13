import numpy as np
from pysph.base.utils import get_particle_array


def create_initial_state(
    x_dim=128,
    y_dim=128,
    rho_max=1.0,
    dt=0.001,
    x_min=-3.0,
    x_max=3.0,
    y_min=-3.0,
    y_max=3.0,
):
    x = np.linspace(x_min, x_max, x_dim)
    y = np.linspace(y_min, y_max, y_dim)
    X_grid, Y_grid = np.meshgrid(x, y)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    x_part = X_grid.ravel()
    y_part = Y_grid.ravel()
    z_part = np.zeros_like(x_part)
    m_part = np.ones_like(x_part) * dx * dx
    h_part = np.ones_like(x_part) * 1.8 * dx
    rho_part = np.ones_like(x_part) * 1.0
    center_x = (x.min() + x.max()) / 2.0
    center_y = (y.min() + y.max()) / 2.0
    dist = np.sqrt((X_grid - center_x) ** 2 + (Y_grid - center_y) ** 2)
    theta = np.arctan2(Y_grid - center_y, X_grid - center_x)

    R_theta = 0.30 + 0.06 * np.cos(8 * theta)
    rho_b = np.exp(-((dist / R_theta) ** 4))

    noise = 0.10 * np.random.randn(*X_grid.shape)
    rho_b += rho_b * noise
    rho_b = np.clip(rho_b, 0, None)
    rho_b_grown_part = np.clip(rho_b.ravel(), 0, rho_max)

    cs_part = np.clip(rho_b_grown_part * 0.1, 1e-9, None)
    c_o_part = np.zeros_like(x_part)
    c_n_part = np.clip(1.0 - 0.8 * rho_b_grown_part, 1e-9, 1.0)

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
        c_o=c_o_part,
        c_n=c_n_part,
        a_rho_b_grown=np.zeros_like(x_part),
        a_c_s=np.zeros_like(x_part),
        a_c_o=np.zeros_like(x_part),
        a_c_n=np.zeros_like(x_part),
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
    h_solid = np.ones_like(x_solid) * 1.8 * dx
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
