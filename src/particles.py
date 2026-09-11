import numpy as np
from pysph.base.utils import get_particle_array


SEED_AMP = 0.06  # P19 testou 0.0 (sem semente) e REPROVOU — ver licao #83
SEED_MODE = 8  # M5 (N=16) nao se sustentou: dedos finos mas ocos, agar cercado 417 vs 181

# SEED_FRAC e a modulacao RELATIVA do raio. A forma antiga (0.30 + 0.06*cos) e
# identica a 0.30*(1 + 0.20*cos), entao o modo "quartica" fica bit-identico. Expressar
# como fracao e o que impede a semente de enfraquecer sozinha quando R_base muda.
SEED_FRAC = SEED_AMP / 0.30

INOC_MODE = "plato"  # "quartica" = perfil historico (preservado, §10)
INOC_A = 0.45  # plato: rho_b uniforme. 0.45 e o canto de TRES restricoes simultaneas —
INOC_R = 0.500  # pin exige rho_b>0.40 (c_n_eq<0.6), fade_rep pica em ~0.5, divisao
INOC_W = 0.145  # rapida quer <0.35. Em 0.45: c_n_eq=0.571, fade_rep=0.957, T_div=12s.
# R=0.500 conserva a biomassa (1.02x) — o I1 reprovou por +150% de massa, nao pela borda.
# W = 1.5h: a borda tem de ser resolvida pelo kernel (Liu §3.3).


def create_initial_state(
    x_dim=128,
    y_dim=128,
    rho_max=1.0,
    dt=0.001,
    x_min=-3.0,
    x_max=3.0,
    y_min=-3.0,
    y_max=3.0,
    seed=None,
):
    if seed is not None:  # reprodutibilidade: mesma condicao inicial entre rotas
        np.random.seed(seed)
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

    # A 4a potencia faz um quase top-hat: rho_b cai de 1.0 a 0.15 entre r=0.3 e 0.45
    # e a ZERO em r=0.6. Isso deixa a borda do nucleo sem biomassa desde t=0 e nada
    # a repoe depois. Suavizar (I1: p=2, R=0.45) CURA a borda mas mata o motor:
    # +2.5x biomassa -> cs mais uniforme -> a_mar 2.30->0.75 e morfologia vira
    # Circular (licao #48). O perfil agudo e o preco da seletividade dendritica.
    if INOC_MODE == "plato":
        # Densidade uniforme com linha de contato — a gota depositada nao tem cauda
        # gaussiana. A quartica poe 1/3 da biomassa numa saia de densidade decrescente
        # e o resto numa cauda de 1e-3 que nao divide (gate 0.05), nao produz cs e tem
        # fade=0. Ver licao #83 e a analise preditiva do inoculo.
        R_theta = INOC_R * (1.0 + SEED_FRAC * np.cos(SEED_MODE * theta))
        t_edge = np.clip((R_theta - dist) / INOC_W, 0.0, 1.0)
        rho_b = INOC_A * t_edge * t_edge * (3.0 - 2.0 * t_edge)
    else:
        R_theta = 0.30 * (1.0 + SEED_FRAC * np.cos(SEED_MODE * theta))
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
