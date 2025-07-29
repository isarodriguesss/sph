import numpy as np
from pysph.base.utils import get_particle_array_wcsph

def create_biomass_surfactant_particles(domain_min_x, domain_max_x, domain_min_y, domain_max_y, h, initial_dx_dy, rho_max, x_dim, y_dim):

    #* Cria uma grade de pontos no domínio 3D definido
    _x = np.linspace(domain_min_x + h/2, domain_max_x - h/2, int((domain_max_x - domain_min_x)/initial_dx_dy[0]))
    _y = np.linspace(domain_min_y + h/2, domain_max_y - h/2, int((domain_max_y - domain_min_y)/initial_dx_dy[1]))
    _X, _Y = np.meshgrid(_x, _y)

    #* Coordenadas das partículas
    x = _X.ravel()
    y = _Y.ravel()
    z = np.zeros_like(x)

    #* Velocidades iniciais das partículas
    u, v, w = np.zeros_like(x), np.zeros_like(x), np.zeros_like(x)

    #* Densidade inicial das partículas
    rho = np.ones_like(x) * 1.0
    #* Massa das partículas
    m = rho * (initial_dx_dy[0] * initial_dx_dy[1])
    #* Concentração inicial de surfactante
    cs = np.zeros_like(x) + 1e-9

    center_x = (domain_min_x + domain_max_x) / 2.0
    center_y = (domain_min_y + domain_max_y) / 2.0

    _X_grid, _Y_grid = np.meshgrid(np.linspace(domain_min_x, domain_max_x, x_dim), np.linspace(domain_min_y, domain_max_y, y_dim))

    #* Perturbação inicial da densidade de biomassa
    gaussian1 = np.exp(-((_X_grid - center_x)**2 + (_Y_grid - center_y)**2) / 0.05)
    perturb_offset_x = initial_dx_dy[0] * 5
    perturb_offset_y = initial_dx_dy[1] * 5
    gaussian2 = np.exp(-((_X_grid - (center_x + perturb_offset_x))**2 + (_Y_grid - (center_y + perturb_offset_y))**2) / 0.01)

    # DEBUG PRINTS PARA ACOMPANHAR OS VALORES
    print(f"DEBUG particles.py: Max da 1a gaussiana (sem fator): {np.max(gaussian1):.4f}")
    print(f"DEBUG particles.py: Max da 2a gaussiana (sem fator): {np.max(gaussian2):.4f}")

    # Ajustando os fatores. Ex: 0.6 para o primeiro pico e 0.2 para a perturbação
    rho_b_grid = 0.6 * gaussian1 # Reduzir o pico inicial para 0.6
    print(f"DEBUG particles.py: Max rho_b_grid após 1a gaussiana (fator 0.6): {np.max(rho_b_grid):.4f}")

    perturb = 0.2 * gaussian2 # Reduzir o pico da perturbação para 0.2
    print(f"DEBUG particles.py: Max da perturbação (fator 0.2): {np.max(perturb):.4f}")

    rho_b_grid += perturb
    print(f"DEBUG particles.py: Max rho_b_grid após somar perturbação (soma de 0.6 e 0.2, esperada 0.8 ou menos): {np.max(rho_b_grid):.4f}")

    # Remova ou comente esta linha para eliminar a contribuição do ruído aleatório
    # rho_b_grid += 0.15 * np.random.randn(*_X.shape)
    # Se você quiser ruído, adicione-o com um fator muito pequeno para não saturar.
    # Ex: rho_b_grid += 0.01 * np.random.randn(*_X.shape)

    rho_b_grid = np.clip(rho_b_grid, 0, rho_max)
    print(f"DEBUG particles.py: Max rho_b_grid após clipagem: {np.max(rho_b_grid):.4f}")
    # --- FIM DA CORREÇÃO ---

    rho_b_grown = rho_b_grid.ravel()

    #* Cria o objeto ParticleArray com as propriedades definidas
    fluid = get_particle_array_wcsph(name='fluid', x=x, y=y, z=z, u=u, v=v, w=w, rho=rho, m=m, h=h, c_s=cs, rho_b_grown=rho_b_grown)

    return [fluid]