import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysph.solver.utils import load
import glob
import os

# Configuracoes
OUTPUT_DIR = "main_output"
PLOT_FILE = os.path.join(OUTPUT_DIR, "trajectories.png")
MAX_TRAJECTORIES = 1500  # Limite de linhas para nao poluir o grafico
MOVE_THRESHOLD = 0.2  # Só plota partículas que andaram mais que esta distância


def main():
    files = sorted(glob.glob(os.path.join(OUTPUT_DIR, "main_*.hdf5")))
    if not files:
        print(f"Nenhum arquivo HDF5 encontrado em {OUTPUT_DIR}/")
        return

    print(f"Encontrados {len(files)} arquivos. Carregando dados de trajetória...")

    X_hist = []
    Y_hist = []

    # 1. Carregar todas as posicoes ao longo do tempo
    for i, f in enumerate(files):
        if i % 20 == 0:
            print(f"Lendo frame {i}/{len(files)}...")
        data = load(f)
        fluid = data["arrays"]["fluid"]
        X_hist.append(fluid.x.copy())
        Y_hist.append(fluid.y.copy())

        # Guardar as propriedades do ultimo frame para o background
        if i == len(files) - 1:
            final_x = fluid.x
            final_y = fluid.y
            final_rho_b = fluid.rho_b_grown

    # Converter para arrays numpy: shape (n_frames, n_particles)
    X_hist = np.array(X_hist)
    Y_hist = np.array(Y_hist)
    n_particles = X_hist.shape[1]

    # 2. Filtrar particulas
    print("Calculando deslocamentos...")
    # Distancia total percorrida entre o inicio e o fim
    dist_moved = np.sqrt((X_hist[-1] - X_hist[0]) ** 2 + (Y_hist[-1] - Y_hist[0]) ** 2)

    # Selecionar apenas particulas que sairam do lugar (exclui o nucleo ancorado)
    moving_mask = dist_moved > MOVE_THRESHOLD
    selected_indices = np.where(moving_mask)[0]

    # Se houver muitas particulas, pegamos uma amostra aleatoria para o plot ficar legivel
    if len(selected_indices) > MAX_TRAJECTORIES:
        np.random.seed(42)  # Semente fixa para reprodutibilidade
        selected_indices = np.random.choice(
            selected_indices, MAX_TRAJECTORIES, replace=False
        )
        print(
            f"Amostrando {MAX_TRAJECTORIES} trajetórias de {np.sum(moving_mask)} ativas."
        )
    else:
        print(f"Plotando todas as {len(selected_indices)} trajetórias ativas.")

    # 3. Gerar o Plot
    print("Gerando a imagem...")
    fig, ax = plt.subplots(figsize=(12, 12), facecolor="black")
    ax.set_facecolor("black")

    # Background: plota a biomassa final (rho_b) com cores escuras para dar contraste
    sc = ax.scatter(
        final_x, final_y, c=final_rho_b, cmap="Purples", s=2.0, alpha=0.3, zorder=1
    )

    # Plotar as trajetorias (linhas)
    # Usando alpha baixo (transparencia) os caminhos que se sobrepoem ficam mais brilhantes
    for idx in selected_indices:
        ax.plot(
            X_hist[:, idx],
            Y_hist[:, idx],
            color="cyan",
            alpha=0.25,
            linewidth=0.8,
            zorder=2,
        )

        # Colocar um ponto vermelho indicando onde a particula parou no final
        ax.scatter(X_hist[-1, idx], Y_hist[-1, idx], color="red", s=4.0, zorder=3)

    # Marcador central (origem)
    ax.plot(0, 0, marker="+", color="white", markersize=15, zorder=4)

    # Estilizacao
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_aspect("equal")
    ax.set_title(
        "Trajetórias Lagrangianas das Bactérias (Swarmers)", color="white", fontsize=16
    )
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_color("gray")

    # Salvar
    plt.tight_layout()
    plt.savefig(PLOT_FILE, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Sucesso! Imagem salva em: {PLOT_FILE}")


if __name__ == "__main__":
    main()
