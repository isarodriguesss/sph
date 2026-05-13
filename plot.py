import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysph.solver.utils import load
import glob
import os

os.makedirs("main_output/movie", exist_ok=True)
files = sorted(glob.glob("main_output/main_*.hdf5"))
n = len(files)
# Key frames: initial, early expansion, peak activity, decline, final
indices = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 90, 110, 130, 150, 180, 210, n - 1]
indices = [i for i in indices if i < n]

for i in indices:
    data = load(files[i])
    fluid = data["arrays"]["fluid"]
    x, y = fluid.x, fluid.y
    rho_b = fluid.rho_b_grown
    cs = fluid.cs
    c_o = fluid.c_o
    au_flag = fluid.au_flag

    fig, axes = plt.subplots(1, 4, figsize=(32, 8))

    # 1: rho_b (biomassa)
    sc1 = axes[0].scatter(x, y, c=rho_b, cmap="viridis", s=1.5, vmin=0, vmax=1)
    axes[0].set_xlim(-3, 3)
    axes[0].set_ylim(-3, 3)
    axes[0].set_aspect("equal")
    axes[0].set_title(f"rho_b — {os.path.basename(files[i])}")
    plt.colorbar(sc1, ax=axes[0])

    # 2: cs (surfactante)
    sc2 = axes[1].scatter(
        x, y, c=cs, cmap="hot", s=1.5, vmin=0, vmax=max(np.percentile(cs, 99.5), 0.1)
    )
    axes[1].set_xlim(-3, 3)
    axes[1].set_ylim(-3, 3)
    axes[1].set_aspect("equal")
    axes[1].set_title(f"cs (surfactant) — {os.path.basename(files[i])}")
    plt.colorbar(sc2, ax=axes[1])

    # 3: c_o (Pass M-A — osmolito / pressao osmotica via van't Hoff)
    # Pi = i·C·R·T → c_o e proporcional a pressao osmotica.
    # Esperado: alto nas baias (agar entre dendritos), ~0 nas pontas (agar virgem).
    sc3 = axes[2].scatter(x, y, c=c_o, cmap="cividis", s=1.5, vmin=0, vmax=1)
    axes[2].set_xlim(-3, 3)
    axes[2].set_ylim(-3, 3)
    axes[2].set_aspect("equal")
    axes[2].set_title(f"c_o (osmolito ∝ Π) — {os.path.basename(files[i])}")
    plt.colorbar(sc3, ax=axes[2])

    # 4: au_flag (motilidade flagelar)
    # vmax=3.0 = f0 (K.22), para ver diretamente se o gate esta saturando
    sc4 = axes[3].scatter(x, y, c=au_flag, cmap="plasma", s=1.5, vmin=0, vmax=3.0)
    axes[3].set_xlim(-3, 3)
    axes[3].set_ylim(-3, 3)
    axes[3].set_aspect("equal")
    axes[3].set_title(f"|a_flag| — {os.path.basename(files[i])}")
    plt.colorbar(sc4, ax=axes[3])

    plt.tight_layout()
    plt.savefig(f"main_output/movie/frame_i7_{i:03d}.png", dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Frame {i:03d} — {os.path.basename(files[i])}")

print(f"Done. {n} total files.")
