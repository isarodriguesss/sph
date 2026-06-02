import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysph.solver.utils import load
import glob
import os

os.makedirs("main_output/movie", exist_ok=True)
all_files = sorted(glob.glob("main_output/main_*.hdf5"))

# Subamostra: nao gera frame para TODO main_*.hdf5. Runs com Pass N tem dt
# pequeno → milhares de snapshots; plotar todos e desnecessario e lento.
# Mantem ~TARGET_FRAMES distribuidos uniformemente + sempre o ultimo (estado final).
TARGET_FRAMES = 40
stride = max(1, len(all_files) // TARGET_FRAMES)
files = all_files[::stride]
if all_files and all_files[-1] not in files:
    files.append(all_files[-1])
n = len(files)

rho0 = 1.0  # densidade de referência SPH

for i, fpath in enumerate(files):
    data = load(fpath)
    fluid = data["arrays"]["fluid"]
    x, y = fluid.x, fluid.y
    rho_b = fluid.rho_b_grown
    cs = fluid.cs
    rho = fluid.rho  # densidade SPH — detecta gaps estruturais

    # Máscara colônia: só partículas com biomassa
    colony = rho_b > 0.05

    fig, axes = plt.subplots(1, 3, figsize=(22, 8))

    # --- Painel 1: rho_b (biomassa) ---
    sc1 = axes[0].scatter(x, y, c=rho_b, cmap="viridis", s=2.0, vmin=0, vmax=1)
    axes[0].set_xlim(-3, 3)
    axes[0].set_ylim(-3, 3)
    axes[0].set_aspect("equal")
    axes[0].set_title(f"rho_b — {os.path.basename(fpath)}")
    plt.colorbar(sc1, ax=axes[0])

    # --- Painel 2: cs (surfactante) ---
    sc2 = axes[1].scatter(
        x, y, c=cs, cmap="hot", s=2.0, vmin=0, vmax=max(np.percentile(cs, 99.5), 0.1)
    )
    axes[1].set_xlim(-3, 3)
    axes[1].set_ylim(-3, 3)
    axes[1].set_aspect("equal")
    axes[1].set_title(f"cs (surfactant) — {os.path.basename(fpath)}")
    plt.colorbar(sc2, ax=axes[1])

    # --- Painel 3: rho SPH normalizado (gaps estruturais) ---
    # Apenas partículas da colônia — no agar rho naturalmente baixo
    rho_norm = rho / rho0
    # Fundo: agar em cinza claro
    axes[2].scatter(x[~colony], y[~colony], c="0.88", s=1.0, linewidths=0)
    # Colônia: colorido por rho/rho0 — vermelho = gap, verde = cheio
    sc3 = axes[2].scatter(
        x[colony],
        y[colony],
        c=rho_norm[colony],
        cmap="RdYlGn",  # vermelho (gap) → amarelo → verde (cheio)
        s=3.0,
        vmin=0.4,
        vmax=1.2,
    )
    axes[2].set_xlim(-3, 3)
    axes[2].set_ylim(-3, 3)
    axes[2].set_aspect("equal")
    axes[2].set_title(f"rho/rho0 (gaps: vermelho) — {os.path.basename(fpath)}")
    cb3 = plt.colorbar(sc3, ax=axes[2])
    cb3.set_label("rho / rho0   [verde=cheio, vermelho=gap]")

    # Contorno dos gaps críticos (rho < 0.7 * rho0)
    gap_mask = colony & (rho_norm < 0.7)
    if gap_mask.any():
        axes[2].scatter(
            x[gap_mask],
            y[gap_mask],
            facecolors="none",
            edgecolors="blue",
            s=12,
            linewidths=0.5,
            label=f"gap crítico ({gap_mask.sum()})",
        )
        axes[2].legend(fontsize=7, loc="upper right")

    plt.tight_layout()
    out = f"main_output/movie/frame_i7_{i:03d}.png"
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Frame {i:03d} — {os.path.basename(fpath)}")

print(f"Done. {n} total files.")
