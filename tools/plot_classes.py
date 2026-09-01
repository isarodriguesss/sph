"""Plota a CLASSE de cada particula, nao o campo.

Vazio (sem particula) e agar (`rho_b=0`) renderizam com a mesma cor escura em
qualquer painel de campo — foi o que inverteu o diagnostico duas vezes (licao #66).
Este painel os separa por construcao.

    python tools/plot_classes.py [dir_do_run] [saida.png]
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from pysph.solver.utils import load

CMAP = ListedColormap(["#ff00ff", "#00ff66", "#ffd400", "#ff3b1f"])
LEGENDA = (
    "MAGENTA = AGAR (rho_b=0)\nVERDE = VIVA >0.1\n"
    "AMARELO = LIMBO\nVERMELHO = FILLER\nPRETO = VAZIO"
)
DX = 0.0538


def classes(pa):
    rb = pa.rho_b_grown
    fil = pa.is_filler > 0.5
    return np.where(fil, 3, np.where(rb >= 0.1, 1, np.where(rb >= 1e-12, 2, 0)))


def main(run="runs/C4", out="plots/classes.png", alvos=(0, 4.4, 8.7, 12.9, 21.5, 29.1, 38.3, 50.0)):
    fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
    ts = [load(f)["solver_data"]["t"] for f in fs]
    ks = [int(np.argmin([abs(t - a) for t in ts])) for a in alvos]

    fig, axes = plt.subplots(2, len(ks), figsize=(3.1 * len(ks), 6.6))
    for c, k in enumerate(ks):
        d = load(fs[k])
        pa = d["arrays"]["fluid"]
        t = d["solver_data"]["t"]
        x, y, rb = pa.x, pa.y, pa.rho_b_grown
        r = np.hypot(x, y)
        R99 = np.percentile(r[rb > 0.1], 99)
        col_all = classes(pa)
        for row, L in ((0, 1.1), (1, max(1.3 * R99, 1.2))):
            s = (np.abs(x) < L) & (np.abs(y) < L)
            ax = axes[row, c]
            ax.scatter(x[s], y[s], c=col_all[s], s=(34 if row == 0 else 3.2),
                       cmap=CMAP, vmin=-0.5, vmax=3.5, linewidths=0)
            for rad, cc in ((0.35, "red"), (0.65, "cyan")):
                ax.add_patch(plt.Circle((0, 0), rad, fill=False, ec=cc, lw=1.2, ls="--"))
            if row == 1:
                ax.add_patch(plt.Circle((0, 0), R99, fill=False, ec="w", lw=1.0, ls="--"))
            ax.set_aspect("equal")
            ax.set_xlim(-L, L)
            ax.set_ylim(-L, L)
            ax.set_facecolor("#000")
            ax.tick_params(labelsize=6)
            ax.set_title(
                f"t={t:.1f}" + ("  (zoom juncao)" if row == 0 else f"  R99={R99:.2f}"),
                fontsize=8.5,
            )
    axes[0, 0].text(0.02, 0.98, LEGENDA, transform=axes[0, 0].transAxes, va="top",
                    fontsize=7.5, color="w", bbox=dict(fc="#000", ec="w", alpha=0.9))
    plt.tight_layout()
    plt.savefig(out, dpi=118)
    print(out)

    print("\nAGAR ENGOLIDO (rho_b=0, is_filler=0, dentro do R99)")
    print(f"{'t':>6} {'R99':>5} {'n agar':>8} {'% area':>7}")
    for k in ks:
        d = load(fs[k])
        pa = d["arrays"]["fluid"]
        t = d["solver_data"]["t"]
        r = np.hypot(pa.x, pa.y)
        rb = pa.rho_b_grown
        R99 = np.percentile(r[rb > 0.1], 99)
        ins = (r < R99) & (pa.is_filler < 0.5) & (rb < 1e-12)
        pct = 100 * ins.sum() * DX * DX / (np.pi * R99 * R99)
        print(f"{t:6.1f} {R99:5.2f} {int(ins.sum()):8d} {pct:6.0f}%")


if __name__ == "__main__":
    main(*(sys.argv[1:3] or ["runs/C4", "plots/classes.png"]))
