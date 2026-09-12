"""Painel por CLASSE em R99 igual: viva, filler, absorvido (phi_osm), limbo, agar.

Vazio e agar renderizam iguais em painel de campo (licao #66); aqui cada classe tem cor.

    python tools/plot_osm.py --raios 1.6 2.4 3.2 --saida plots/osm_classes.png RUN...
"""
import argparse
import glob
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CORES = {"agar": "#e9e4f0", "limbo": "#f2c14e", "absorvido": "#4aa3df",
         "filler": "#d1495b", "viva": "#2a9d3f"}


def frame(run, Ra):
    for f in sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5"))):
        with h5py.File(f, "r") as h:
            a = h["particles"]["fluid"]["arrays"]
            rb = np.asarray(a["rho_b_grown"], float)
            x, y = np.asarray(a["x"], float), np.asarray(a["y"], float)
            port = rb > 0.1
            if port.sum() < 20 or np.percentile(np.hypot(x, y)[port], 99) < Ra:
                continue
            t = float(h["solver_data"].attrs["t"])
            isf = np.asarray(a["is_filler"], float) > 0.5
            osm = np.asarray(a["phi_osm"], float) > 0 if "phi_osm" in a else np.zeros(len(x), bool)
            return t, x, y, rb, isf, osm
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--raios", nargs="+", type=float, default=[1.6, 2.4, 3.2])
    ap.add_argument("--saida", default="plots/osm_classes.png")
    args = ap.parse_args()
    nr, nc = len(args.runs), len(args.raios)
    fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 4.2 * nr), squeeze=False)
    for i, run in enumerate(args.runs):
        for j, Ra in enumerate(args.raios):
            ax = axs[i, j]
            ax.set_xticks([])
            ax.set_yticks([])
            F = frame(run, Ra)
            if F is None:
                ax.set_title(f"{os.path.basename(run)}: nao chegou a R99={Ra}", fontsize=9)
                continue
            t, x, y, rb, isf, osm = F
            L = 1.25 * Ra
            w = (np.abs(x) < L) & (np.abs(y) < L)
            cls = {"agar": (rb < 1e-12) & ~isf & ~osm,
                   "limbo": (rb >= 1e-12) & (rb < 0.1) & ~isf & ~osm,
                   "absorvido": osm & (rb < 0.1),
                   "filler": isf,
                   "viva": (rb >= 0.1) & ~isf}
            s = 2.2 * (1.25 * 2.4 / L) ** 2
            for k in ("agar", "limbo", "absorvido", "filler", "viva"):
                m = cls[k] & w
                ax.scatter(x[m], y[m], s=s, c=CORES[k], lw=0, label=k)
            ax.set_xlim(-L, L)
            ax.set_ylim(-L, L)
            ax.set_aspect("equal")
            ax.set_facecolor("black")
            ax.set_title(f"{os.path.basename(run)}  R99={Ra}  t={t:.1f}", fontsize=9)
    axs[0, 0].legend(loc="lower left", fontsize=7, markerscale=3)
    fig.tight_layout()
    fig.savefig(args.saida, dpi=110)
    print(args.saida)


if __name__ == "__main__":
    main()
