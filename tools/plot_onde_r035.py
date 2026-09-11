import argparse
import glob
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

R_CLIFF = 0.35
R_BAND = 0.65
RHO_B_FLOOR = 1e-6


def fade_att(rho_b):
    """Ramo ATRATIVO (coesao) da BiomassEOS — src/equations.py:828-843."""
    f = np.zeros_like(rho_b)
    mid = (rho_b >= 0.1) & (rho_b < 0.5)
    t = (rho_b[mid] - 0.1) / 0.4
    f[mid] = t * t * (3.0 - 2.0 * t)
    f[rho_b >= 0.5] = 1.0
    return f


def load(path):
    with h5py.File(path, "r") as h:
        t = float(h["solver_data"].attrs["t"])
        a = h["particles"]["fluid"]["arrays"]
        d = {k: np.asarray(a[k]) for k in ("x", "y", "rho_b_grown")}
    d["t"] = t
    return d


def frames(run):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
    ts = []
    for f in fs:
        with h5py.File(f, "r") as h:
            ts.append(float(h["solver_data"].attrs["t"]))
    return fs, np.array(ts)


def pick(fs, ts, alvo):
    return fs[int(np.argmin(np.abs(ts - alvo)))]


def r99(d):
    r = np.hypot(d["x"], d["y"])
    m = d["rho_b_grown"] > 0.1
    return float(np.percentile(r[m], 99)) if m.sum() > 10 else 1.0


def painel(ax, d, campo, R, cmap, norm, titulo):
    rb = d["rho_b_grown"]
    lim = 1.3 * R
    vis = (np.abs(d["x"]) < lim) & (np.abs(d["y"]) < lim)
    val = campo[vis]
    ordem = np.argsort(val)
    ax.set_facecolor("#1a1a1a")
    ax.scatter(
        d["x"][vis][ordem],
        d["y"][vis][ordem],
        c=val[ordem],
        s=2.0,
        cmap=cmap,
        norm=norm,
        linewidths=0,
    )
    for rr, cor, ls in ((R_CLIFF, "red", "-"), (R_BAND, "cyan", "-"), (R, "white", "--")):
        ax.add_patch(plt.Circle((0, 0), rr, fill=False, ec=cor, ls=ls, lw=1.4))
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(titulo, fontsize=8)
    _ = rb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--rotulos", nargs="+", required=True)
    ap.add_argument("--tempos", nargs="+", type=float, default=[13.0, 29.0, 40.0])
    ap.add_argument("--out", default="plots/onde_r035_serieP.png")
    args = ap.parse_args()

    nt = len(args.tempos)
    fig, axes = plt.subplots(
        len(args.runs), 2 * nt, figsize=(3.05 * 2 * nt, 3.35 * len(args.runs))
    )
    axes = np.atleast_2d(axes)

    norm_rb = LogNorm(vmin=RHO_B_FLOOR, vmax=1.0)
    norm_fd = plt.Normalize(0.0, 1.0)

    for i, (run, rot) in enumerate(zip(args.runs, args.rotulos)):
        fs, ts = frames(run)
        for j, alvo in enumerate(args.tempos):
            d = load(pick(fs, ts, alvo))
            R = r99(d)
            rb = np.clip(d["rho_b_grown"], RHO_B_FLOOR, 1.0)
            painel(
                axes[i, j],
                d,
                rb,
                R,
                "viridis",
                norm_rb,
                f"rho_b  t={d['t']:.0f}  R99={R:.2f}",
            )
            painel(
                axes[i, nt + j],
                d,
                fade_att(d["rho_b_grown"]),
                R,
                "inferno",
                norm_fd,
                f"COESAO (fade)  t={d['t']:.0f}",
            )
        axes[i, 0].set_ylabel(rot, fontsize=11, fontweight="bold")
        axes[i, 0].set_yticks([])

    axes[0, 0].text(
        0.03,
        0.97,
        "vermelho r=0.35 (penhasco)\nciano r=0.65\nbranco R99",
        transform=axes[0, 0].transAxes,
        va="top",
        fontsize=7,
        color="white",
        bbox=dict(fc="#00000099", ec="none", pad=3),
    )
    fig.suptitle(
        "ONDE ESTA O PENHASCO — rho_b (esquerda) x COESAO da EOS (direita).  "
        "Escala do painel = 1.3*R99, entao os circulos fixos encolhem com a colonia.",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(args.out, dpi=110)
    print("escrito:", args.out)


if __name__ == "__main__":
    main()
