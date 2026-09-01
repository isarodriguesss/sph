"""Paineis de `rho_b` + COESAO da EOS, no estilo de `plots/onde_r035.png`.

Duas linhas por rodada: em cima o campo de biomassa (viridis), embaixo o `fade` da
`BiomassEOS` — que e o que decide se a particula transmite esforco. A distincao importa
porque material com `rho_b` entre 0 e 0.1 aparece no painel de biomassa e tem `fade = 0`:
esta la e nao participa de nada mecanico (licao #53). O painel de coesao mostra a colonia
que a FISICA enxerga.

Circulos de referencia: vermelho em `r=0.35` (borda do nucleo denso), ciano em `r=0.65`
(onde a biomassa viva do C4 vai a zero) e branco tracejado em `R99`.

    python tools/plot_coesao.py runs/C4 [saida.png] [t1,t2,...]
    python tools/plot_coesao.py runs/C4 out.png 12.9,29.1,50
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pysph.solver.utils import load

MARCOS = ((0.35, "red", 2.2), (0.65, "cyan", 1.6))
FUNDO = "0.12"


def fade(v):
    """fade_rep = fade_att da BiomassEOS: smoothstep em rho_b [0.1, 0.5]."""
    s = np.zeros_like(v)
    m = (v >= 0.1) & (v < 0.5)
    t = (v[m] - 0.1) / 0.4
    s[m] = t * t * (3.0 - 2.0 * t)
    s[v >= 0.5] = 1.0
    return s


def painel(runs, out="plots/coesao.png", alvos=(12.9, 29.1, 50.0), zoom=2.2):
    nlin = 2 * len(runs)
    fig, axes = plt.subplots(nlin, len(alvos), figsize=(5.5 * len(alvos), 5.5 * nlin),
                             squeeze=False)
    for ri, run in enumerate(runs):
        fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
        ts = [load(f)["solver_data"]["t"] for f in fs]
        nome = run.rstrip("/").split("/")[-1]
        for c, a in enumerate(alvos):
            k = int(np.argmin([abs(t - a) for t in ts]))
            d = load(fs[k])
            pa = d["arrays"]["fluid"]
            t = d["solver_data"]["t"]
            x, y, rb = pa.x, pa.y, pa.rho_b_grown
            r = np.hypot(x, y)
            R99 = float(np.percentile(r[rb > 0.1], 99)) if (rb > 0.1).any() else 1.0
            fd = fade(rb)
            L = zoom * max(R99, 1.0)
            s = r < L
            for sub, campo, ttl, cm in ((0, rb, "rho_b", "viridis"),
                                        (1, fd, "COESAO da EOS (fade)", "inferno")):
                ax = axes[2 * ri + sub, c]
                ax.scatter(x[s], y[s], c=campo[s], s=7, cmap=cm,
                           vmin=0, vmax=1, linewidths=0)
                for rad, col, lw in MARCOS:
                    ax.add_patch(plt.Circle((0, 0), rad, fill=False, ec=col, lw=lw))
                ax.add_patch(plt.Circle((0, 0), R99, fill=False, ec="white",
                                        lw=1.0, ls="--"))
                ax.set_aspect("equal")
                ax.set_xlim(-L, L)
                ax.set_ylim(-L, L)
                ax.set_facecolor(FUNDO)
                ax.set_title(f"{nome} — {ttl} — t={t:.1f}   R99={R99:.2f}",
                             fontsize=10, color="k")
                ax.tick_params(labelsize=7)
    axes[0, 0].text(
        0.02, 0.97,
        "vermelho r=0.35 (nucleo)\nciano r=0.65 (fim da biomassa viva)\nbranco R99",
        transform=axes[0, 0].transAxes, va="top", fontsize=9, color="w",
        bbox=dict(fc="0.15", ec="none", alpha=0.85),
    )
    plt.tight_layout()
    plt.savefig(out, dpi=115)
    print(out)

    print(f"\n{'run':>16} {'t':>6} {'R99':>6} {'n rho_b>0.1':>12} {'fade med':>9} {'% fade>0.5':>11}")
    for run in runs:
        fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
        ts = [load(f)["solver_data"]["t"] for f in fs]
        for a in alvos:
            k = int(np.argmin([abs(t - a) for t in ts]))
            d = load(fs[k])
            pa = d["arrays"]["fluid"]
            rb = pa.rho_b_grown
            r = np.hypot(pa.x, pa.y)
            R99 = float(np.percentile(r[rb > 0.1], 99))
            m = (rb > 0.1) & (r < R99)
            fd = fade(rb[m])
            print(f"{run.rstrip('/').split('/')[-1]:>16} {d['solver_data']['t']:6.1f} "
                  f"{R99:6.2f} {int(m.sum()):12d} {np.median(fd):9.3f} "
                  f"{100 * np.mean(fd > 0.5):10.0f}%")


if __name__ == "__main__":
    args = sys.argv[1:]
    runs = [args[0]] if args else ["runs/C4"]
    out = args[1] if len(args) > 1 else "plots/coesao.png"
    alvos = tuple(float(v) for v in args[2].split(",")) if len(args) > 2 else (12.9, 29.1, 50.0)
    painel(runs, out, alvos)
