"""Trajetoria da expansao: instantaneos + series temporais, nao so o estado final.

A linha de cima e `rho_b` em LOG (sem log, material entre 1e-3 e 0.1 le como fundo —
licoes #49/#56). A de baixo tras o criterio **C5 (§2.2)** — `R_conn/R99` e `frac_conn`,
que medem se a colonia esta conexa ao proprio centro em CADA instante — junto da
ocupacao da BANDA DE JUNCAO r [0.6, 1.8) e da mediana do limbo, o material que ainda
nao cruzou o quorum.

C5 e criterio de TRAJETORIA: uma rodada pode terminar aceitavel e viola-lo o tempo
todo. No E5 o raio conexo estagna em 0.66-0.82 enquanto `R99` vai a 4.13.

    python tools/plot_trajetoria.py runs/E5_hillK025 [saida.png]
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from pysph.solver.utils import load
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

DX = 0.0538
JU = (0.6, 1.8)
LINK = 1.05  # escala de CONTATO, em dx (a escala e parte do criterio — licao #68)


def serie(run):
    fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
    out = []
    for f in fs:
        d = load(f)
        pa = d["arrays"]["fluid"]
        rb = pa.rho_b_grown
        fil = pa.is_filler > 0.5
        env = (pa.is_env > 0.5) if "is_env" in pa.properties else np.zeros(len(rb), bool)
        r = np.hypot(pa.x, pa.y)
        corpo = (rb >= 0.1) & (~env)
        if corpo.sum() < 20:
            continue
        mat = corpo | fil | env
        i = np.where(mat)[0]
        P = np.column_stack([pa.x[i], pa.y[i]])
        tr = cKDTree(P)
        pp = tr.query_pairs(LINK * DX, output_type="ndarray")
        nn = len(i)
        _, lab = connected_components(
            coo_matrix((np.ones(len(pp)), (pp[:, 0], pp[:, 1])), shape=(nn, nn)),
            directed=False)
        nuc = lab == lab[int(np.argmin(np.hypot(P[:, 0], P[:, 1])))]
        R99 = float(np.percentile(r[corpo], 99))
        b = (r >= JU[0]) & (r < JU[1])
        lim = b & (rb > 1e-12) & (rb < 0.1) & (~fil) & (~env)
        area = np.pi * (JU[1] ** 2 - JU[0] ** 2) / DX / DX
        out.append(dict(
            t=d["solver_data"]["t"], f=f, R99=R99,
            c5a=float(np.hypot(P[nuc, 0], P[nuc, 1]).max()) / R99,
            c5b=100.0 * float(np.count_nonzero(nuc)) / nn,
            ocup=float(np.count_nonzero(b & (corpo | fil | env)) / area),
            nq=int(np.count_nonzero(b & corpo & (~fil))),
            lim=float(np.median(rb[lim])) if lim.sum() else 0.0,
        ))
    return out


def painel(run, out="plots/trajetoria.png", nsnap=4):
    S = serie(run)
    ts = np.array([s["t"] for s in S])
    alvos = np.linspace(ts.min(), ts.max(), nsnap)
    idx = [int(np.argmin(np.abs(ts - a))) for a in alvos]
    Rmax = max(s["R99"] for s in S)
    fig = plt.figure(figsize=(5.0 * nsnap, 9.6))
    gs = fig.add_gridspec(2, nsnap, height_ratios=[1.5, 1.0], hspace=0.22)
    for c, k in enumerate(idx):
        pa = load(S[k]["f"])["arrays"]["fluid"]
        rb = pa.rho_b_grown
        r = np.hypot(pa.x, pa.y)
        L = 1.25 * Rmax
        s = r < L
        ax = fig.add_subplot(gs[0, c])
        ax.set_facecolor("0.10")
        vv = np.maximum(rb, 1e-6)
        o = np.argsort(vv[s])
        ax.scatter(pa.x[s][o], pa.y[s][o], c=vv[s][o], s=9, cmap="viridis",
                   norm=LogNorm(vmin=1e-3, vmax=1), linewidths=0)
        for rad, cor in ((JU[0], "#ff3b8d"), (JU[1], "#ff3b8d"), (S[k]["R99"], "w")):
            ax.add_patch(plt.Circle((0, 0), rad, fill=False, ec=cor, ls="--", lw=1.0,
                                    alpha=0.55))
        ax.set_title(f't = {S[k]["t"]:.0f} s\nR99={S[k]["R99"]:.2f}   '
                     f'ocup juncao={S[k]["ocup"]:.2f}', color="w", fontsize=13)
        ax.set_xlim(-L, L)
        ax.set_ylim(-L, L)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
    for c, (ch, lb, cor) in enumerate((("R99", "R99 (raio da colonia)", "#4ec9f5"),
                                       ("c5a", "C5a  R_conn / R99   (alvo 1.0)", "#ff5555"),
                                       ("c5b", "C5b  frac_conn %   (alvo 100)", "#ff9f43"),
                                       ("lim", "rho_b mediano do limbo", "#f5a623"))):
        ax = fig.add_subplot(gs[1, c])
        v = [s[ch] for s in S]
        ax.plot(ts, v, "-o", color=cor, ms=3.5, lw=1.8)
        if ch == "c5a":
            ax.axhline(1.0, color="k", ls="--", lw=1, alpha=0.5)
            ax.set_ylim(0, 1.08)
        if ch == "c5b":
            ax.set_ylim(0, 105)
        if ch == "lim":
            ax.axhline(0.1, color="w", ls="--", lw=1, alpha=0.6)
            ax.text(ts[0], 0.103, "quorum 0.1", color="w", fontsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.set_title(lb, fontsize=12)
        ax.set_xlabel("t (s)")
        ax.grid(alpha=0.25)
    plt.savefig(out, dpi=118, facecolor="0.09", bbox_inches="tight")
    print("->", out)


if __name__ == "__main__":
    a = sys.argv[1:]
    painel(a[0] if a else ".", a[1] if len(a) > 1 else "plots/trajetoria.png")
