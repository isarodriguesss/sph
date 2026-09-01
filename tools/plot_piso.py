"""Compara rodadas em `rho_b` (linear + log) e por CLASSE de particula.

Com `--fill[=N]` aplica PREENCHIMENTO TOPOLOGICO na renderizacao: rasteriza o corpo,
fecha vaos de ate N*dx (padrao 3.5) e inunda a partir de FORA — o que a inundacao nao
alcanca e buraco, e e desenhado como colonia. Baia e ligada ao exterior por construcao,
entao NUNCA e preenchida; nao ha raio de busca para vazar.

E pos-processamento puro: nao toca no solver. Medido no solver (runs/E10_topo) o mesmo
mecanismo custava `a_pressure` mediana 2.64 -> 3.24 e `mean_v` a menor da serie; aqui
custa zero. E a licao #38 — vacuo VISUAL se resolve na renderizacao, nunca adicionando
particula ao solver.

Tres linhas por coluna: `rho_b` em escala LINEAR (como o viewer mostra), `rho_b` em
escala LOG (sem a qual material entre 1e-3 e 0.1 le como fundo — licoes #49/#56/#72-G)
e a CLASSE de cada particula, porque agar e vazio renderizam iguais em qualquer painel
de campo (licao #66).

As classes separam o piso (`is_env`) da matriz (`is_filler`) e da biomassa viva, o que
importa porque toda metrica normalizada por "quem e colonia" muda de significado quando
uma rota cria portadoras (licao #62).

    python tools/plot_piso.py runs/E5_hillK025 . [saida.png] [t]
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from pysph.solver.utils import load
from scipy import ndimage as ndi

FUNDO = "0.10"
CORES = (("agar", "#2b2b38"), ("limbo", "#4a4a5e"), ("preenchido", "#ff3b8d"),
         ("matriz/filler", "#f5a623"), ("biomassa VIVA", "#2ecc71"))


def preenche(pa, corpo, cand, fecha, cell=0.5):
    """Buracos do corpo, por inundacao a partir do exterior. So renderizacao."""
    dx = cell * float(pa.h[0]) / 1.8
    rc = 1.2 * float(np.percentile(np.hypot(pa.x[corpo], pa.y[corpo]), 99))
    nc = int(2.0 * rc / dx) + 1
    ix = np.clip(((pa.x + rc) / dx).astype(int), 0, nc - 1)
    iy = np.clip(((pa.y + rc) / dx).astype(int), 0, nc - 1)
    A = np.zeros((nc, nc), dtype=bool)
    A[ix[corpo], iy[corpo]] = True
    k = int(np.ceil(fecha / cell))
    yy, xx = np.mgrid[-k:k + 1, -k:k + 1]
    st = (xx * xx + yy * yy) <= k * k  # disco: quadrado deixa borda em degrau
    buraco = ndi.binary_fill_holes(ndi.binary_closing(A, structure=st)) & (~A)
    return cand & buraco[ix, iy]


def carrega(run, alvo):
    fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
    ts = [load(f)["solver_data"]["t"] for f in fs]
    d = load(fs[int(np.argmin([abs(t - alvo) for t in ts]))])
    return d["arrays"]["fluid"], d["solver_data"]["t"]


def painel(runs, out="plots/piso.png", alvo=50.0, zoom=1.5, fecha=None):
    fig, axes = plt.subplots(3, len(runs), figsize=(6.7 * len(runs), 19.6), squeeze=False)
    for c, run in enumerate(runs):
        pa, t = carrega(run, alvo)
        nome = run.rstrip("/").split("/")[-1] or "atual"
        x, y, rb = pa.x, pa.y, pa.rho_b_grown
        r = np.hypot(x, y)
        env = (pa.is_env > 0.5) if "is_env" in pa.properties else np.zeros(len(x), bool)
        fil = pa.is_filler > 0.5
        if fecha:
            env = env | preenche(pa, (rb >= 0.1) | fil, (rb < 0.1) & (~fil) & (~env),
                                 fecha)
            rb = np.where(env & (rb < 0.1), 0.1, rb)
        corpo = (rb >= 0.1) & (~env)
        R99 = float(np.percentile(r[corpo], 99))
        L = zoom * R99
        s = r < L

        ax = axes[0, c]
        ax.set_facecolor(FUNDO)
        o = np.argsort(rb[s])
        ax.scatter(x[s][o], y[s][o], c=rb[s][o], s=13, cmap="viridis", vmin=0, vmax=1,
                   linewidths=0)
        ax.set_title(f"{nome}   t={t:.1f}\nrho_b LINEAR", color="w", fontsize=13)

        ax = axes[1, c]
        ax.set_facecolor(FUNDO)
        vv = np.maximum(rb, 1e-6)
        o = np.argsort(vv[s])
        ax.scatter(x[s][o], y[s][o], c=vv[s][o], s=13, cmap="viridis",
                   norm=LogNorm(vmin=1e-3, vmax=1), linewidths=0)
        ax.set_title("rho_b LOG (1e-3 .. 1)", color="w", fontsize=13)

        ax = axes[2, c]
        ax.set_facecolor(FUNDO)
        grupos = (s & (rb <= 1e-12) & (~env), s & (rb < 0.1) & (rb > 1e-12) & (~env),
                  s & env, s & corpo & fil, s & corpo & (~fil))
        hs = []
        for m, (lb, cor) in zip(grupos, CORES):
            if m.any():
                ax.scatter(x[m], y[m], c=cor, s=13, linewidths=0)
                hs.append(Line2D([], [], marker="o", ls="", color=cor,
                                 label=f"{lb}  ({int(m.sum())})"))
        ax.legend(handles=hs, loc="upper right", fontsize=9, framealpha=0.9)
        ax.set_title("classes", color="w", fontsize=13)

        for row in range(3):
            a = axes[row, c]
            a.add_patch(plt.Circle((0, 0), R99, fill=False, ec="w", ls="--", lw=1.1,
                                   alpha=0.55))
            a.set_xlim(-L, L)
            a.set_ylim(-L, L)
            a.set_aspect("equal")
            a.set_xticks([])
            a.set_yticks([])
    fig.patch.set_facecolor("0.06")
    plt.tight_layout()
    plt.savefig(out, dpi=124, facecolor="0.06")
    print("->", out)


if __name__ == "__main__":
    a = sys.argv[1:]
    fecha = None
    for w in [w for w in a if w.startswith("--fill")]:
        a.remove(w)
        fecha = float(w.split("=")[1]) if "=" in w else 3.5
    t = 50.0
    if a and a[-1].replace(".", "").isdigit():
        t = float(a.pop())
    out = "plots/piso.png"
    if a and a[-1].endswith(".png"):
        out = a.pop()
    painel(a or ["."], out, t, fecha=fecha)
