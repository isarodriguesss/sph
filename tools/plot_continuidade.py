"""Continuidade: quem tem caminho de particulas ate o centro, e em que escala.

Responde a pergunta "os bracos estao unidos ao nucleo?" sem depender da cor do
painel de `rho_b` — que confunde colonia com limbo (licao #53: material em
(0, 0.1) e mecanicamente invisivel e quimicamente ativo).

Duas colunas deliberadamente separadas:
  COLONIA   = `rho_b >= 0.1` ou `is_filler`  (definicao da §2.2)
  + LIMBO   = idem, mais `0 < rho_b < 0.1`   (o que a cor SUGERE ser colonia)

A diferenca entre as duas e o tamanho da ilusao. A escala de ligacao e varrida
(1.05 a 2.7 dx) porque conectividade medida numa escala so e circular — licao #68.

    python tools/plot_continuidade.py runs/P2_fillerdonor . [saida.png]
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pysph.solver.utils import load
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

DX = 0.0538
ESCALAS = [1.05, 1.2, 1.4, 1.7, 2.0, 2.7]


def ultimo(run):
    b = None
    for f in sorted(glob.glob(f"{run}/main_output/*.hdf5")):
        d = load(f)
        t = d["solver_data"]["t"]
        if b is None or t > b[0]:
            b = (t, d["arrays"]["fluid"])
    return b


def conexo(P, link):
    tr = cKDTree(P)
    pp = tr.query_pairs(link * DX, output_type="ndarray")
    n = len(P)
    if len(pp) == 0:
        return np.zeros(n, bool)
    _, lab = connected_components(
        coo_matrix((np.ones(len(pp)), (pp[:, 0], pp[:, 1])), shape=(n, n)),
        directed=False)
    return lab == lab[int(np.argmin(np.hypot(P[:, 0], P[:, 1])))]


runs = [a for a in sys.argv[1:] if not a.endswith(".png")]
saida = next((a for a in sys.argv[1:] if a.endswith(".png")), "plots/continuidade.png")

fig = plt.figure(figsize=(6.0 * 2 + 6.4, 5.6 * len(runs)))
gs = fig.add_gridspec(len(runs), 3, width_ratios=[1, 1, 1.05], wspace=0.12, hspace=0.15)
serie = {}

for i, run in enumerate(runs):
    t, pa = ultimo(run)
    rb = pa.rho_b_grown
    fil = pa.is_filler > 0.5
    r = np.hypot(pa.x, pa.y)
    col = (rb >= 0.1) | fil
    limbo = (rb > 1e-6) & (rb < 0.1) & (~fil)
    R99 = float(np.percentile(r[col], 99))
    nome = run.rstrip("/").split("/")[-1] or "atual"

    for j, (rotulo, mask) in enumerate([("COLONIA  (rho_b>=0.1 ou filler)", col),
                                        ("COLONIA + LIMBO", col | limbo)]):
        ax = fig.add_subplot(gs[i, j])
        ax.set_facecolor("#0c0c10")
        idx = np.where(mask)[0]
        P = np.column_stack([pa.x[idx], pa.y[idx]])
        c = conexo(P, 1.05)
        ax.scatter(P[~c, 0], P[~c, 1], s=2.2, c="#d0342c", lw=0)
        ax.scatter(P[c, 0], P[c, 1], s=2.2, c="#2ecc71", lw=0)
        ax.add_patch(plt.Circle((0, 0), R99, fill=False, ec="w", ls="--", lw=1, alpha=.5))
        ax.set_xlim(-1.15 * R99, 1.15 * R99)
        ax.set_ylim(-1.15 * R99, 1.15 * R99)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("%s\nligado ao centro: %d de %d = %.1f%%"
                     % (rotulo, c.sum(), len(idx), 100 * c.mean()),
                     fontsize=10.5, color="w")
        if j == 0:
            ax.set_ylabel("%s\nt=%.0f" % (nome, t), fontsize=12)

    ax = fig.add_subplot(gs[i, 2])
    for rotulo, mask, cor in [("colonia", col, "#1f77b4"), ("+ limbo", col | limbo, "#ff7f0e")]:
        idx = np.where(mask)[0]
        P = np.column_stack([pa.x[idx], pa.y[idx]])
        fr, re = [], []
        for L in ESCALAS:
            c = conexo(P, L)
            fr.append(100 * c.mean())
            re.append(np.percentile(np.hypot(P[c, 0], P[c, 1]), 99) / R99 if c.sum() > 5 else 0)
        serie.setdefault(nome, {})[rotulo] = (fr, re)
        ax.plot(ESCALAS, fr, "o-", color=cor, label="%s  frac ligada" % rotulo)
        ax.plot(ESCALAS, [100 * v for v in re], "s--", color=cor, alpha=.55,
                label="%s  alcance R_conn/R99" % rotulo)
    ax.axvline(1.05, color="k", ls=":", lw=1)
    ax.text(1.07, 94, "contato", fontsize=8.5, rotation=90, va="top")
    ax.set_xlabel("escala de ligacao do grafo (dx)")
    ax.set_ylabel("%")
    ax.set_ylim(0, 100)
    ax.grid(alpha=.3)
    ax.legend(fontsize=8, loc="lower right")
    ax.set_title("varredura de escala (licao #68)", fontsize=10.5)

fig.suptitle("Continuidade — verde = tem caminho de particulas ate o centro a 1.05 dx (contato); "
             "vermelho = nao tem", fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(saida, dpi=115, facecolor="w")
print("->", saida)
for nome, d in serie.items():
    for rot, (fr, re) in d.items():
        print("  %-18s %-9s frac %s" % (nome, rot, " ".join("%5.1f" % v for v in fr)))
