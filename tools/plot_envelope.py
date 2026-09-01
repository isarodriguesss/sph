"""Expansao como CORPO CONTINUO — envelope da biomassa, nao particulas soltas.

Uso:  python tools/plot_envelope.py runs/C4 [runs/N2_t50 ...] [--times 15 30 48]
                                            [--out FIG.png]

Renderizar `rho_b` por particulas mostra um esqueleto esparso e sugere uma colonia
cheia de buracos; renderizar por campo Shepard mostra um corpo solido e esconde a
fragmentacao. Nenhum dos dois responde a pergunta que interessa: **a colonia e um
corpo unico que se expande, ou pedacos desconexos?**

Aqui a colonia e definida por TOPOLOGIA, nao por um disco:

  1. particulas com `rho_b > 0.1` (limiar de quorum do Hill — abaixo disso a
     particula nao produz surfactante, nao passa no gate flagelar e nao tem
     pressao na EOS, logo nao e colonia em nenhum sentido funcional);
  2. grafo de vizinhanca com aresta se a distancia < 2h (suporte do kernel);
  3. componentes conexas. A MAIOR e o corpo; as demais sao fragmentos.
  4. envelope = regiao a menos de 1.5h de uma particula do corpo.

O disco R99 usado nas metricas antigas conta as BAIAS como colonia e infla a area
~3x numa morfologia dendritica (C4: 66.7 contra 20.6). O envelope segue os braços.

Painel esquerdo por instante: corpo (escuro) e fragmentos desconexos (claro).
Painel da direita: contornos do envelope superpostos no tempo — a expansao vista
como sequencia de frentes, que e a leitura de "expansao continua".
"""

import sys
import os
import re
import csv
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
import h5py
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

DX = 14.0 / 260
H = 1.8 * DX
NGRID = 340
QUORUM = 0.1


def load(run):
    rows = list(csv.DictReader(open(os.path.join(run, "log.csv"))))
    it2t = {int(r["iteration"]): float(r["t"]) for r in rows}
    its = np.array(sorted(it2t))
    ts = np.array([it2t[i] for i in its])
    out = []
    for fn in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        it = int(re.search(r"main_(\d+)", fn).group(1))
        out.append((float(np.interp(it, its, ts)), fn))
    return out


def bodies(fn):
    """Corpo (maior componente conexa da biomassa) e fragmentos."""
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        x, y = np.array(a["x"]), np.array(a["y"])
        rb = np.array(a["rho_b_grown"])
    bio = rb > QUORUM
    if bio.sum() < 10:
        return np.zeros((0, 2)), np.zeros((0, 2)), 0
    P = np.column_stack([x, y])
    pairs = cKDTree(P[bio]).query_pairs(2 * H, output_type="ndarray")
    n = int(bio.sum())
    if len(pairs) == 0:
        return P[bio], np.zeros((0, 2)), n
    adj = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    ncomp, lab = connected_components(adj, directed=False)
    big = np.bincount(lab).argmax()
    idx = np.where(bio)[0]
    return P[idx[lab == big]], P[idx[lab != big]], ncomp


def envelope(pts, lim):
    """Mascara booleana da regiao a <1.5h de alguma particula do corpo."""
    g = np.linspace(-lim, lim, NGRID)
    GX, GY = np.meshgrid(g, g)
    if len(pts) == 0:
        return GX, GY, np.zeros(GX.shape, bool)
    d, _ = cKDTree(pts).query(np.column_stack([GX.ravel(), GY.ravel()]))
    return GX, GY, (d < 1.5 * H).reshape(GX.shape)


def main(argv):
    times, out, runs = [15.0, 30.0, 48.0], "runs/envelope.png", []
    i = 0
    while i < len(argv):
        if argv[i] == "--times":
            times, i = [], i + 1
            while i < len(argv) and not argv[i].startswith("--"):
                times.append(float(argv[i]))
                i += 1
        elif argv[i] == "--out":
            out, i = argv[i + 1], i + 2
        else:
            runs.append(argv[i])
            i += 1
    if not runs:
        raise SystemExit(__doc__)

    ncol = len(times) + 1
    fig, axes = plt.subplots(len(runs), ncol, figsize=(4.3 * ncol, 4.5 * len(runs)))
    axes = np.atleast_2d(axes)

    for ri, run in enumerate(runs):
        series = load(run)
        lim = 5.0
        cmapt = cm.viridis(np.linspace(0.15, 0.95, len(series)))
        for ci, tq in enumerate(times):
            t, fn = min(series, key=lambda s: abs(s[0] - tq))
            body, frag, ncomp = bodies(fn)
            GX, GY, env = envelope(body, lim)
            ax = axes[ri, ci]
            ax.set_facecolor("#f7f7f2")
            ax.contourf(
                GX, GY, env.astype(float), levels=[0.5, 1.5], colors=["#1b5e20"]
            )
            if len(frag):
                ax.scatter(frag[:, 0], frag[:, 1], s=7, c="#ef6c00", lw=0, zorder=3)
            area = env.sum() * (2 * lim / NGRID) ** 2
            ax.set_title(
                f"t = {t:.0f} s   ·   {ncomp} componentes\n"
                f"corpo = {100 * len(body) / max(len(body) + len(frag), 1):.0f}% "
                f"da biomassa   ·   área {area:.1f}",
                fontsize=9,
            )
            ax.set_aspect("equal")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_xticks([])
            ax.set_yticks([])
            if ci == 0:
                ax.set_ylabel(
                    os.path.basename(run.rstrip("/")), fontsize=13, fontweight="bold"
                )

        ax = axes[ri, -1]
        ax.set_facecolor("#f7f7f2")
        for k, (t, fn) in enumerate(series):
            body, _, _ = bodies(fn)
            if len(body) < 10:
                continue
            GX, GY, env = envelope(body, lim)
            ax.contour(
                GX,
                GY,
                env.astype(float),
                levels=[0.5],
                colors=[cmapt[k]],
                linewidths=1.3,
            )
        ax.set_title(
            "frentes do envelope no tempo\n(escuro → claro = t crescente)", fontsize=9
        )
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        "Colonia por TOPOLOGIA — corpo = maior componente conexa de ρ_b > 0.1\n"
        "verde = envelope do corpo · laranja = fragmentos desconexos",
        fontsize=14,
        fontweight="bold",
    )
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    plt.savefig(out, dpi=115, bbox_inches="tight")
    print(f"figura: {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
