"""Campo `rho_b_vis` — representacao continua de `rho_b`, SEM tocar na dinamica.

O `rho_b` do C4 tem 77.6% de zeros EXATOS dentro do disco da colonia: a gaussiana
inicial `exp(-(r/R)**4)` faz underflow, e o `BiomassGrowth` e multiplicativo, entao zero
e estado absorvente (licao #48). O anel roxo no painel nao corresponde a ausencia de
colonia — corresponde a ausencia de representacao.

Este modulo reconstroi, SO PARA VISUALIZACAO, o envelope da colonia: particula com
`rho_b < 0.1` recebe a densidade do material de braco vizinho,

    rho_b_vis = (sum_j V_j rho_b_j^2 W_ij) / (sum_j V_j rho_b_j W_ij)

Media ponderada por BIOMASSA (nao Shepard): num campo esparso a Shepard nivela para baixo
e devolve ~0.03, sub-quorum (licao #59). A ponderada devolve a densidade do DOADOR.

E POS-PROCESSAMENTO PURO. `rho_b_vis` e variavel local desta funcao: nao existe no
particle_array, nao e lida por nenhuma equacao, nao e integrada e nao vai ao HDF5. A
dinamica e o C4 bit a bit.

MODOS (parametro `modo`):

  "iterativo"  (padrao) — apos o 1o passe, particula CERCADA (material cheio em >=6 dos
                8 setores angulares a <2h) tambem e preenchida, repetindo ate convergir.
                Fecha os buracos internos (1810 -> 46 em t=50) MAS tambem preenche as
                baias, porque assim que a borda da baia enche o interior dela vira
                cercado e a propagacao cascateia. O painel resultante e o ENVELOPE
                EXTERNO da colonia, nao a morfologia dendritica.

  "buracos"    — fecha apenas componentes conexas de nao-preenchidas com <= HOLE_MAX
                particulas. Medido em t=50: 99.3% do nao-preenchido e UMA componente de
                21 812 particulas (as baias sao um dominio unico ligado ao agar externo),
                entao este modo mexe em ~152 particulas e PRESERVA a morfologia
                dendritica no painel.

A escolha entre os dois nao e tecnica, e de qual figura se quer: envelope ou dendrito.

    python tools/rho_b_vis.py runs/C4 [saida.png] [modo]
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

NB = 72
FLOOR = 1e-6
MAX_ITER = 12
HOLE_MAX = 60
SETORES_CERCO = 6  # de 8


def r_base(x, y, rb):
    """Raio onde os dendritos se destacam: p25 do perfil azimutal do raio da colonia."""
    th = np.arctan2(y, x)
    r = np.hypot(x, y)
    body = rb > 0.1
    e = np.linspace(-np.pi, np.pi, NB + 1)
    rad = np.full(NB, np.nan)
    for i in range(NB):
        m = body & (th >= e[i]) & (th < e[i + 1])
        if m.sum() > 2:
            rad[i] = np.percentile(r[m], 95)
    ok = ~np.isnan(rad)
    rad = np.interp(np.arange(NB), np.where(ok)[0], rad[ok], period=NB)
    return float(np.percentile(rad, 25))


def preencher(pa, h_mult=1.0, modo="iterativo"):
    """Devolve (rho_b_vis, R_base). NAO modifica pa."""
    x, y, rb = pa.x, pa.y, pa.rho_b_grown
    h = float(pa.h[0])
    V = pa.m / pa.rho
    vis = rb.copy()
    R = r_base(x, y, rb)

    def transferir(mask_doador, alvos, exigir_cerco):
        dj = np.where(mask_doador)[0]
        if len(dj) == 0 or len(alvos) == 0:
            return 0
        tree = cKDTree(np.column_stack([x[dj], y[dj]]))
        raio = (2.0 if exigir_cerco else h_mult) * h
        n = 0
        for k, j in enumerate(tree.query_ball_point(np.column_stack([x[alvos], y[alvos]]), raio)):
            if not j:
                continue
            j = dj[np.asarray(j, dtype=int)]
            dx, dy = x[j] - x[alvos[k]], y[j] - y[alvos[k]]
            if exigir_cerco:
                setores = set(((np.arctan2(dy, dx) + np.pi) / (2 * np.pi) * 8).astype(int) % 8)
                if len(setores) < SETORES_CERCO:
                    continue
            w = np.exp(-((np.hypot(dx, dy) / h) ** 2))
            den = np.sum(V[j] * vis[j] * w)
            if den > 0:
                vis[alvos[k]] = np.sum(V[j] * vis[j] * vis[j] * w) / den
                n += 1
        return n

    transferir(rb > 0.1, np.where(rb < 0.1)[0], exigir_cerco=False)

    if modo == "iterativo":
        for _ in range(MAX_ITER):
            if transferir(vis >= 0.1, np.where(vis < 0.1)[0], exigir_cerco=True) == 0:
                break
    elif modo == "buracos":
        falta = np.where(vis < 0.1)[0]
        if len(falta) > 1:
            pares = cKDTree(np.column_stack([x[falta], y[falta]])).query_pairs(
                1.5 * h, output_type="ndarray"
            )
            g = coo_matrix(
                (np.ones(len(pares)), (pares[:, 0], pares[:, 1])),
                shape=(len(falta), len(falta)),
            )
            _, lab = connected_components(g, directed=False)
            tam = np.bincount(lab)
            transferir(vis >= 0.1, falta[tam[lab] <= HOLE_MAX], exigir_cerco=False)
    else:
        raise ValueError(f"modo desconhecido: {modo}")
    return vis, R


def main(run="runs/C4", out="plots/rho_b_vis.png", modo="iterativo"):
    f = sorted(glob.glob(f"{run}/main_output/*.hdf5"))[-1]
    d = load(f)
    pa = d["arrays"]["fluid"]
    t = d["solver_data"]["t"]
    vis, R = preencher(pa, modo=modo)
    x, y, rb = pa.x, pa.y, pa.rho_b_grown
    r = np.hypot(x, y)
    ch = (rb < 0.1) & (vis >= 0.1)
    Rout = np.percentile(r[rb > 0.1], 99)

    print(f"{run}  t={t:.1f}  modo={modo}  (R_base={R:.2f}, so referencia)")
    print(f"  preenchidas (rho_b<0.1 -> vis>=0.1) : {int(ch.sum())}")
    print(f"  rho_b_vis nas preenchidas           : mediana {np.median(vis[ch]):.3f}")
    print(f"  vazamento p/ agar externo (r>1.15R) : {int((ch & (rb < 1e-12) & (r > 1.15 * Rout)).sum())}")
    print(f"  rho_b REAL (nao tocado)             : soma {np.sum(rb):.6f}")

    L = 1.25 * Rout
    s = (np.abs(x) < L) & (np.abs(y) < L)
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 6.4))
    for a, campo, ttl in ((ax[0], rb, "rho_b REAL (dinamico)"),
                          (ax[1], vis, f"rho_b_vis ({modo})")):
        sc = a.scatter(x[s], y[s], c=np.maximum(campo[s], FLOOR), s=5,
                       cmap="viridis", norm=LogNorm(FLOOR, 1.0), linewidths=0)
        a.set_aspect("equal")
        a.set_xlim(-L, L)
        a.set_ylim(-L, L)
        a.set_facecolor("#101010")
        a.set_title(f"{ttl} — t={t:.1f}", fontsize=10)
        plt.colorbar(sc, ax=a, fraction=0.046)
    plt.tight_layout()
    plt.savefig(out, dpi=118)
    print(f"\n{out}")


if __name__ == "__main__":
    main(*(sys.argv[1:4] or ["runs/C4", "plots/rho_b_vis.png", "iterativo"]))
