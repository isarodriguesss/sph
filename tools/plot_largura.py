"""Diagnostico da largura do braco: onde engorda, onde afina, e por que.

A largura local vem da transformada de distancia (EDT) dentro do corpo — o raio
do maior circulo inscrito, que e a definicao correta de espessura e nao depende
de escolher setor angular (licao #67-F: AR por setor confunde dedo com cone).

O painel de velocidade POR CLASSE e o que fecha o diagnostico: a biomassa viva
carrega o motor e acelera para fora, o limbo e passivo e so e arrastado. Onde a
razao entre as duas dispara, o limbo estagna e o material que vem atras empilha.

    python tools/plot_largura.py runs/P5_csmin010_REPROVADO [saida.png]
"""
import glob
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pysph.solver.utils import load
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree

DX = 0.0538
run = sys.argv[1] if len(sys.argv) > 1 else "runs/P5_csmin010_REPROVADO"
saida = sys.argv[2] if len(sys.argv) > 2 else "plots/largura.png"

b = None
for f in sorted(glob.glob(f"{run}/main_output/*.hdf5")):
    d = load(f)
    if b is None or d["solver_data"]["t"] > b[0]:
        b = (d["solver_data"]["t"], d["arrays"]["fluid"])
t, pa = b
rb = pa.rho_b_grown
fil = pa.is_filler > 0.5
r = np.hypot(pa.x, pa.y)
vel = np.hypot(pa.u, pa.v)
col = (rb >= 0.1) | fil
lim = (rb > 1e-6) & (rb < 0.1) & (~fil)
viva = (rb >= 0.1) & (~fil)
R99 = float(np.percentile(r[col], 99))

res = DX / 2
g = np.arange(-1.08 * R99, 1.08 * R99, res)
X, Y = np.meshgrid(g, g)
rg = np.hypot(X, Y)
th = np.arctan2(Y, X)
Pg = np.column_stack([X.ravel(), Y.ravel()])
mk = cKDTree(np.column_stack([pa.x[col | lim], pa.y[col | lim]])).query(Pg)[0].reshape(X.shape) < 0.75 * DX
edt = distance_transform_edt(mk) * res / DX * 2

fig = plt.figure(figsize=(19, 9.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.25, 1, 1], hspace=0.28, wspace=0.24)

ax = fig.add_subplot(gs[:, 0])
ax.set_facecolor("#0b0b0f")
W = np.where(mk, edt, np.nan)
im = ax.imshow(W, origin="lower", extent=[g[0], g[-1], g[0], g[-1]],
               cmap="turbo", vmin=0, vmax=8)
for f_, c_, l_ in ((0.65, "#00e5ff", "cintura  r/R99=0.65"), (0.83, "#ff2d95", "barriga  r/R99=0.83")):
    ax.add_patch(plt.Circle((0, 0), f_ * R99, fill=False, ec=c_, ls="--", lw=1.6, label=l_))
ax.legend(fontsize=9, loc="upper right", framealpha=.85)
ax.set_xticks([]); ax.set_yticks([])
ax.set_title("largura local do corpo (dx) — vermelho = gordo, azul = fino", fontsize=12)
plt.colorbar(im, ax=ax, fraction=0.045, label="largura (dx)")

# perfis por braco
NB = 360
occ = np.zeros(NB, bool)
for k in range(NB):
    lo = -np.pi + k * 2 * np.pi / NB
    occ[k] = (mk & (th >= lo) & (th < lo + 2 * np.pi / NB) & (rg > 0.85 * R99)).sum() > 2
lab = np.zeros(NB, int)
st = int(np.argmax(~occ)); c = 0; prev = False
for j in range(NB):
    k = (st + j) % NB
    if occ[k]:
        if not prev:
            c += 1
        lab[k] = c
    prev = occ[k]

rr = np.arange(0.40, 1.00, 0.05)
ax = fig.add_subplot(gs[0, 1])
prof_all = []
for i in range(1, c + 1):
    sel = np.zeros_like(mk)
    for k in np.where(lab == i)[0]:
        lo = -np.pi + k * 2 * np.pi / NB
        sel |= (th >= lo) & (th < lo + 2 * np.pi / NB)
    p = []
    for x in rr:
        m = mk & sel & (rg >= (x - .025) * R99) & (rg < (x + .025) * R99)
        p.append(np.median(edt[m]) if m.sum() > 4 else np.nan)
    prof_all.append(p)
    ax.plot(rr, p, lw=.9, alpha=.45, color="#888")
med = np.nanmedian(np.array(prof_all), axis=0)
ax.plot(rr, med, lw=2.6, color="#d62828", label="mediana dos %d bracos" % c)
ax.axvline(0.65, color="#00b8d4", ls="--", lw=1.2)
ax.axvline(0.83, color="#ff2d95", ls="--", lw=1.2)
ax.set_xlabel("r / R99"); ax.set_ylabel("largura (dx)")
ax.legend(fontsize=9); ax.grid(alpha=.3)
ax.set_title("perfil de largura: cintura em 0.65, barriga em 0.83", fontsize=11)

ax = fig.add_subplot(gs[0, 2])
bins = np.arange(0.35, 1.02, 0.07)
for sel, cor, lb in ((viva, "#2ecc71", "biomassa viva"), (fil, "#3a86ff", "filler"),
                     (lim, "#ff9f1c", "limbo")):
    y = []
    for a, b2 in zip(bins[:-1], bins[1:]):
        m = sel & (r >= a * R99) & (r < b2 * R99)
        y.append(np.median(vel[m]) if m.sum() > 3 else np.nan)
    ax.semilogy(0.5 * (bins[:-1] + bins[1:]), y, "o-", color=cor, label=lb, lw=1.8, ms=4)
ax.axvline(0.65, color="#00b8d4", ls="--", lw=1.2)
ax.axvline(0.83, color="#ff2d95", ls="--", lw=1.2)
ax.set_xlabel("r / R99"); ax.set_ylabel("|v| mediano")
ax.legend(fontsize=9); ax.grid(alpha=.3, which="both")
ax.set_title("a viva ACELERA, o limbo ESTAGNA", fontsize=11)

ax = fig.add_subplot(gs[1, 1])
y1, y2 = [], []
for a, b2 in zip(bins[:-1], bins[1:]):
    m = (col | lim) & (r >= a * R99) & (r < b2 * R99)
    y1.append(100 * np.mean(viva[m]) if m.sum() else np.nan)
    y2.append(100 * np.mean(((rb >= 0.2) & (rb <= 0.6) & (~fil))[m]) if m.sum() else np.nan)
xm = 0.5 * (bins[:-1] + bins[1:])
ax.plot(xm, y1, "o-", color="#2ecc71", label="% biomassa viva")
ax.plot(xm, y2, "s-", color="#8338ec", label="% no gate flagelar [0.2,0.6]")
ax.axvline(0.83, color="#ff2d95", ls="--", lw=1.2)
ax.set_xlabel("r / R99"); ax.set_ylabel("% do corpo")
ax.legend(fontsize=9); ax.grid(alpha=.3)
ax.set_title("quem pode se mover some antes da barriga", fontsize=11)

ax = fig.add_subplot(gs[1, 2])
razao = []
for a, b2 in zip(bins[:-1], bins[1:]):
    mv = viva & (r >= a * R99) & (r < b2 * R99)
    ml = lim & (r >= a * R99) & (r < b2 * R99)
    razao.append(np.median(vel[mv]) / np.median(vel[ml])
                 if mv.sum() > 3 and ml.sum() > 3 else np.nan)
ax.plot(xm, razao, "o-", color="#d62828", lw=2)
ax.axhline(1, color="k", ls=":", lw=1)
ax.axvline(0.83, color="#ff2d95", ls="--", lw=1.2)
ax.set_xlabel("r / R99"); ax.set_ylabel("|v| viva / |v| limbo")
ax.grid(alpha=.3)
ax.set_title("descolamento: 2x no corpo, 77x na ponta", fontsize=11)

fig.suptitle("%s  t=%.0f — a barriga e um ENGARRAFAMENTO: o limbo estagna onde a viva some"
             % (run.rstrip("/").split("/")[-1], t), fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.955])
fig.savefig(saida, dpi=118, facecolor="w")
print("->", saida)
