"""Expansao vista pelos tres campos escalares — CAMPO CONTINUO, nao pontos.

Uso:  python tools/plot_fields.py runs/C4_t100 [--times 15 40 55 99]
                                               [--out FIG.png] [--n 420]

Os campos sao reconstruidos por interpolacao SPH Shepard-normalizada
(Price 2007, `splash`, PASA 24:159):

    A(x) = Σ_j V_j W(|x−x_j|, h) A_j  /  Σ_j V_j W(|x−x_j|, h)

A normalizacao pelo denominador e o que torna o campo continuo mesmo onde a
amostragem e irregular — e a mesma razao pela qual `sigma_a` mede consistencia.
Renderizar por pontos mostra as PARTICULAS; renderizar o campo mostra o FLUIDO,
que e o objeto fisico.

Linhas = instantes; colunas = os tres escalares que governam a dinamica:

  rho_b   biomassa local (0-1). Contorno em 0.8 = limiar do hard pin (K.17):
          dentro dele o nucleo e imovel.
  c_s     surfactante. A forca de Marangoni e -beta*grad(c_s) — o que move a
          colonia e o GRADIENTE. O halo alem da biomassa e a assinatura do
          painel (b) de Trinschek [T1]; contorno branco = borda da biomassa.
  c_n     nutriente. Reservatorio FINITO (inicia em 1.0, sem fonte). Contorno
          em 0.4 = piso do gate de crescimento; dentro dele a colonia PARA.

Particulas inseridas (`is_filler`) sao quimicamente transparentes: seu `c_s` fica
congelado e nao representa o campo, por isso sao EXCLUIDAS da reconstrucao de c_s.
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
import h5py
from scipy.spatial import cKDTree

KNN = 40

# (chave, rotulo, cmap, vmin, vmax, contornos, cor do contorno, nota)
FIELDS = [
    ("rho_b_grown", "biomassa  ρ_b", "YlGn", 0.0, 0.9,
     [0.1, 0.8], "#1b5e20",
     "contorno 0.8 = núcleo pinado (imóvel)"),
    ("cs", "surfactante  c_s", "YlOrRd", 0.0, 0.35,
     [], "#4a148c",
     "motor = −β·∇c_s  (o gradiente, não o valor)"),
    ("c_n", "nutriente  c_n", "Blues", 0.0, 1.0,
     [0.4], "#b71c1c",
     "contorno 0.4 = piso do gate; dentro dele o crescimento para"),
]


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


def spline(q):
    """Forma radial do cubic spline (a constante cancela na razao Shepard)."""
    w = np.zeros_like(q)
    m1 = q < 1.0
    m2 = (q >= 1.0) & (q < 2.0)
    w[m1] = 1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3
    w[m2] = 0.25 * (2.0 - q[m2]) ** 3
    return w


def shepard(px, py, vals, vol, h, grid_pts, shape):
    """Campo Shepard-normalizado numa grade regular."""
    tree = cKDTree(np.column_stack([px, py]))
    k = min(KNN, len(px))
    dist, idx = tree.query(grid_pts, k=k)
    if k == 1:
        dist, idx = dist[:, None], idx[:, None]
    w = spline(dist / h) * vol[idx]
    den = np.sum(w, axis=1)
    num = np.sum(w * vals[idx], axis=1)
    return (num / (den + 1e-12)).reshape(shape)


def main(argv):
    times = [15.0, 40.0, 55.0, 99.0]
    out = "runs/campos.png"
    n = 420
    runs = []
    i = 0
    while i < len(argv):
        if argv[i] == "--times":
            times = []
            i += 1
            while i < len(argv) and not argv[i].startswith("--"):
                times.append(float(argv[i]))
                i += 1
        elif argv[i] == "--out":
            out = argv[i + 1]
            i += 2
        elif argv[i] == "--n":
            n = int(argv[i + 1])
            i += 2
        else:
            runs.append(argv[i])
            i += 1
    run = runs[0] if runs else "runs/C4_t100"
    series = load(run)

    nr, nc = len(times), len(FIELDS)
    fig, axes = plt.subplots(nr, nc, figsize=(5.4 * nc, 5.2 * nr))
    if nr == 1:
        axes = axes[None, :]

    lim = None
    grid_pts = gx = gy = None
    for ri, tq in enumerate(times):
        t, fn = min(series, key=lambda s: abs(s[0] - tq))
        f = h5py.File(fn, "r")
        a = f["particles"]["fluid"]["arrays"]
        d = {k: np.array(a[k]) for k in a if np.array(a[k]).size}
        x, y, h = d["x"], d["y"], d["h"]
        vol = d["m"] / np.maximum(d["rho"], 1e-9)
        if lim is None:
            lim = float(np.max(np.abs(np.concatenate([x, y]))))
            g = np.linspace(-lim, lim, n)
            gx, gy = np.meshgrid(g, g)
            grid_pts = np.column_stack([gx.ravel(), gy.ravel()])
        isf = d.get("is_filler")
        real = (isf < 0.5) if isf is not None and isf.size == x.size else np.ones_like(x, bool)

        rb_field = shepard(x, y, d["rho_b_grown"], vol, float(h[0]), grid_pts, gx.shape)

        for ci, (key, label, cmap, vmin, vmax, levels, lc, note) in enumerate(FIELDS):
            ax = axes[ri, ci]
            if key not in d:
                ax.axis("off")
                continue
            if key == "cs":  # filler tem cs congelado — nao e campo
                fld = shepard(x[real], y[real], d[key][real], vol[real],
                              float(h[0]), grid_pts, gx.shape)
            elif key == "rho_b_grown":
                fld = rb_field
            else:
                fld = shepard(x, y, d[key], vol, float(h[0]), grid_pts, gx.shape)

            im = ax.imshow(fld, origin="lower", extent=[-lim, lim, -lim, lim],
                           cmap=cmap, vmin=vmin, vmax=vmax, interpolation="bilinear")
            if levels:
                ax.contour(gx, gy, fld, levels=levels, colors=lc,
                           linewidths=1.1, alpha=0.85)
            if key == "cs":  # borda da biomassa, p/ enxergar o halo alem dela
                ax.contour(gx, gy, rb_field, levels=[0.1], colors="#1b5e20",
                           linewidths=1.2, alpha=0.9)
            cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
            cb.set_label(label, fontsize=9)
            ax.set_aspect("equal")
            if ri == 0:
                ax.set_title(f"{label}\n{note}", fontsize=10, fontweight="bold")
            if ci == 0:
                ax.set_ylabel(f"t = {t:.0f} s", fontsize=13, fontweight="bold")
            ax.set_xticks([])
            ax.set_yticks([])

    fig.suptitle(
        f"Expansão pelos campos escalares (interpolação SPH Shepard) — "
        f"{os.path.basename(run)}\n"
        "ρ_b estrutura a colônia · c_s a move (via ∇c_s) · c_n a limita",
        fontsize=15, fontweight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    plt.savefig(out, dpi=100, bbox_inches="tight")
    print(f"figura salva em {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
