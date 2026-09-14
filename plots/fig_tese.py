"""Figura de tese de UM run, no formato do painel (b) de Trinschek et al. 2018 [T1].

    python plots/fig_tese.py runs/E11_t100 [--t 50] [--tempos 10 20 30 40 50]
                             [--out plots/fig_E11] [--limbo] [--fecha 3.5] [--n 700]

Mesma disposicao e mesmas cores de Trinschek: (a) colonia no instante final, oliva sobre
azul-claro; (b) contornos da borda nos `--tempos`, em preto; (c) raios maximo (vermelho) e
minimo (azul-petroleo) contra t, com os instantes de (b) marcados; (d) surfactante `c_s`
azul-branco-vermelho com o contorno da colonia — o halo alem da biomassa. Diferenca
deliberada: `c_s` em escala LOG fixa (Trinschek usa linear), sem a qual o halo some.

COLONIA = `rho_b >= 0.1` ou filler (definicao do §2.2); `--limbo` soma `0 < rho_b < 0.1`.
O campo e a fracao de colonia reconstruida por Shepard (Price 2007, PASA 24:159) com
h = dx, que preserva bracos de 1-2 dx. Buracos FECHADOS de ate `--fecha`*dx sao
preenchidos por inundacao a partir de fora — so na renderizacao (licao #72-K); baia e
ligada ao exterior e nunca e preenchida. `c_s` exclui o filler, cujo valor e congelado
(licao #39). Escalas fixas: figuras de runs diferentes sao comparaveis.
"""

import argparse
import glob
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import MultipleLocator
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

DX = 0.0538
KNN = 32
CS_MIN, CS_MAX = 1e-3, 0.5
# paleta do painel (b) de Trinschek et al. 2018 (assets/reference_result.png)
CMAP_COL = LinearSegmentedColormap.from_list("colonia", ["#b7dde8", "#f3f1dc", "#9a9855"])
CMAP_CS = LinearSegmentedColormap.from_list(
    "cs", ["#156f8c", "#6cb3c8", "#f5f5f5", "#d23a3a", "#8e0000"])
COR_BORDA = "#2b2b2b"
COR_MAX = "#c0272d"
COR_MIN = "#1b7a98"

plt.rcParams.update({
    "font.family": "serif", "mathtext.fontset": "cm", "font.size": 9,
    "axes.linewidth": 0.8, "xtick.direction": "out", "ytick.direction": "out",
})


def frames(run):
    out = []
    for f in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        with h5py.File(f, "r") as h:
            out.append((float(h["solver_data"].attrs["t"]), f))
    return out


def carrega(f):
    with h5py.File(f, "r") as h:
        a = h["particles"]["fluid"]["arrays"]
        d = {k: np.asarray(a[k], float) for k in
             ("x", "y", "m", "rho", "rho_b_grown", "is_filler", "cs")}
    return d


def colonia(d, limbo):
    rb = d["rho_b_grown"]
    c = (rb >= 0.1) | (d["is_filler"] > 0.5)
    if limbo:
        c |= rb > 1e-12
    return c


def spline(q):
    w = np.zeros_like(q)
    m1 = q < 1.0
    m2 = (q >= 1.0) & (q < 2.0)
    w[m1] = 1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3
    w[m2] = 0.25 * (2.0 - q[m2]) ** 3
    return w


def shepard(x, y, val, vol, pts, h):
    dist, idx = cKDTree(np.c_[x, y]).query(pts, k=KNN)
    w = spline(dist / h) * vol[idx]
    den = w.sum(1)
    return (w * val[idx]).sum(1) / np.maximum(den, 1e-12), den


def campo(d, L, n, limbo, fecha):
    """Colonia na grade.

    Devolve (g, img, borda, area_preenchida): `img` e a colonia para exibir (buracos de
    ate `fecha`*dx fechados, borda suavizada em ~0.6 dx); `borda` tem TODOS os buracos
    fechados — so o contorno externo, que e o que (b) e (c) medem.
    """
    g = np.linspace(-L, L, n)
    GX, GY = np.meshgrid(g, g)
    pts = np.c_[GX.ravel(), GY.ravel()]
    sel = (np.abs(d["x"]) < 1.2 * L) & (np.abs(d["y"]) < 1.2 * L)
    vol = d["m"] / d["rho"]
    chi, den = shepard(d["x"][sel], d["y"][sel], colonia(d, limbo)[sel].astype(float),
                       vol[sel], pts, DX)
    corpo = (np.where(den > 1e-9, chi, 0.0) > 0.5).reshape(n, n)
    cel = 2 * L / (n - 1)
    k = max(1, int(round(fecha * DX / cel)))
    yy, xx = np.mgrid[-k:k + 1, -k:k + 1]
    buraco = ndi.binary_fill_holes(
        ndi.binary_closing(corpo, structure=(xx * xx + yy * yy) <= k * k)) & ~corpo
    img = ndi.gaussian_filter((corpo | buraco).astype(float), 0.6 * DX / cel)
    borda = ndi.binary_fill_holes(img > 0.5)
    return g, img, borda, buraco.sum() * cel * cel


def raios(g, borda, nth=360):
    """Raios maximo e minimo de R(theta), o raio mais externo em cada direcao."""
    cel = g[1] - g[0]
    r = np.arange(0.0, g[-1], 0.5 * cel)
    th = np.linspace(-np.pi, np.pi, nth, endpoint=False)
    X = r[None, :] * np.cos(th[:, None])
    Y = r[None, :] * np.sin(th[:, None])
    dentro = ndi.map_coordinates(borda.astype(float), [(Y - g[0]) / cel, (X - g[0]) / cel],
                                 order=0) > 0.5
    Rt = np.where(dentro.any(1), r[dentro.shape[1] - 1 - np.argmax(dentro[:, ::-1], 1)], 0.0)
    return Rt.max(), Rt.min()


def eixos(ax, L, rot=True):
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_aspect("equal")
    ax.xaxis.set_major_locator(MultipleLocator(2))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.set_xlabel(r"$x/L$")
    if rot:
        ax.set_ylabel(r"$y/L$")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--t", type=float, default=50.0)
    ap.add_argument("--tempos", type=float, nargs="+", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--limbo", action="store_true")
    ap.add_argument("--fecha", type=float, default=3.5)
    ap.add_argument("--n", type=int, default=700)
    ap.add_argument("--L", type=float, default=None, help="meia-largura da janela")
    a = ap.parse_args()

    fr = [(t, f) for t, f in frames(a.run) if t <= a.t + 1e-6]
    t_fim, f_fim = fr[-1]
    tempos = a.tempos or list(np.linspace(t_fim / 5, t_fim, 5))
    esc = [min(fr, key=lambda p: abs(p[0] - tt)) for tt in tempos]
    esc = sorted(set(esc))
    nome = os.path.basename(a.run.rstrip("/"))
    out = a.out or f"plots/fig_{nome}"

    d_fim = carrega(f_fim)
    c = colonia(d_fim, a.limbo)
    r_col = np.hypot(d_fim["x"][c], d_fim["y"][c])
    L = a.L or float(np.ceil(1.12 * np.percentile(r_col, 99.9) * 2) / 2)

    fig, ax = plt.subplots(2, 2, figsize=(6.3, 6.0), constrained_layout=True)

    g, img, _, area_buraco = campo(d_fim, L, a.n, a.limbo, a.fecha)
    ext = [-L, L, -L, L]
    ax[0, 0].imshow(img, origin="lower", extent=ext, cmap=CMAP_COL, vmin=0, vmax=1,
                    interpolation="bilinear")
    eixos(ax[0, 0], L)
    ax[0, 0].set_title(rf"(a) colônia, $t={t_fim:.0f}$", loc="left")

    n_lo = max(300, a.n // 2)
    serie = {}
    for tt, f in fr:
        gg, _, borda, _ = campo(carrega(f), L, n_lo, a.limbo, a.fecha)
        serie[tt] = (gg, borda)

    for tt, _ in esc:
        gg, borda = serie[tt]
        ax[0, 1].contour(gg, gg, borda.astype(float), levels=[0.5], colors=COR_BORDA,
                         linewidths=0.8)
    eixos(ax[0, 1], L)
    ax[0, 1].set_title("(b) borda da colônia", loc="left")

    ts, rmax, rmin = [], [], []
    for tt, _ in fr:
        R1, R0 = raios(*serie[tt])
        ts.append(tt)
        rmax.append(R1)
        rmin.append(R0)
    ax[1, 0].plot(ts, rmax, color=COR_MAX, lw=1.4, label="máx")
    ax[1, 0].plot(ts, rmin, color=COR_MIN, lw=1.4, label="mín")
    for tt, _ in esc:
        i = ts.index(tt)
        ax[1, 0].plot(tt, rmax[i], "o", color="k", ms=3)
        ax[1, 0].plot(tt, rmin[i], "s", color="k", ms=3)
    ax[1, 0].set_xlabel(r"$t$")
    ax[1, 0].set_ylabel(r"raio $/L$")
    ax[1, 0].set_xlim(0, t_fim)
    ax[1, 0].set_ylim(0, None)
    ax[1, 0].legend(fontsize=7, frameon=False, loc="upper left")
    ax[1, 0].set_title("(c) raios extremos", loc="left")
    ax[1, 0].set_box_aspect(1)

    gx, gy = np.meshgrid(g, g)
    pts = np.c_[gx.ravel(), gy.ravel()]
    vol = d_fim["m"] / d_fim["rho"]
    viva = d_fim["is_filler"] < 0.5
    cs, den = shepard(d_fim["x"][viva], d_fim["y"][viva], d_fim["cs"][viva], vol[viva],
                      pts, 1.8 * DX)
    cs = np.where(den > 1e-9, cs, CS_MIN).reshape(a.n, a.n)
    im = ax[1, 1].imshow(np.clip(cs, CS_MIN, CS_MAX), origin="lower", extent=ext,
                         cmap=CMAP_CS, norm=LogNorm(CS_MIN, CS_MAX),
                         interpolation="bilinear")
    ax[1, 1].contour(g, g, img, levels=[0.5], colors="k", linewidths=0.6)
    eixos(ax[1, 1], L)
    ax[1, 1].set_title(rf"(d) surfactante $c_s$, $t={t_fim:.0f}$", loc="left")
    cb = fig.colorbar(im, ax=ax[1, 1], shrink=0.85, pad=0.02)
    cb.set_label(r"$c_s$", rotation=0, labelpad=8)

    for ext_ in ("pdf", "png"):
        fig.savefig(f"{out}.{ext_}", dpi=300)
    print(f"{out}.pdf / .png  (L={L}, buracos preenchidos na renderizacao: "
          f"{area_buraco / DX ** 2:.0f} dx^2, frames: {[round(t, 1) for t, _ in esc]})")


if __name__ == "__main__":
    main()
