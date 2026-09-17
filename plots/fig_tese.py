"""Figura de tese de UM run, no formato do painel (b) de Trinschek et al. 2018 [T1].

    python plots/fig_tese.py runs/E11_t100 [--t 50] [--tempos 10 20 30 40 50]
                             [--out plots/fig_E11] [--limbo] [--fecha 0] [--n 700]

Mesma disposicao e mesmas cores de Trinschek: (a) colonia no instante final, oliva sobre
azul-claro; (b) contornos da borda nos `--tempos`, em preto; (c) raios maximo (vermelho) e
minimo (azul-petroleo) contra t, com os instantes de (b) marcados; (d) surfactante `c_s`
azul-branco-vermelho com o contorno da colonia — o halo alem da biomassa. Diferenca
deliberada: `c_s` em escala LOG fixa (Trinschek usa linear), sem a qual o halo some.

COLONIA = `rho_b >= 0.1` ou filler (definicao do §2.2); `--limbo` soma `0 < rho_b < 0.1`.
O painel (a) e a regiao SEM agar preenchida (colonia, limbo e o vazio entre elas); cada
particula de agar ocupa a sua celula da rede e nenhum buraco e preenchido por padrao — todo
buraco contem agar (licao #92; `--fecha` > 0 volta ao preenchimento da licao #72-K). Agar
isolado (sem outro agar a < `--iso` dx, padrao 2) nao e desenhado: e um buraco de uma particula,
abaixo do suporte do kernel. `c_s` (Shepard, Price 2007) exclui o filler, cujo valor e congelado
(licao #39); onde nao ha particula nao-filler no suporte (zona so de filler mais larga que 2h),
recebe o `c_s` da nao-filler mais proxima (decisao da usuaria, 2026-09-15) — extrapolacao visual. Escalas fixas: figuras de runs diferentes sao comparaveis.
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
CMAP_COL = LinearSegmentedColormap.from_list("colonia", ["#b7dde8", "#9a9855"])
CMAP_COL_CONTORNO = LinearSegmentedColormap.from_list(
    "colonia_contorno", ["#b7dde8", "#f3f1dc", "#9a9855"])
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


def campo(d, L, n, limbo, fecha=0.0, r_agar=0.8, ponte_agar=0.0, area_min=4.0,
          migalha=0.0, agar_iso=0.0):
    """Colonia na grade: a regiao SEM agar, preenchida.

    Agar = `rho_b` = 0 e nao filler. Um ponto e colonia se a particula de agar mais proxima
    esta a mais de `r_agar`*dx (a rede do agar deixa no maximo 0.71 dx ate um no): entram
    colonia, limbo e o vazio entre elas. A mascara do agar recebe um fechamento de
    `ponte_agar`*dx, que costura linhas de agar deslocado — sem isso bracos vizinhos se unem
    pelas frestas; depois, grao de agar solto com menos de `migalha` dx^2 (1-2 particulas
    deslocadas dentro do rastro) sai. Manchas de colonia menores que `area_min` dx^2 saem.
    Nenhum buraco maior e preenchido por padrao (licao #92); `fecha` > 0 volta a preencher
    buracos cercados de ate `fecha`*dx; `limbo` fica pela assinatura. `agar_iso` > 0 ignora a
    particula de agar sem nenhuma outra de agar a menos de `agar_iso`*dx (buraco de uma
    particula, abaixo do suporte do kernel 2h = 3.6 dx); as figuras usam 2, as metricas 0.

    Devolve (g, img, borda, area_preenchida): `img` para exibir (borda suavizada em ~0.6 dx);
    `borda` com TODOS os buracos fechados — so o contorno externo, que (b) e (c) medem.
    """
    g = np.linspace(-L, L, n)
    GX, GY = np.meshgrid(g, g)
    pts = np.c_[GX.ravel(), GY.ravel()]
    sel = (np.abs(d["x"]) < 1.2 * L) & (np.abs(d["y"]) < 1.2 * L)
    agar = (d["is_filler"] < 0.5) & (d["rho_b_grown"] <= 1e-12)
    xa, ya = d["x"][sel & agar], d["y"][sel & agar]
    if agar_iso > 0:
        d2, _ = cKDTree(np.c_[xa, ya]).query(np.c_[xa, ya], k=2)
        manter = d2[:, 1] <= agar_iso * DX
        xa, ya = xa[manter], ya[manter]
    da, _ = cKDTree(np.c_[xa, ya]).query(pts)
    cel = 2 * L / (n - 1)
    mask_agar = (da <= r_agar * DX).reshape(n, n)
    if ponte_agar > 0:
        kp = max(1, int(round(ponte_agar * DX / cel)))
        yy, xx = np.mgrid[-kp:kp + 1, -kp:kp + 1]
        # borda estendida com agar: o fechamento nao pode corroer a moldura da janela
        pad = np.pad(mask_agar, kp + 1, constant_values=True)
        mask_agar = ndi.binary_closing(pad, structure=(xx * xx + yy * yy) <= kp * kp)[
            kp + 1:-kp - 1, kp + 1:-kp - 1]
    lab, n_lab = ndi.label(mask_agar)
    if n_lab:
        area = ndi.sum(mask_agar, lab, index=np.arange(1, n_lab + 1)) * cel * cel
        mask_agar = (lab > 0) & (area >= migalha * DX * DX)[lab - 1]
    corpo = ~mask_agar
    lab, n_lab = ndi.label(corpo)
    if n_lab:
        area = ndi.sum(corpo, lab, index=np.arange(1, n_lab + 1)) * cel * cel
        corpo = (lab > 0) & (area >= area_min * DX * DX)[lab - 1]
    buraco = np.zeros_like(corpo)
    if fecha > 0:
        k = max(1, int(round(fecha * DX / cel)))
        yy, xx = np.mgrid[-k:k + 1, -k:k + 1]
        fechado = ndi.binary_closing(corpo, structure=(xx * xx + yy * yy) <= k * k)
        lab, n_lab = ndi.label(ndi.binary_fill_holes(corpo) & ~corpo)
        if n_lab:
            vedado = ndi.sum(~fechado, lab, index=np.arange(1, n_lab + 1)) == 0
            buraco = (lab > 0) & vedado[lab - 1]
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
    ap.add_argument("--fecha", type=float, default=0.0)
    ap.add_argument("--n", type=int, default=700)
    ap.add_argument("--L", type=float, default=None, help="meia-largura da janela")
    ap.add_argument("--contorno", action="store_true", help="borda clara no painel (a)")
    ap.add_argument("--borda-preta", action="store_true", help="linha preta na borda do painel (a)")
    ap.add_argument("--iso", type=float, default=2.0, help="ignora agar isolado (dx); 0 desliga")
    ap.add_argument("--cs-filler", action="store_true", help="inclui o filler no painel (d) (runs com FILLER_CS_CONDUZ)")
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

    g, img, _, area_buraco = campo(d_fim, L, a.n, a.limbo, a.fecha, agar_iso=a.iso)
    ext = [-L, L, -L, L]
    ax[0, 0].imshow(img, origin="lower", extent=ext,
                    cmap=CMAP_COL_CONTORNO if a.contorno else CMAP_COL, vmin=0, vmax=1,
                    interpolation="bilinear")
    if a.borda_preta:
        ax[0, 0].contour(g, g, img, levels=[0.5], colors=COR_BORDA, linewidths=0.6)
    eixos(ax[0, 0], L)
    ax[0, 0].set_title(rf"(a) colônia, $t={t_fim:.0f}$", loc="left")

    n_lo = max(300, a.n // 2)
    serie = {}
    for tt, f in fr:
        gg, _, borda, _ = campo(carrega(f), L, n_lo, a.limbo, a.fecha, agar_iso=a.iso)
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
    viva = (d_fim["is_filler"] < 0.5) | a.cs_filler
    cs, den = shepard(d_fim["x"][viva], d_fim["y"][viva], d_fim["cs"][viva], vol[viva],
                      pts, 1.8 * DX)
    vazio = den <= 1e-9
    if vazio.any():
        _, j = cKDTree(np.c_[d_fim["x"][viva], d_fim["y"][viva]]).query(pts[vazio])
        cs[vazio] = d_fim["cs"][viva][j]
    cs = cs.reshape(a.n, a.n)
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
