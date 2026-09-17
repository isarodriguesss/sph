"""Compara RUNS com as REFERENCIAS do §2.2 medindo os dois com o MESMO algoritmo.

    python tools/lit_rank.py runs/P2R10_filler_conduz runs/P2R11_filler_conduz_lento ...

As referencias entram como IMAGEM: assets/reference_result.png painel (b) de Trinschek
(mapa de h = colonia, e mapa de Gamma, que e invertido pela propria barra de cores da figura,
entao sai na unidade do artigo, Gamma_max = 0.5) e assets/reference.jpg paineis A e C (PA14,
limiar no vale do histograma). Os runs entram pelo campo da colonia e pelo campo de `c_s` do
plots/fig_tese.py no instante da figura (t = 50). Tudo adimensional, dividido por Rmax ou pelo
teto do campo (Gamma_max la, `cs_max` aqui), entao os numeros sao comparaveis.

Forma:  nucleo = maior raio em que a colonia cobre >= 98% das direcoes (o nucleo solido);
        baia = Rmin/Rmax (o fundo da baia mais profunda); n dedos (picos de R(theta) com
        relevo >= 6% de Rmax e separacao >= 10 graus); area/disco de Rmax;
        AR = (Rmax - nucleo)/largura do braco (EDT), so ordem de grandeza nas fotos.
Campo:  mediana no corpo e nas baias sobre o teto; valor em 0.9 Rmax e em Rmax;
        r(50%) = raio onde a mediana azimutal cruza metade do teto, sobre Rmax;
        L = decaimento fora da colonia medido pela DISTANCIA a borda (imune a forma).

Componentes: nas fotos/figuras conta so o maior (tira ruido); nos runs, todos com area >= 1%
do maior — e o que o tools/rank_runs e o plots/fig_tese medem (bracos soltos entram).
"""
import sys

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.signal import find_peaks

NTH = 720
TETO = 0.5  # Gamma_max de Trinschek (b) e cs_max do nosso modelo


def rtheta(mask, centro=None, nth=NTH, so_maior=True):
    """so_maior=True (fotos/figuras: descarta ruido) | False (runs: todos os pedacos >= 1%,
    como o tools/rank_runs e o plots/fig_tese, que medem o campo inteiro)."""
    lab, n = ndi.label(mask)
    area = ndi.sum(mask, lab, np.arange(1, n + 1))
    keep = (area == area.max()) if so_maior else (area >= 0.01 * area.max())
    m = ndi.binary_fill_holes((lab > 0) & keep[lab - 1])
    cy, cx = ndi.center_of_mass(m) if centro is None else centro[::-1]
    r = np.arange(0.0, 0.995 * min(m.shape) / 2, 0.5)
    th = np.linspace(-np.pi, np.pi, nth, endpoint=False)
    X = cx + r[None, :] * np.cos(th[:, None])
    Y = cy + r[None, :] * np.sin(th[:, None])
    dentro = ndi.map_coordinates(m.astype(float), [Y, X], order=0) > 0.5
    Rt = np.where(dentro.any(1), r[dentro.shape[1] - 1 - np.argmax(dentro[:, ::-1], 1)], 0.0)
    ocup = dentro.mean(0)                    # fracao de direcoes com colonia no raio r
    i = np.argmax(ocup < 0.98)               # nucleo solido: cobre todas as direcoes
    Rc = r[i] if i > 0 else (r[-1] if ocup.all() else 0.0)
    return Rt, (cx, cy), m, Rc


def dedos(Rt):
    s = ndi.uniform_filter1d(np.r_[Rt, Rt, Rt], 5)[len(Rt):2 * len(Rt)]
    pk, _ = find_peaks(np.r_[s, s], prominence=0.06 * s.max(), distance=NTH / 36)
    return int(np.sum((pk >= len(s) // 2) & (pk < 3 * len(s) // 2)))


def forma(mask, so_maior=True):
    Rt, centro, m, Rc = rtheta(mask, so_maior=so_maior)
    R1, R0 = Rt.max(), Rt.min()
    edt = ndi.distance_transform_edt(m)
    gy, gx = np.mgrid[:m.shape[0], :m.shape[1]]
    rr = np.hypot(gx - centro[0], gy - centro[1])
    braco = m & (rr > R0) & (rr < 0.95 * R1)
    larg = 2 * np.median(edt[braco]) if braco.any() else np.nan
    return dict(R1=R1, R0=R0, rmin_rmax=R0 / R1, nucleo=Rc / R1, n=dedos(Rt),
                ocup=m.sum() / (np.pi * R1 ** 2), AR=(R1 - R0) / larg if larg else np.nan,
                centro=centro, mask=m)


def campo_stats(F, mask_col, centro, R1, teto=TETO):
    baia = ndi.binary_fill_holes(ndi.binary_closing(mask_col, np.ones((int(R1 / 2) | 1,) * 2)))
    baia &= ~ndi.binary_dilation(mask_col, iterations=2)
    r = np.arange(0.0, 0.995 * min(F.shape) / 2, 0.5)
    th = np.linspace(-np.pi, np.pi, NTH, endpoint=False)
    X = centro[0] + r[None, :] * np.cos(th[:, None])
    Y = centro[1] + r[None, :] * np.sin(th[:, None])
    perf = np.median(ndi.map_coordinates(F, [Y, X], order=1, mode="nearest"), 0)
    meio = (perf < 0.50 * teto) & (r > 0.5 * R1)
    r50 = r[np.argmax(meio)] if meio.any() else np.nan
    # decaimento medido pela DISTANCIA a borda da colonia (nao pelo raio): imune a forma
    dout = ndi.distance_transform_edt(~mask_col)
    bins = np.arange(0.0, 0.5 * R1, max(1.0, R1 / 40))
    med = np.array([np.median(F[(dout > a) & (dout <= b)]) if ((dout > a) & (dout <= b)).any()
                    else np.nan for a, b in zip(bins[:-1], bins[1:])])
    cen = 0.5 * (bins[:-1] + bins[1:])
    ok = np.isfinite(med) & (med > 0.02 * teto) & (med < 0.9 * teto)
    L = -1.0 / np.polyfit(cen[ok], np.log(med[ok]), 1)[0] if ok.sum() > 3 else np.nan
    return dict(corpo=np.median(F[mask_col]) / teto,
                baia=np.median(F[baia]) / teto if baia.any() else np.nan,
                p90=np.interp(0.9 * R1, r, perf) / teto,
                ponta=np.interp(R1, r, perf) / teto, halo=r50 / R1, L_R=L / R1)


def linha(nome, f, c=None):
    s = (f"{nome:24s} nucleo {f['nucleo']:.2f} R | baia {f['rmin_rmax']:.2f} R | dedos "
         f"{f['n']:3d} | area/disco {f['ocup']:.2f} | AR {f['AR']:.1f}")
    if c:
        s += (f" || campo: corpo {c['corpo']:.2f} baia {c['baia']:.2f} 0.9R {c['p90']:.2f}"
              f" ponta {c['ponta']:.2f} | r(50%) {c['halo']:.2f} R | L {c['L_R']:.2f} R")
    print(s)


# ---------------- referencias ----------------
im = np.asarray(Image.open("assets/reference_result.png").convert("RGB"), float)
h_box = im[124:305, 1719:1906]
oliva, azul = np.array([154, 152, 85]), np.array([183, 221, 232])
col_T = (np.linalg.norm(h_box - oliva, axis=2) < np.linalg.norm(h_box - azul, axis=2))
fT = forma(col_T)

g_box = im[411:591, 1719:1906]
bar = im[412:589, 1922:1932].mean(1)
vals = np.linspace(TETO, 0.0, len(bar))
j = np.argmin(((g_box.reshape(-1, 1, 3) - bar[None]) ** 2).sum(2), axis=1)
G = vals[j].reshape(g_box.shape[:2])
zy, zx = np.array(g_box.shape[:2]) / np.array(h_box.shape[:2])
maskG = ndi.zoom(fT["mask"].astype(float), (zy, zx), order=0) > 0.5
cT = campo_stats(G, maskG, (fT["centro"][0] * zx, fT["centro"][1] * zy), fT["R1"] * zx)
linha("Trinschek (b) fingering", fT, cT)

ref = np.asarray(Image.open("assets/reference.jpg").convert("RGB"), float)
def vale(v):
    h, be = np.histogram(v, 64, (0, 256))
    h = ndi.gaussian_filter1d(h.astype(float), 1.5)
    c = (be[:-1] + be[1:]) / 2
    pico = int(np.argmax(h[:32]))
    mins = [i for i in np.where((h[1:-1] < h[:-2]) & (h[1:-1] < h[2:]))[0] + 1 if i > pico]
    return c[mins[0]] if mins else c[pico] + 40


for nome, (x0, x1) in (("A", (3, 158)), ("C", (322, 476))):
    p = ref[4:152, x0:x1].max(2)
    m = ndi.binary_opening(p > vale(p), np.ones((2, 2)))
    linha(f"PA14 painel {nome}", forma(m))

# ---------------- runs ----------------
if len(sys.argv) > 1:
    sys.path.insert(0, "plots")
    import fig_tese as FT
    from scipy.spatial import cKDTree
    for run in sys.argv[1:]:
        t, f = min(FT.frames(run), key=lambda p: abs(p[0] - 50))
        d = FT.carrega(f)
        rc = np.hypot(d["x"], d["y"])[FT.colonia(d, False)]
        L = float(np.ceil(1.12 * np.percentile(rc, 99.9) * 2) / 2)
        n = 500
        g, img, borda, _ = FT.campo(d, L, n, False, agar_iso=2.0)
        fR = forma(borda, so_maior=False)
        vol = d["m"] / d["rho"]
        gx, gy = np.meshgrid(g, g)
        pts = np.c_[gx.ravel(), gy.ravel()]
        viva = d["is_filler"] < 0.5
        cs, den = FT.shepard(d["x"][viva], d["y"][viva], d["cs"][viva], vol[viva], pts, 1.8 * FT.DX)
        vazio = den <= 1e-9
        if vazio.any():
            _, jj = cKDTree(np.c_[d["x"][viva], d["y"][viva]]).query(pts[vazio])
            cs[vazio] = d["cs"][viva][jj]
        cR = campo_stats(cs.reshape(n, n), fR["mask"], fR["centro"], fR["R1"])
        # mesmo campo, agora COM o filler (runs que conduzem): e o painel (d) --cs-filler
        cs2, den2 = FT.shepard(d["x"], d["y"], d["cs"], vol, pts, 1.8 * FT.DX)
        cR2 = campo_stats(cs2.reshape(n, n), fR["mask"], fR["centro"], fR["R1"])
        nome = run.split("/")[-1][:22]
        linha(nome, fR, cR)
        print(f"{'':24s} (com filler)                                             "
              f"   || campo: corpo {cR2['corpo']:.2f} baia {cR2['baia']:.2f} 0.9R {cR2['p90']:.2f}"
              f" ponta {cR2['ponta']:.2f} | r(50%) {cR2['halo']:.2f} R | L {cR2['L_R']:.2f} R")
