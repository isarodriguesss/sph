"""Largura do braco (EDT) e ciclo de trabalho azimutal, sobre a MASCARA rasterizada.

A v1 desta ferramenta amostrava as particulas numa fatia fina (|r-R| < 0.6 dx) e
contava arcos ocupados. Com ~170 particulas cobrindo 720 bins a ocupacao crua era
0.17 e o detector achava 93 "bracos" — fragmentacao de amostragem, nao geometria.
A largura saia 2x menor que a real (licao #101, correcao de 2026-09-16).

Aqui a colonia e rasterizada em dx/2, fechada e preenchida; a largura local e
2x a transformada de distancia (regua da licao #75-A, independente de setor
angular) e o ciclo e a ocupacao angular DA MASCARA.

    python tools/ciclo_largura.py runs/A runs/B [--t 50]
"""

import argparse
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, "plots")
import fig_tese as FT  # noqa: E402
from pysph.solver.utils import load  # noqa: E402

FRACS = (0.3, 0.45, 0.6, 0.75, 0.9)
NBIN = 720


def mascara(fl):
    rb = fl.rho_b_grown
    fil = fl.is_filler > 0.5
    r = np.hypot(fl.x, fl.y)
    col = (rb >= 0.1) | fil
    Rmax = np.percentile(r[col], 99.5)
    pix = FT.DX / 2.0
    L = Rmax * 1.25
    n = int(2 * L / pix)
    ix = ((fl.x[col] + L) / pix).astype(int)
    iy = ((fl.y[col] + L) / pix).astype(int)
    ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
    m = np.zeros((n, n), bool)
    m[iy[ok], ix[ok]] = True
    m = ndi.binary_closing(m, ndi.generate_binary_structure(2, 2), iterations=2)
    m = ndi.binary_fill_holes(m)
    return m, Rmax, pix, L


def main():
    p = argparse.ArgumentParser()
    p.add_argument("runs", nargs="+")
    p.add_argument("--t", type=float, default=50.0)
    a = p.parse_args()

    print(f"{'run':22s} {'Rmax':>5s} {'ocup':>5s} | largura EDT (dx) | ciclo angular")
    cab = f"{'':22s} {'':5s} {'':5s} |"
    for f in FRACS:
        cab += f" {f'{f:g}R':>5s}"
    cab += " |"
    for f in FRACS:
        cab += f" {f'{f:g}R':>5s}"
    print(cab)

    for run in a.runs:
        t, fname = min(FT.frames(run), key=lambda q: abs(q[0] - a.t))
        fl = load(fname)["arrays"]["fluid"]
        m, Rmax, pix, L = mascara(fl)
        edt = ndi.distance_transform_edt(m) * pix
        n = m.shape[0]
        yy, xx = np.mgrid[0:n, 0:n]
        rr = np.hypot(xx * pix - L, yy * pix - L)
        tt = np.arctan2(yy * pix - L, xx * pix - L)
        larg, ciclo = [], []
        for f in FRACS:
            anel = np.abs(rr - f * Rmax) < 1.0 * FT.DX
            sel = m & anel
            larg.append(2 * np.median(edt[sel]) / FT.DX if sel.sum() > 20 else np.nan)
            b = ((tt[anel] + np.pi) / (2 * np.pi) * NBIN).astype(int) % NBIN
            occ = np.zeros(NBIN, bool)
            occ[b[m[anel]]] = True
            ciclo.append(occ.mean())
        ocup = m.sum() * pix**2 / (np.pi * Rmax**2)
        linha = f"{run.split('/')[-1][:22]:22s} {Rmax:5.2f} {ocup:5.2f} |"
        linha += "".join(f" {v:5.1f}" for v in larg) + " |"
        linha += "".join(f" {v:5.2f}" for v in ciclo)
        print(linha)


if __name__ == "__main__":
    main()
