"""Monitora se a mitose esta nascendo onde o vazio esta, e a que custo de pressao.

As filhas sao apensadas ao fim do array, entao `indice >= n_anterior` as isola sem
precisar de rastreamento. Filha = nao-filler (mitose); filler = wake/insert.

O vazio e definido por RASTERIZACAO (celula dentro de R99 sem particula a < 1.0 dx) e as
filhas pelo evento de divisao — parametros DIFERENTES, para nao repetir a metrica circular
da licao #68, em que o alcance do mecanismo era a propria regua.
"""

import argparse
import glob
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

DX = 0.0538


def frames(run):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
    out = []
    for f in fs:
        with h5py.File(f, "r") as h:
            out.append((float(h["solver_data"].attrs["t"]), f))
    return out


def le(f):
    with h5py.File(f, "r") as h:
        a = h["particles"]["fluid"]["arrays"]
        return {k: np.asarray(a[k]) for k in ("x", "y", "rho_b_grown", "is_filler", "u", "v")}


def vazio(x, y, corpo, R99, passo=0.5, fecha=3.5):
    """Buraco TOPOLOGICO: fecha vaos de ate `fecha`*dx e inunda a partir de FORA.

    O que a inundacao nao alcanca e buraco. A baia e ligada ao exterior por construcao,
    entao NUNCA e marcada — sem isso a metrica conta baia como defeito e premia o disco
    (licao #67-E), que e o modo de falha do P8.
    """
    from scipy import ndimage

    if corpo.sum() == 0:
        return np.empty((0, 2))
    lim = R99 * 1.25
    g = np.arange(-lim, lim + passo * DX, passo * DX)
    GX, GY = np.meshgrid(g, g)
    d, _ = cKDTree(np.c_[x[corpo], y[corpo]]).query(np.c_[GX.ravel(), GY.ravel()])
    ocup = (d <= 1.0 * DX).reshape(GX.shape)

    raio = max(1, int(round(fecha / (2 * passo))))
    yy, xx = np.ogrid[-raio : raio + 1, -raio : raio + 1]
    disco = xx * xx + yy * yy <= raio * raio
    fechado = ndimage.binary_closing(ocup, structure=disco)

    livre = ~fechado
    rot, _ = ndimage.label(livre)
    exterior = set(np.unique(np.r_[rot[0], rot[-1], rot[:, 0], rot[:, -1]])) - {0}
    buraco = livre & ~np.isin(rot, list(exterior))
    return np.c_[GX[buraco], GY[buraco]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", default=".")
    ap.add_argument("--ate", type=float, default=1e9)
    args = ap.parse_args()

    fs = frames(args.run)
    prev_n = None
    print(f"{'t':>6} {'filhas':>7} {'deposito':>9} | {'d(filha,vazio) p50':>18} {'<2dx':>6} |"
          f" {'vazio dx2':>10} {'%p!=0':>7} {'|v| anel':>9}")
    for t, f in fs:
        if t > args.ate:
            break
        d = le(f)
        x, y, rb, isf = d["x"], d["y"], d["rho_b_grown"], d["is_filler"]
        n = len(x)
        r = np.hypot(x, y)
        corpo = (rb >= 0.1) | (isf > 0.5)
        R99 = float(np.percentile(r[rb > 0.1], 99)) if (rb > 0.1).sum() > 10 else 1.0
        vz = vazio(x, y, corpo, R99)
        area = len(vz) * (0.5 * DX) ** 2 / DX**2

        anel = (r > 0.15 * R99) & (r < 0.30 * R99) & corpo
        vmed = np.median(np.hypot(d["u"][anel], d["v"][anel])) if anel.sum() else np.nan

        nf, dep, dm, frac2 = 0, 0, np.nan, np.nan
        if prev_n is not None and n > prev_n:
            novas = np.arange(prev_n, n)
            fil = novas[isf[novas] < 0.5]
            nf, dep = len(fil), len(novas) - len(fil)
            if len(fil) and len(vz):
                dd, _ = cKDTree(vz).query(np.c_[x[fil], y[fil]])
                dm, frac2 = np.median(dd) / DX, 100 * np.mean(dd < 2 * DX)
        prev_n = n
        print(f"{t:6.1f} {nf:7d} {dep:9d} | {dm:18.2f} {frac2:5.0f}% |"
              f" {area:10.0f} {'':7} {vmed:9.2e}")


if __name__ == "__main__":
    main()
