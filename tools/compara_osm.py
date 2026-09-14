"""OSM1 (absorcao osmotica) contra OSM0 (sem deposicao) e E11 — criterios pre-registrados
em runs/OSM1_osmotico/CRITERIOS.md.

COLONIA = rho_b>=0.1 ou filler ou phi_osm>0 (o absorvido e fase passiva, [T2]); a coluna
+limbo soma 0<rho_b<0.1 (§2.2). Forma (amplitude, dedos) sobre a COLONIA, em R99 igual;
C5 em media e desvio sobre t in [35, 50]. "Agar limpo" = fracao de agar puro (rho_b<1e-6,
nao filler, nao absorvido) no anel 0.5-0.9 R99, onde ficam as baias.

    python tools/compara_osm.py runs/E11_t100 runs/OSM0_semdep runs/OSM1_osmotico
"""
import glob
import os
import sys

import h5py
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

DX = 0.0538
T_MAX = 50.5
RAIOS = (1.2, 1.6, 2.0, 2.4, 2.8, 3.2)


def c5(x, y, sel):
    i = np.where(sel)[0]
    P = np.c_[x[i], y[i]]
    pp = cKDTree(P).query_pairs(1.05 * DX, output_type="ndarray")
    n = len(i)
    _, lab = connected_components(
        coo_matrix((np.ones(len(pp)), (pp[:, 0], pp[:, 1])), shape=(n, n)), directed=False)
    nuc = lab == lab[int(np.argmin(np.hypot(P[:, 0], P[:, 1])))]
    return float(np.hypot(P[nuc, 0], P[nuc, 1]).max()), 100.0 * nuc.sum() / n


def forma(r, th, sel, nb=128):
    ib = np.clip(((th[sel] + np.pi) / (2 * np.pi) * nb).astype(int), 0, nb - 1)
    Rt = np.zeros(nb)
    np.maximum.at(Rt, ib, r[sel])
    amp = Rt.std() / max(Rt.mean(), 1e-12)
    ded = int(((Rt > np.roll(Rt, 1)) & (Rt > np.roll(Rt, -1))).sum())
    return amp, ded


def serie(run):
    out = []
    for f in sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5"))):
        with h5py.File(f, "r") as h:
            t = float(h["solver_data"].attrs["t"])
            if t > T_MAX:
                break
            a = h["particles"]["fluid"]["arrays"]
            x, y, rb, isf = (np.asarray(a[k], float) for k in ("x", "y", "rho_b_grown", "is_filler"))
            osm = np.asarray(a["phi_osm"], float) > 0 if "phi_osm" in a else np.zeros(len(x), bool)
        port = rb > 0.1
        if port.sum() < 20:
            continue
        r, th = np.hypot(x, y), np.arctan2(y, x)
        fil = isf > 0.5
        viva = (~fil) & port
        col = (rb >= 0.1) | fil | osm
        coll = col | (rb > 1e-12)
        R99 = float(np.percentile(r[port], 99))
        amp, ded = forma(r, th, col)
        ampv, dedv = forma(r, th, viva)
        ra, fb = c5(x, y, col)
        ral, fbl = c5(x, y, coll)
        an = (r >= 0.5 * R99) & (r < 0.9 * R99)
        limpo = an & (rb < 1e-6) & (~fil) & (~osm)
        out.append(dict(t=t, R99=R99, amp=amp, ded=ded, ampv=ampv, dedv=dedv,
                        c5a=ra / R99, c5b=fb, c5bl=fbl, limpo=limpo.sum() / max(an.sum(), 1),
                        viva=int(viva.sum()), fil=int(fil.sum()), osm=int(osm.sum())))
    return out


def main():
    runs = sys.argv[1:]
    S = {run: serie(run) for run in runs}
    print("\n== C5 e composicao, media +- desvio em t in [35, 50] ==")
    print(f"{'run':>22} {'n':>2} {'C5b col':>13} {'C5b +limbo':>13} {'C5a':>11} "
          f"{'agar limpo':>11} {'vivas':>6} {'filler':>6} {'absorv':>6} {'R99 fim':>7}")
    for run in runs:
        w = [s for s in S[run] if 35.0 <= s["t"] <= T_MAX]
        if not w:
            print(f"{os.path.basename(run):>22}  (sem frames em [35,50])")
            continue
        m = lambda k: (np.mean([s[k] for s in w]), np.std([s[k] for s in w]))  # noqa: E731
        print(f"{os.path.basename(run):>22} {len(w):2d} {m('c5b')[0]:6.1f}+-{m('c5b')[1]:4.1f} "
              f"{m('c5bl')[0]:6.1f}+-{m('c5bl')[1]:4.1f} {m('c5a')[0]:5.2f}+-{m('c5a')[1]:4.2f} "
              f"{100 * m('limpo')[0]:10.1f}% {w[-1]['viva']:6d} {w[-1]['fil']:6d} {w[-1]['osm']:6d} "
              f"{w[-1]['R99']:7.2f}")

    print("\n== forma em R99 IGUAL (colonia | so vivas) ==")
    print(f"{'run':>22} {'R99':>4} {'t':>5} {'amp col':>7} {'dedos':>5} | {'amp viva':>8} "
          f"{'dedos':>5} | {'C5b':>5} {'limpo':>6}")
    for run in runs:
        s = S[run]
        R = np.array([q["R99"] for q in s])
        for Ra in RAIOS:
            if R.max() < Ra:
                break
            i = int(np.argmax(R >= Ra))
            j = s[max(0, i - 1):i + 2]
            md = lambda k: np.median([q[k] for q in j])  # noqa: E731
            print(f"{os.path.basename(run):>22} {Ra:4.1f} {s[i]['t']:5.1f} {md('amp'):7.3f} "
                  f"{md('ded'):5.0f} | {md('ampv'):8.3f} {md('dedv'):5.0f} | {md('c5b'):5.1f} "
                  f"{100 * md('limpo'):5.1f}%")


if __name__ == "__main__":
    main()
