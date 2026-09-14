"""Viscosidade artificial de Monaghan como o solver faz (c_ij = cs medio, o SURFACTANTE) contra
como deveria (c_ij = c0), por populacao, no estado real do P2."""
import sys

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "tools")
import diag_forcas as DF  # noqa: E402

ALPHA, C0, GAMMA = 0.12, 0.35, 60.0
DX = 14.0 / 260
rng = np.random.default_rng(0)


def a_visc_art(d, sel):
    x, y = d["x"], d["y"]
    H = float(np.median(d["h"]))
    viz = cKDTree(np.c_[x, y]).query_ball_point(np.c_[x[sel], y[sel]], 2.0 * H)
    ii = np.concatenate([np.full(len(v), k) for k, v in enumerate(viz)])
    jj = np.concatenate(viz).astype(int)
    di = sel[ii]
    keep = jj != di
    ii, jj, di = ii[keep], jj[keep], di[keep]
    xij, yij = x[di] - x[jj], y[di] - y[jj]
    r = np.hypot(xij, yij)
    hij = 0.5 * (d["h"][di] + d["h"][jj])
    g = DF.dw(r, hij) / np.maximum(r, 1e-12)
    DXg, DYg = g * xij, g * yij
    mj, rbar = d["m"][jj], 0.5 * (d["rho"][di] + d["rho"][jj])
    vdotr = (d["u"][di] - d["u"][jj]) * xij + (d["v"][di] - d["v"][jj]) * yij
    mu = hij * vdotr / (r * r + 0.01 * hij * hij)
    base = np.where(vdotr < 0.0, -ALPHA * mu / rbar, 0.0)
    n = len(sel)
    out = {}
    for nome, cij in (("atual", 0.5 * (d["cs"][di] + d["cs"][jj])), ("c0", np.full(len(di), C0))):
        pi = base * cij
        ax = np.bincount(ii, weights=-mj * pi * DXg, minlength=n)
        ay = np.bincount(ii, weights=-mj * pi * DYg, minlength=n)
        out[nome] = np.hypot(ax, ay)
    return out


for run, tt in (("runs/P2_t100", 50), ("runs/P2_t100", 75), ("runs/P2_t100", 85)):
    t, d = DF.carrega(run, tt)
    rb, isf, isw = d["rho_b_grown"], d["is_filler"], d["is_wake"]
    r = np.hypot(d["x"], d["y"])
    R99 = np.percentile(r[(rb > 0.1) | (isf > 0.5)], 99)
    pops = {
        "viva": (isf < 0.5) & (rb >= 0.1),
        "limbo": (isf < 0.5) & (rb > 1e-12) & (rb < 0.1),
        "filler wake": (isf > 0.5) & (isw > 0.5),
        "agar ate 1.2 R99": (isf < 0.5) & (rb <= 1e-12) & (r < 1.2 * R99),
    }
    print(f"\n=== {run} t={t:.1f}  R99={R99:.2f}")
    print(f"{'populacao':18s} {'n':>6s} {'cs p50':>7s} {'c0/cs':>6s} | {'|a_av| atual':>12s} "
          f"{'|a_av| c0':>10s} {'|a_drag|':>9s} {'|a_mar|':>8s}")
    for nome, m in pops.items():
        idx = np.where(m)[0]
        if len(idx) == 0:
            continue
        s = rng.choice(idx, min(3000, len(idx)), replace=False)
        A = a_visc_art(d, s)
        cs = np.median(d["cs"][s])
        ge = GAMMA + 1.5 * GAMMA * rb[s] ** 2
        drag = np.median(ge * np.hypot(d["u"][s], d["v"][s]))
        print(f"{nome:18s} {m.sum():6d} {cs:7.3f} {C0 / max(cs, 1e-9):6.1f} | "
              f"{np.median(A['atual']):12.2e} {np.median(A['c0']):10.2e} {drag:9.2e} "
              f"{np.median(d['au_mar'][s]):8.2e}")
    pares = np.array(list(cKDTree(np.c_[d["x"], d["y"]]).query_pairs(0.05 * DX)))
    if len(pares):
        u = np.unique(pares.ravel())
        print(f"pares < 0.05 dx: {len(pares)}  particulas {len(u)}  cs p50 {np.median(d['cs'][u]):.3f} "
              f"(c0/cs = {C0 / np.median(d['cs'][u]):.2f})  r/R99 p50 {np.median(r[u] / R99):.2f}")
