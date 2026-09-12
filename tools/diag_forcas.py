"""Balanco de forcas por particula, reconstruido offline a partir do HDF5.

Pergunta: por que as vivas do meio do corpo andam 5-10x mais devagar que as da frente?
Cada termo do `equations_main` e recalculado com a MESMA formula do solver (CubicSpline 2D,
h_ij medio, KGC de Bonet-Lok na Marangoni) e projetado na direcao radial. `au_flag` e
`au_mar` gravados no HDF5 servem de validacao da reconstrucao; o fechamento
`sum(forcas) ~ gamma_eff*u_r` valida o conjunto.
"""

import argparse
import glob
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

MU, GAMMA, BETA, F0, C0, ALPHA = 0.020, 60.0, 5.0, 3.0, 0.35, 0.12
FLAG_LO, FLAG_HI = 0.2, 0.6
MAR_LO, MAR_HI = 0.15, 0.5
KGC_DET_MIN = 0.25
RHO0 = 1.0


def dw(r, h):
    """dW/dr do CubicSpline 2D do PySPH."""
    s = 10.0 / (7.0 * np.pi * h * h)
    q = r / h
    d = np.where(q < 1.0, -3.0 * q + 2.25 * q * q, np.where(q < 2.0, -0.75 * (2.0 - q) ** 2, 0.0))
    return s * d / h


def smooth(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def carrega(run, alvo):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
    ts = []
    for f in fs:
        with h5py.File(f, "r") as h:
            ts.append(float(h["solver_data"].attrs["t"]))
    f = fs[int(np.argmin(np.abs(np.array(ts) - alvo)))]
    with h5py.File(f, "r") as h:
        t = float(h["solver_data"].attrs["t"])
        a = h["particles"]["fluid"]["arrays"]
        ks = ("x", "y", "u", "v", "m", "rho", "p", "h", "cs", "c_n", "rho_b_grown",
              "is_filler", "is_matrix", "is_env", "is_wake", "au_flag", "au_mar")
        d = {k: np.asarray(a[k], dtype=float) for k in ks}
    return t, d


def forcas(d, sel):
    """Aceleracoes em `sel` (indices), vizinhos em todo o fluido."""
    x, y = d["x"], d["y"]
    T = cKDTree(np.c_[x, y])
    H = float(np.median(d["h"]))
    viz = T.query_ball_point(np.c_[x[sel], y[sel]], 2.0 * H)
    ii = np.concatenate([np.full(len(v), k) for k, v in enumerate(viz)])
    jj = np.concatenate(viz).astype(int)
    di = sel[ii]
    keep = jj != di
    ii, jj, di = ii[keep], jj[keep], di[keep]

    xij, yij = x[di] - x[jj], y[di] - y[jj]
    r = np.hypot(xij, yij)
    hij = 0.5 * (d["h"][di] + d["h"][jj])
    g = dw(r, hij) / np.maximum(r, 1e-12)
    DX, DY = g * xij, g * yij                      # grad_i W_ij
    mj, rj, ri = d["m"][jj], d["rho"][jj], d["rho"][di]
    vol = mj / rj
    n = len(sel)
    acc = lambda w: np.bincount(ii, weights=w, minlength=n)  # noqa: E731

    # pressao (MomentumEquation) + viscosidade artificial de Monaghan
    pi, pj = d["p"][di], d["p"][jj]
    fp = -mj * (pi / ri**2 + pj / rj**2)
    uij, vij = d["u"][di] - d["u"][jj], d["v"][di] - d["v"][jj]
    vdotr = uij * xij + vij * yij
    muij = hij * vdotr / (r * r + 0.01 * hij * hij)
    # o PySPH usa `cij = 0.5*(cs_i + cs_j)` e `cs` aqui e o SURFACTANTE (colisao de nome)
    cij = 0.5 * (d["cs"][di] + d["cs"][jj])
    piij = np.where(vdotr < 0.0, -ALPHA * cij * muij / (0.5 * (ri + rj)), 0.0)
    ap = (acc(fp * DX), acc(fp * DY))
    aav = (acc(-mj * piij * DX), acc(-mj * piij * DY))

    # viscosidade fisica (ViscousForce)
    dot = xij * DX + yij * DY
    mue = MU * np.minimum(0.5 * (ri + rj), 1.0)
    cf = 2.0 * mue * vol * dot / (r * r + 0.01 * d["h"][di] ** 2)
    avis = (acc(cf * uij), acc(cf * vij))

    # gradiente de rho_b (BiomassGradient) -> gate da Marangoni
    rb = np.where(d["is_env"] > 0.5, 0.0, d["rho_b_grown"])
    drb = rb[jj] - rb[di]
    grb = np.hypot(acc(vol * drb * DX), acc(vol * drb * DY))

    # KGC (Bonet-Lok): M = sum V (x_j - x_i) (x) grad W ; L = M^-1
    Mxx, Mxy = acc(-vol * xij * DX), acc(-vol * xij * DY)
    Myx, Myy = acc(-vol * yij * DX), acc(-vol * yij * DY)
    det = Mxx * Myy - Mxy * Myx
    ok = np.abs(det) >= KGC_DET_MIN
    sd = np.where(ok, det, 1.0)
    Lxx = np.where(ok, Myy / sd, 1.0)
    Lxy = np.where(ok, -Mxy / sd, 0.0)
    Lyx = np.where(ok, -Myx / sd, 0.0)
    Lyy = np.where(ok, Mxx / sd, 1.0)

    # Marangoni: fonte/destino fora se filler nao-matriz
    fonte = (d["is_filler"][jj] < 0.5) | (d["is_matrix"][jj] > 0.5)
    dest = (d["is_filler"][sel] < 0.5) | (d["is_matrix"][sel] > 0.5)
    dcs = (d["cs"][jj] - d["cs"][di]) * fonte
    gcx, gcy = acc(vol * dcs * DX), acc(vol * dcs * DY)
    gate_m = np.where(grb >= MAR_LO, smooth((grb - MAR_LO) / (MAR_HI - MAR_LO)), 0.0) * dest
    ncx = Lxx * gcx + Lxy * gcy
    ncy = Lyx * gcx + Lyy * gcy
    amar = (-BETA * gate_m * ncx, -BETA * gate_m * ncy)

    # flagelar (gradiente SEM KGC, como no solver)
    rbs = d["rho_b_grown"][sel]
    mid = 0.5 * (FLAG_LO + FLAG_HI)
    tt = np.where(rbs < mid, (rbs - FLAG_LO) / (mid - FLAG_LO), (FLAG_HI - rbs) / (FLAG_HI - mid))
    gate_f = np.where((rbs >= FLAG_LO) & (rbs <= FLAG_HI), smooth(tt), 0.0)
    gate_f *= d["is_filler"][sel] < 0.5
    mag = np.hypot(gcx, gcy) + 1e-9
    aflag = (-F0 * gate_f * gcx / mag, -F0 * gate_f * gcy / mag)

    # arrasto
    ge = GAMMA + 1.5 * GAMMA * rbs**2
    adrag = (-ge * d["u"][sel], -ge * d["v"][sel])

    return dict(p=ap, av=aav, vis=avis, mar=amar, flag=aflag, drag=adrag,
                gate_m=gate_m, gate_f=gate_f, grb=grb, ge=ge)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--t", type=float, nargs="+", default=[20.0, 40.0])
    args = ap.parse_args()

    for alvo in args.t:
        t, d = carrega(args.run, alvo)
        x, y, rb = d["x"], d["y"], d["rho_b_grown"]
        r = np.hypot(x, y)
        R99 = np.percentile(r[rb > 0.1], 99)
        viva = (d["is_filler"] < 0.5) & (rb > 0.1)
        pin = ((rb >= 0.8) | (d["c_n"] < 0.6) | (d["is_filler"] > 0.5)) & (d["is_wake"] < 0.5)
        sel = np.where(viva & (r < 1.3 * R99))[0]
        F = forcas(d, sel)
        ex, ey = x[sel] / np.maximum(r[sel], 1e-12), y[sel] / np.maximum(r[sel], 1e-12)
        rad = {k: F[k][0] * ex + F[k][1] * ey for k in ("p", "av", "vis", "mar", "flag", "drag")}
        tot = F["p"][0] * 0
        motor = rad["mar"] + rad["flag"]
        soma = sum(rad[k] for k in ("p", "av", "vis", "mar", "flag"))
        ur = d["u"][sel] * ex + d["v"][sel] * ey
        _ = tot

        # validacao contra o que o solver gravou
        mm = np.hypot(*F["mar"])
        mf = np.hypot(*F["flag"])
        vm = np.corrcoef(mm, d["au_mar"][sel])[0, 1]
        vf = np.corrcoef(mf, d["au_flag"][sel])[0, 1]
        print(f"\n=== t={t:.1f}  R99={R99:.2f}  vivas analisadas={len(sel)}")
        print(f"validacao: corr(|a_mar| recalc, au_mar)={vm:.3f}  mediana {np.median(mm):.3f} vs "
              f"{np.median(d['au_mar'][sel]):.3f} | corr(|a_flag|, au_flag)={vf:.3f}  "
              f"mediana {np.median(mf):.3f} vs {np.median(d['au_flag'][sel]):.3f}")

        print(f"{'r/R99':>9} {'n':>5} {'%pin':>5} {'rho_b':>5} {'gateF>0':>7} {'gateM>0':>7} "
              f"{'|grb|':>6} | {'a_mar':>7} {'a_flag':>7} {'a_p':>8} {'a_av':>8} {'a_vis':>8} | "
              f"{'soma':>7} {'g*u_r':>7} {'u_r':>8}")
        for lo, hi in ((0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.0), (1.0, 1.3)):
            b = (r[sel] >= lo * R99) & (r[sel] < hi * R99)
            if b.sum() < 3:
                continue
            md = lambda a: np.median(a[b])  # noqa: E731
            livre = b & ~pin[sel]
            print(f"{lo:4.2f}-{hi:4.2f} {b.sum():5d} {100 * np.mean(pin[sel][b]):4.0f}% "
                  f"{md(rb[sel]):5.2f} {100 * np.mean(F['gate_f'][b] > 0):6.0f}% "
                  f"{100 * np.mean(F['gate_m'][b] > 0):6.0f}% {md(F['grb']):6.2f} | "
                  f"{md(rad['mar']):7.3f} {md(rad['flag']):7.3f} {md(rad['p']):8.2e} "
                  f"{md(rad['av']):8.2e} {md(rad['vis']):8.2e} | {md(soma):7.3f} "
                  f"{np.median((F['ge'] * ur)[b]):7.3f} {np.median(ur[livre]) if livre.sum() else np.nan:8.2e}")

        # decomposicao por gate flagelar (a hipotese): quem esta fora do gate anda devagar?
        print("  por classe de gate (bandas 0.3-0.85 R99, nao pinadas):")
        meio = (r[sel] >= 0.3 * R99) & (r[sel] < 0.85 * R99) & ~pin[sel]
        for nome, m in (("gate flagelar = 0", meio & (F["gate_f"] == 0)),
                        ("gate flagelar > 0", meio & (F["gate_f"] > 0)),
                        ("gate Marangoni = 0", meio & (F["gate_m"] == 0)),
                        ("gate Marangoni > 0", meio & (F["gate_m"] > 0))):
            if m.sum() < 3:
                print(f"    {nome:20s} n={m.sum()}")
                continue
            print(f"    {nome:20s} n={m.sum():5d} ({100 * m.sum() / meio.sum():4.0f}%)  "
                  f"rho_b {np.median(rb[sel][m]):.2f}  a_motor_r {np.median(motor[m]):7.3f}  "
                  f"u_r {np.median(ur[m]):8.2e}")


if __name__ == "__main__":
    main()
