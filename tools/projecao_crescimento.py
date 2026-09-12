"""Saida C (licao #86): velocidade de EXPANSAO POR CRESCIMENTO, por projecao.

Numa colonia incompressivel que so ganha volume por proliferacao, `div u = Gamma`, com
`Gamma = d ln m/dt = G*r_growth*rho_b*(1-rho_b)*f(c_n)` nas vivas que crescem. Escrevendo
`u = grad(phi)`, isso vira a equacao de Poisson `lap(phi) = Gamma` na colonia, com `phi = 0`
na superficie livre (o agar, sem pressao). E a etapa de projecao do ISPH de Cummins &
Rudman 1999 (Violeau §5.4) restrita ao termo-fonte de crescimento.

Laplaciano de Brookshaw (o mesmo operador da `SurfactantEquation`) montado como matriz
esparsa; gradiente com KGC de Bonet-Lok (Liu §3.3). O agar entra como no de Dirichlet.

Uso:  --teste   valida contra o disco analitico (u = Gamma*r/2)
      RUN --t ...  mede u_g num run e compara com o movimento real
"""

import argparse
import glob
import os

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree

R_GROWTH = 0.02


def dw(r, h):
    s = 10.0 / (7.0 * np.pi * h * h)
    q = r / h
    d = np.where(q < 1.0, -3.0 * q + 2.25 * q * q, np.where(q < 2.0, -0.75 * (2.0 - q) ** 2, 0.0))
    return s * d / h


def fcn(c):
    t = np.clip((c - 0.4) / 0.4, 0.0, 1.0)
    f = t * t * (3.0 - 2.0 * t)
    f[c < 0.4] = 0.0
    return f


def projecao(x, y, V, h, colonia, gamma):
    """Resolve lap(phi)=gamma na colonia (phi=0 fora). Devolve phi e u=grad(phi)."""
    P = np.c_[x, y]
    T = cKDTree(P)
    idx_col = np.where(colonia)[0]
    viz = T.query_ball_point(P[idx_col], 2.0 * h)
    ii = np.concatenate([np.full(len(v), k) for k, v in enumerate(viz)])
    jj = np.concatenate(viz).astype(int)
    di = idx_col[ii]
    m = jj != di
    ii, jj, di = ii[m], jj[m], di[m]

    xij, yij = x[di] - x[jj], y[di] - y[jj]
    r = np.hypot(xij, yij)
    g = dw(r, h) / np.maximum(r, 1e-12)
    DX, DY = g * xij, g * yij
    c = 2.0 * V[jj] * (xij * DX + yij * DY) / (r * r + 0.01 * h * h)  # < 0

    n = len(idx_col)
    loc = -np.ones(len(x), dtype=int)
    loc[idx_col] = np.arange(n)
    diag = np.bincount(ii, weights=c, minlength=n)
    col_j = loc[jj]
    viz_col = col_j >= 0
    rows = np.r_[np.arange(n), ii[viz_col]]
    cols = np.r_[np.arange(n), col_j[viz_col]]
    vals = np.r_[diag, -c[viz_col]]
    A = csr_matrix((vals, (rows, cols)), shape=(n, n))
    phi_c = spsolve(A, gamma[idx_col])
    phi = np.zeros(len(x))
    phi[idx_col] = phi_c

    # gradiente KGC
    acc = lambda w: np.bincount(ii, weights=w, minlength=n)  # noqa: E731
    Vj = V[jj]
    Mxx, Mxy = acc(-Vj * xij * DX), acc(-Vj * xij * DY)
    Myx, Myy = acc(-Vj * yij * DX), acc(-Vj * yij * DY)
    det = Mxx * Myy - Mxy * Myx
    ok = np.abs(det) > 0.25
    sd = np.where(ok, det, 1.0)
    Lxx, Lxy = np.where(ok, Myy / sd, 1.0), np.where(ok, -Mxy / sd, 0.0)
    Lyx, Lyy = np.where(ok, -Myx / sd, 0.0), np.where(ok, Mxx / sd, 1.0)
    dphi = phi[jj] - phi[di]
    gx, gy = acc(Vj * dphi * DX), acc(Vj * dphi * DY)
    ux = np.zeros(len(x))
    uy = np.zeros(len(x))
    ux[idx_col] = Lxx * gx + Lxy * gy
    uy[idx_col] = Lyx * gx + Lyy * gy
    return phi, ux, uy


def teste():
    dx = 0.0538
    h = 1.8 * dx
    g = np.arange(-2.0, 2.0 + dx / 2, dx)
    X, Y = np.meshgrid(g, g)
    x, y = X.ravel(), Y.ravel()
    r = np.hypot(x, y)
    R0, G0 = 1.0, 0.01
    col = r < R0
    gam = np.where(col, G0, 0.0)
    V = np.full(len(x), dx * dx)
    phi, ux, uy = projecao(x, y, V, h, col, gam)
    ur = (ux * x + uy * y) / np.maximum(r, 1e-12)
    print("disco R0=1, Gamma=0.01 uniforme  (exato: u_r = Gamma*r/2)")
    for lo, hi in ((0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 0.95)):
        b = col & (r >= lo) & (r < hi)
        rm = np.median(r[b])
        print(f"  r~{rm:.2f}: u_r = {np.median(ur[b]):.5f}  exato {G0 * rm / 2:.5f}  "
              f"erro {100 * (np.median(ur[b]) / (G0 * rm / 2) - 1):+.1f}%")


def frame(run, alvo):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
    ts = []
    for f in fs:
        with h5py.File(f, "r") as hh:
            ts.append(float(hh["solver_data"].attrs["t"]))
    return fs[int(np.argmin(np.abs(np.array(ts) - alvo)))]


def medir(run, alvos, G):
    for alvo in alvos:
        f = frame(run, alvo)
        with h5py.File(f, "r") as hh:
            t = float(hh["solver_data"].attrs["t"])
            a = hh["particles"]["fluid"]["arrays"]
            d = {k: np.asarray(a[k], dtype=float) for k in
                 ("x", "y", "u", "v", "m", "rho", "h", "rho_b_grown", "is_filler", "c_n")}
            d["is_env"] = np.asarray(a["is_env"], dtype=float) if "is_env" in a else np.zeros_like(d["x"])
        x, y, rb, isf = d["x"], d["y"], d["rho_b_grown"], d["is_filler"]
        r = np.hypot(x, y)
        R99 = np.percentile(r[rb > 0.1], 99)
        colonia = ((rb > 0.0) | (isf > 0.5)) & (r < 1.3 * R99)
        cresce = (isf < 0.5) & (d["is_env"] < 0.5) & (rb > 1e-12) & (rb < 0.8)
        gam = np.where(cresce, G * R_GROWTH * rb * (1 - rb) * fcn(d["c_n"]), 0.0)
        V = d["m"] / d["rho"]
        phi, ux, uy = projecao(x, y, V, float(np.median(d["h"])), colonia, gam)
        ex, ey = x / np.maximum(r, 1e-12), y / np.maximum(r, 1e-12)
        ug = ux * ex + uy * ey
        ur = d["u"] * ex + d["v"] * ey
        viva = (isf < 0.5) & (rb > 0.1)
        print(f"\n=== {run}  t={t:.1f}  R99={R99:.2f}  colonia={colonia.sum()}  "
              f"int(Gamma dV)={np.sum(gam * V):.4f}  area colonia={np.sum(V[colonia]):.3f}")
        print(f"{'r/R99':>9} {'n viva':>6} {'Gamma med':>9} | {'u_g (C)':>9} {'u_r real':>9} {'razao':>6}")
        for lo, hi in ((0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.0), (1.0, 1.3)):
            b = viva & (r >= lo * R99) & (r < hi * R99)
            if b.sum() < 3:
                continue
            print(f"{lo:4.2f}-{hi:4.2f} {b.sum():6d} {np.median(gam[b]):9.5f} | "
                  f"{np.median(ug[b]):9.5f} {np.median(ur[b]):9.5f} "
                  f"{np.median(ur[b]) / max(np.median(ug[b]), 1e-12):5.1f}x")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?")
    ap.add_argument("--t", type=float, nargs="+", default=[20.0, 50.0, 80.0])
    ap.add_argument("--ganho", type=float, default=1.0)
    ap.add_argument("--teste", action="store_true")
    args = ap.parse_args()
    if args.teste:
        teste()
    if args.run:
        medir(args.run, args.t, args.ganho)


if __name__ == "__main__":
    main()
