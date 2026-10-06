"""Relacao de dispersao do modelo: sigma(m), a taxa de crescimento por modo azimutal.

    python tools/dispersao.py runs/swarm/P2R23_conduz [--out plots/dispersao.png]

Mede R(theta) direto das PARTICULAS (nao do raster): colonia = rho_b>=0.1 ou filler (§2.2),
R(theta) = maior raio com colonia no bin angular. Espectro por rfft, amplitude RELATIVA
a_m = |R_m| / R_medio. sigma_m por ajuste de ln(a_m) contra t na janela pedida.

Reporta junto os dois cortes de comprimento de onda curto, que sao o ponto da medicao:
  m_dx    = pi*R/dx          — Nyquist das particulas na borda
  m_port  = N_portadores/2   — amostragem pelas celulas VIVAS na frente (licoes #98/#106)
"""
import argparse
import glob
import os

import h5py
import numpy as np

NTH = 720


def frame(fn):
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        d = {k: a[k][:] for k in ("x", "y", "rho_b_grown", "is_filler")}
        t = float(f.attrs["solver_data"]["t"]) if "solver_data" in f.attrs else None
        if t is None:
            t = float(f["solver_data"].attrs["t"])
    return d, t


def r_theta(d, nth=NTH):
    col = (d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)
    x, y = d["x"][col], d["y"][col]
    r = np.hypot(x, y)
    th = np.arctan2(y, x)
    ib = ((th + np.pi) / (2 * np.pi) * nth).astype(int) % nth
    R = np.zeros(nth)
    np.maximum.at(R, ib, r)
    vaz = R == 0.0
    if vaz.any() and not vaz.all():  # bin sem colonia: interpola circularmente
        idx = np.arange(nth)
        R[vaz] = np.interp(idx[vaz], idx[~vaz], R[~vaz], period=nth)
    return R, int(col.sum())


def n_portadores(d, R, frac=0.85):
    """Celulas VIVAS (nao filler, acima do quorum) na frente — o gargalo da licao #98."""
    viva = (d["rho_b_grown"] >= 0.1) & (d["is_filler"] < 0.5)
    rr = np.hypot(d["x"], d["y"])
    return int(np.sum(viva & (rr >= frac * np.median(R))))


def espectro(R):
    A = np.abs(np.fft.rfft(R)) / len(R)
    A[1:] *= 2.0
    return A / A[0]  # amplitude relativa ao raio medio


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--out", default="plots/dispersao.png")
    p.add_argument("--mmax", type=int, default=40)
    p.add_argument("--t0", type=float, default=0.0)
    p.add_argument("--t1", type=float, default=25.0)
    p.add_argument("--dx", type=float, default=14.0 / 260.0)
    args = p.parse_args()

    fs = sorted(glob.glob(os.path.join(args.run, "main_output", "*.hdf5")))
    T, A, RM, NP, NC = [], [], [], [], []
    for fn in fs:
        d, t = frame(fn)
        R, ncol = r_theta(d)
        if R.max() == 0:
            continue
        T.append(t); A.append(espectro(R)); RM.append(R.mean())
        NP.append(n_portadores(d, R)); NC.append(ncol)
    T = np.array(T); A = np.array(A); RM = np.array(RM)
    NP = np.array(NP); NC = np.array(NC)

    print("run:", args.run)
    print("%6s %8s %9s %9s %8s %8s" % ("t", "R_medio", "n_col", "n_viva_frente", "m_dx", "m_port"))
    for i in range(len(T)):
        print("%6.2f %8.3f %9d %13d %8.0f %8.1f"
              % (T[i], RM[i], NC[i], NP[i], np.pi * RM[i] / args.dx, NP[i] / 2.0))

    jan = (T >= args.t0) & (T <= args.t1)
    tt = T[jan]
    ms = np.arange(1, args.mmax + 1)
    sig = np.full(len(ms), np.nan)
    for k, m in enumerate(ms):
        y = A[jan, m]
        ok = y > 0
        if ok.sum() >= 3:
            sig[k] = np.polyfit(tt[ok], np.log(y[ok]), 1)[0]

    print("\njanela do ajuste: t em [%.1f, %.1f], %d instantes" % (tt.min(), tt.max(), len(tt)))
    print("\n%4s %10s %10s %10s" % ("m", "a_m(t0)", "a_m(t1)", "sigma_m"))
    for k, m in enumerate(ms):
        if m <= 30:
            print("%4d %10.2e %10.2e %10.4f" % (m, A[jan][0, m], A[jan][-1, m], sig[k]))

    fin = np.isfinite(sig)
    if fin.any():
        kmax = int(np.nanargmax(sig))
        print("\nmaximo de sigma(m): m* = %d, sigma = %.4f" % (ms[kmax], sig[kmax]))
        print("sigma(m) e monotonico crescente ate m=%d? %s"
              % (args.mmax, bool(np.all(np.diff(sig[fin]) > -1e-9))))

    np.savez(os.path.splitext(args.out)[0] + ".npz",
             T=T, A=A, RM=RM, NP=NP, NC=NC, ms=ms, sig=sig)
    print("\ndados em", os.path.splitext(args.out)[0] + ".npz")


if __name__ == "__main__":
    main()
