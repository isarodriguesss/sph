"""Compara rodadas em R99 IGUAL, nao em t igual.

Quando uma alavanca muda a velocidade da frente (B1-motor, licao #86), comparar no mesmo
instante confunde "outra forma" com "outro tamanho". Aqui cada rodada e amostrada no
primeiro frame em que R99 cruza cada alvo, e as metricas de forma (amplitude, dedos,
R99/R90) vem SEMPRE junto das de biologia (licao #76: continuidade e ocupacao premiam o disco).

Razao frente/proliferacao = 2*(dR/dt)/(R*k): expansao exigida pela area sobre a taxa
especifica de ganho de massa das vivas, k = G*r_growth*rho_b*(1-rho_b)*f(c_n).
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from plot_expansao_comp import serie  # noqa: E402

R_GROWTH = 0.02


def fcn(c):
    t = np.clip((c - 0.4) / 0.4, 0.0, 1.0)
    f = t * t * (3.0 - 2.0 * t)
    f[c < 0.4] = 0.0
    return f


def bio(f, G):
    with h5py.File(f, "r") as h:
        a = h["particles"]["fluid"]["arrays"]
        x, y, rb, isf, cn = (np.asarray(a[k]) for k in ("x", "y", "rho_b_grown", "is_filler", "c_n"))
        isw = np.asarray(a["is_wake"]) if "is_wake" in a else np.zeros_like(x)
    r = np.hypot(x, y)
    R99 = np.percentile(r[rb > 0.1], 99)
    viva = (isf < 0.5) & (rb > 0.1)
    cresce = viva & (rb < 0.8)
    k = np.median(G * R_GROWTH * rb[cresce] * (1 - rb[cresce]) * fcn(cn[cresce]))
    frente = r >= 0.85 * R99
    nv, nf = (viva & frente).sum(), ((isf > 0.5) & frente).sum()
    pin = ((rb >= 0.8) | (cn < 0.6)) & (isw < 0.5)
    return dict(k=k, viva_frente=int(nv), frac_frente=nv / max(nv + nf, 1),
                cn_viva=np.median(cn[viva]), pin_viva=np.mean(pin[viva]), n_viva=int(viva.sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--rotulos", nargs="+", required=True)
    ap.add_argument("--ganho", nargs="+", type=float, required=True, help="GROWTH_MASS_GAIN de cada run")
    ap.add_argument("--raios", nargs="+", type=float, default=[1.2, 1.6, 2.0, 2.4])
    args = ap.parse_args()

    print(f"{'run':>10} {'R99':>4} {'t':>5} {'dR/dt':>6} {'amp':>6} {'dedos':>5} {'R99/R90':>7} | "
          f"{'vivas':>5} {'k_med':>7} {'razao':>6} | {'vivas>.85R':>10} {'frac viva':>9} "
          f"{'c_n viva':>8} {'%pin':>5}")
    for run, rot, G in zip(args.runs, args.rotulos, args.ganho):
        S = serie(run)
        if len(S) == 0:
            continue
        fs = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
        ts = []
        for f in fs:
            with h5py.File(f, "r") as h:
                ts.append(float(h["solver_data"].attrs["t"]))
        ts = np.array(ts)
        t, R99, R90, amp, ded = S.T
        for Ra in args.raios:
            if R99.max() < Ra:
                print(f"{rot:>10} {Ra:4.1f}  ainda nao alcancou (R99 atual {R99[-1]:.2f}, t={t[-1]:.1f})")
                break
            i = int(np.argmax(R99 >= Ra))
            w = np.abs(t - t[i]) <= 4.0
            rate = np.polyfit(t[w], R99[w], 1)[0] if w.sum() >= 3 else np.nan
            j = slice(max(0, i - 1), i + 2)
            B = bio(fs[int(np.argmin(np.abs(ts - t[i])))], G)
            razao = 2 * rate / (R99[i] * B["k"]) if B["k"] > 0 else np.nan
            print(f"{rot:>10} {Ra:4.1f} {t[i]:5.1f} {rate:6.3f} {np.median(amp[j]):6.3f} "
                  f"{np.median(ded[j]):5.0f} {R99[i] / R90[i]:7.2f} | {B['n_viva']:5d} {B['k']:7.4f} "
                  f"{razao:5.1f}x | {B['viva_frente']:10d} {100 * B['frac_frente']:8.0f}% "
                  f"{B['cn_viva']:8.2f} {100 * B['pin_viva']:4.0f}%")


if __name__ == "__main__":
    main()
