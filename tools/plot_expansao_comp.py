"""Trajetoria da expansao comparando rodadas.

`R99` sozinho nao decide nada (licao #67-J: o p99 pode ser fixado por poucas espiculas, e
#76: toda metrica interna e maximizada pelo disco). Entao o painel traz SEMPRE, lado a
lado: raio, ritmo, e as duas metricas de FORMA que o blob reprova — dedos e amplitude.
"""

import argparse
import glob
import os

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

CORES = ["#888888", "#2e86de", "#e17055", "#00b894", "#d63031", "#6c5ce7"]


def serie(run, nb=128):
    out = []
    for f in sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5"))):
        with h5py.File(f, "r") as h:
            t = float(h["solver_data"].attrs["t"])
            a = h["particles"]["fluid"]["arrays"]
            x, y, rb, isf = (np.asarray(a[k]) for k in ("x", "y", "rho_b_grown", "is_filler"))
        port = rb > 0.1
        if port.sum() < 20:
            continue
        r = np.hypot(x, y)
        th = np.arctan2(y, x)
        corpo = (rb >= 0.1) | (isf > 0.5)
        R99 = float(np.percentile(r[port], 99))
        R90 = float(np.percentile(r[port], 90))
        bins = np.linspace(-np.pi, np.pi, nb + 1)
        ib = np.digitize(th[corpo], bins) - 1
        rc = r[corpo]
        Rt = np.zeros(nb)
        for k in range(nb):
            m = ib == k
            if m.sum():
                Rt[k] = rc[m].max()
        amp = Rt.std() / max(Rt.mean(), 1e-12)
        ded = int(((Rt > np.roll(Rt, 1)) & (Rt > np.roll(Rt, -1))).sum())
        out.append((t, R99, R90, amp, ded))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="+", required=True)
    ap.add_argument("--rotulos", nargs="+", required=True)
    ap.add_argument("--tmax", type=float, default=40.0)
    ap.add_argument("--out", default="plots/expansao_comp.png")
    args = ap.parse_args()

    fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))
    for i, (run, rot) in enumerate(zip(args.runs, args.rotulos)):
        S = serie(run)
        S = S[S[:, 0] <= args.tmax]
        c = CORES[i % len(CORES)]
        t, R99, R90, amp, ded = S.T
        ax[0, 0].plot(t, R99, "-o", color=c, ms=3, label=rot)
        ax[0, 0].plot(t, R90, "--", color=c, lw=0.8, alpha=0.6)
        # ritmo por ajuste linear numa janela deslizante de ~8 s: np.gradient nos
        # extremos produz artefato de borda que le como pico de expansao
        tm, rate = [], []
        for k in range(len(t)):
            w = np.abs(t - t[k]) <= 4.0
            if w.sum() >= 3:
                tm.append(t[k])
                rate.append(np.polyfit(t[w], R99[w], 1)[0])
        ax[0, 1].plot(tm, rate, "-", color=c, lw=1.6, label=rot)
        # dedos: mediana movel de 3, com o bruto em transparencia
        med = np.array([np.median(ded[max(0, k - 1) : k + 2]) for k in range(len(ded))])
        ax[1, 0].plot(t, ded, "-", color=c, lw=0.6, alpha=0.3)
        ax[1, 0].plot(t, med, "-o", color=c, ms=3, label=rot)
        ax[1, 1].plot(t, amp, "-o", color=c, ms=3, label=rot)

    for a, tit, yl in (
        (ax[0, 0], "RAIO — cheio R99, tracejado R90 (licao #67-J: nunca ler R99 sozinho)", "R"),
        (ax[0, 1], "RITMO  dR99/dt", "dR/dt"),
        (ax[1, 0], "FORMA — numero de dedos", "dedos"),
        (ax[1, 1], "FORMA — amplitude  std(R)/mean(R)", "amplitude"),
    ):
        a.set_title(tit, fontsize=10)
        a.set_xlabel("t")
        a.set_ylabel(yl)
        a.grid(alpha=0.25)
        a.legend(fontsize=8)
    fig.suptitle(
        "Trajetoria da expansao — raio e ritmo SEMPRE ao lado das metricas de forma "
        "(licao #76: continuidade e ocupacao sao maximizadas pelo disco)",
        fontsize=11,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(args.out, dpi=115)
    print("escrito:", args.out)


if __name__ == "__main__":
    main()
