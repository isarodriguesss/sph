"""Expansao sem obstaculo x com os obstaculos do R1: a MESMA expansao, organizacao diferente.

Por que nao duas curvas de R99(t): a metrica primaria e media azimutal e e CEGA ao efeito
principal (CRITERIOS_RUGOSIDADE.md §12) — no R1 a modulacao vale 22% e a media difere 3%.
Uma figura com so as duas curvas mentiria por omissao. Entao a figura mostra o nulo E o
efeito no mesmo quadro:

  (a),(b) a FORMA em t=50, com a rede desenhada nos DOIS paineis — tenue em (a), onde ela
          nao existe, como referencia geometrica.
  (c)     R99(t) das duas — quase sobrepostas. E o nulo, e aparece honestamente.
  (d)     os modos da REDE ao longo do tempo. Em t=50 o travamento de m=6 ja decaiu (0.223
          em t=25.6 -> 0.013), mas a assinatura NAO sumiu: ela migrou para o harmonico
          m=12 = 2x6, que vale 0.200 no R1 contra 0.064 no liso. E o modo PROPRIO do
          modelo (m=10-11 no liso, licao #107) desaparece do R1. Por isso (d) e serie
          temporal e nao espectro num instante: a historia e a troca de m=6 por m=12.

    python plots/fig_obstaculo.py                    # usa os caminhos padrao
    python plots/fig_obstaculo.py --rmed 1.50 --out plots/fig_obstaculo.png
"""

import argparse
import glob
import os
import sys

import h5py
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Circle

sys.path.insert(0, "plots")
import fig_tese as FT  # noqa: E402

DX = 14.0 / 260.0
NTH = 720
LISO = "runs/swarm/P2R23_conduz"
RUG = "runs/rugosidade/R1_lambda28"
COR_LISO = "#1b7a98"
COR_RUG = "#9a9855"


def r_theta(d):
    col = (d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)
    r = np.hypot(d["x"][col], d["y"][col])
    th = np.arctan2(d["y"][col], d["x"][col])
    ib = ((th + np.pi) / (2 * np.pi) * NTH).astype(int) % NTH
    R = np.zeros(NTH)
    np.maximum.at(R, ib, r)
    v = R == 0.0
    if v.any() and not v.all():
        i = np.arange(NTH)
        R[v] = np.interp(i[v], i[~v], R[~v], period=NTH)
    return np.roll(R, -NTH // 2)


def serie(run):
    out = []
    for f in sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5"))):
        with h5py.File(f, "r") as h:
            t = float(h["solver_data"].attrs["t"])
        d = FT.carrega(f)
        R = r_theta(d)
        out.append(dict(f=f, t=t, R=R, rmed=R.mean(),
                        r99=float(np.percentile(np.hypot(
                            d["x"][(d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)],
                            d["y"][(d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)]), 99))))
    return out


def am(R, m):
    F = np.fft.rfft(R)
    return 2 * np.abs(F[m]) / np.abs(F[0])


def a6(R):
    F = np.fft.rfft(R)
    return (2 * np.abs(F[6]) / np.abs(F[0]),
            np.degrees(-np.angle(F[6]) / 6) % 60.0)


def painel(ax, q, pil, L, n, iso, titulo, tenue):
    g, img, _, _ = FT.campo(FT.carrega(q["f"]), L, n, False, 0.0, agar_iso=iso)
    ax.imshow(img, origin="lower", extent=[-L, L, -L, L], cmap=FT.CMAP_COL,
              vmin=0, vmax=1, interpolation="bilinear")
    if pil is not None:
        C, rp = pil
        dentro = np.hypot(C[:, 0], C[:, 1]) < L * 1.45
        for c in C[dentro]:
            ax.add_patch(Circle(c, rp, fill=not tenue,
                                facecolor="#3a3a3a" if not tenue else "none",
                                edgecolor="#6b6b6b", lw=0.5,
                                ls=":" if tenue else "-", alpha=0.45 if tenue else 1.0,
                                zorder=3))
    for k in range(6):  # direcoes das GARGANTAS da 1a coroa do R1
        a = np.radians(30 + 60 * k)
        ax.plot([0, L * 1.42 * np.cos(a)], [0, L * 1.42 * np.sin(a)],
                color=FT.COR_BORDA, lw=0.6, ls=(0, (6, 4)), alpha=0.5, zorder=4)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(titulo, fontsize=10, pad=6)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--liso", default=LISO)
    p.add_argument("--rug", default=RUG)
    p.add_argument("--t", type=float, default=50.0)
    p.add_argument("--n", type=int, default=600)
    p.add_argument("--iso", type=float, default=2.0)
    p.add_argument("--out", default="plots/fig_obstaculo.png")
    a = p.parse_args()

    SL, SR = serie(a.liso), serie(a.rug)
    ql = min(SL, key=lambda q: abs(q["t"] - a.t))
    qr = min(SR, key=lambda q: abs(q["t"] - a.t))
    pil = FT.carrega_pilares(qr["f"])
    L = 4.9
    al, fl = a6(ql["R"]); ar, fr = a6(qr["R"])
    print(f"liso: t={ql['t']:.1f} R_med={ql['rmed']:.3f} a6={al:.3f} fase={fl:.1f}")
    print(f"R1  : t={qr['t']:.1f} R_med={qr['rmed']:.3f} a6={ar:.3f} fase={fr:.1f}")

    fig, ax = plt.subplots(2, 2, figsize=(9.6, 9.0))
    painel(ax[0, 0], ql, pil, L, a.n, a.iso,
           f"(a) sem obstaculo — t = {ql['t']:.0f} s\n"
           f"$R_{{99}}$ = {ql['r99']:.2f}", tenue=True)
    painel(ax[0, 1], qr, pil, L, a.n, a.iso,
           f"(b) com obstaculos (R1, $\\Lambda$ = 28 dx) — t = {qr['t']:.0f} s\n"
           f"$R_{{99}}$ = {qr['r99']:.2f}  (−5.3%)", tenue=False)
    ax[0, 0].text(0.02, 0.02, "rede desenhada como referencia\n(ausente nesta rodada)",
                  transform=ax[0, 0].transAxes, fontsize=7, color="#555",
                  va="bottom", ha="left")

    b = ax[1, 0]
    for S, cor, lbl in ((SL, COR_LISO, "sem obstaculo"), (SR, COR_RUG, "com obstaculos (R1)")):
        b.plot([q["t"] for q in S], [q["r99"] for q in S], "-o", color=cor, ms=3.5,
               lw=1.6, label=lbl)
    b.axvline(qr["t"], color="#999", lw=0.8, ls=":")
    b.set_xlabel("t (s)"); b.set_ylabel(r"$R_{99}$")
    b.set_title("(c) expansao — media azimutal:\n"
                r"$dR/dt$ 0.0756 $\to$ 0.0733 (−3.1%)", fontsize=10, pad=6)
    b.legend(fontsize=8, frameon=False, loc="upper left")
    b.grid(alpha=0.25, lw=0.5)

    c = ax[1, 1]
    for S, cor, lbl in ((SL, COR_LISO, "sem obstaculo"), (SR, COR_RUG, "com obstaculos (R1)")):
        tt = [q["t"] for q in S]
        c.plot(tt, [a6(q["R"])[0] for q in S], "-o", color=cor, ms=3.5, lw=1.8,
               label=f"{lbl} — $m$=6")
        c.plot(tt, [am(q["R"], 12) for q in S], "--s", color=cor, ms=3.0, lw=1.3,
               alpha=0.75, label=f"{lbl} — $m$=12")
    c.axvline(qr["t"], color="#999", lw=0.8, ls=":")
    c.annotate("frente atravessa a 1a coroa", xy=(25.6, 0.223), xytext=(31.5, 0.205),
               fontsize=7.5, color="#6b6b45",
               arrowprops=dict(arrowstyle="->", color="#6b6b45", lw=0.8))
    c.set_ylim(-0.012, 0.255)
    c.set_xlabel("t (s)")
    c.set_ylabel("amplitude relativa do modo")
    c.set_title("(d) os modos da REDE (6 e seu harmonico 12)\n"
                "o pico de $m$=6 migra para $m$=12", fontsize=10, pad=6)
    c.legend(fontsize=7, frameon=False, loc="upper left", ncol=1)
    c.grid(alpha=0.25, lw=0.5)

    fig.suptitle("A mesma expansao, com organizacao diferente", fontsize=12.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(a.out, dpi=190)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
