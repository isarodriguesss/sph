"""Expansao sem obstaculo x com os obstaculos do R1: a MESMA expansao, organizacao diferente.

Por que nao duas curvas de R99(t): a metrica primaria e media azimutal e e CEGA ao efeito
principal (CRITERIOS_RUGOSIDADE.md §12) — no R1 a modulacao vale 22% e a media difere 3%.
Uma figura com so as duas curvas mentiria por omissao. Entao a figura mostra o nulo E o
efeito no mesmo quadro:

  (a),(b) a FORMA, em R_med CASADO (o tamanho nao e confundidor), com a rede desenhada nos
          DOIS paineis — tenue em (a), onde ela nao existe, como referencia geometrica. E
          isso que torna o travamento de fase visivel em vez de espectral: os lobos do R1
          caem nas GARGANTAS (30/90/150 graus), os do liso nao tem relacao com elas.
  (c)     R99(t) das duas — quase sobrepostas. E o nulo, e aparece honestamente.
  (d)     R(theta)/R_medio desenrolado no instante do pico, com faixas nos azimutes das
          gargantas: seis maximos sobre seis faixas no R1, nada no liso.

O instante da linha de cima e fixado pela GEOMETRIA, nao pelo resultado: o a_6 do R1 pica
quando a frente atravessa a primeira coroa (R_med ~ 1.50 = 0.99 do raio da coroa), que e a
janela pre-registrada na §12.4. Em t=50 o a_6 ja decaiu a 0.013 e a figura nao mostraria nada.

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
    p.add_argument("--rmed", type=float, default=1.50)
    p.add_argument("--n", type=int, default=600)
    p.add_argument("--iso", type=float, default=2.0)
    p.add_argument("--out", default="plots/fig_obstaculo.png")
    a = p.parse_args()

    SL, SR = serie(a.liso), serie(a.rug)
    ql = min(SL, key=lambda q: abs(q["rmed"] - a.rmed))
    qr = min(SR, key=lambda q: abs(q["rmed"] - a.rmed))
    pil = FT.carrega_pilares(qr["f"])
    L = 2.6
    al, fl = a6(ql["R"]); ar, fr = a6(qr["R"])
    print(f"liso: t={ql['t']:.1f} R_med={ql['rmed']:.3f} a6={al:.3f} fase={fl:.1f}")
    print(f"R1  : t={qr['t']:.1f} R_med={qr['rmed']:.3f} a6={ar:.3f} fase={fr:.1f}")

    fig, ax = plt.subplots(2, 2, figsize=(9.6, 9.0))
    painel(ax[0, 0], ql, pil, L, a.n, a.iso,
           f"(a) sem obstaculo — t = {ql['t']:.0f} s\n"
           f"$a_6$ = {al:.3f}   fase = {fl:.0f}°", tenue=True)
    painel(ax[0, 1], qr, pil, L, a.n, a.iso,
           f"(b) com obstaculos (R1, $\\Lambda$ = 28 dx) — t = {qr['t']:.0f} s\n"
           f"$a_6$ = {ar:.3f}   fase = {fr:.0f}°", tenue=False)
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
    th = np.arange(NTH) * 360.0 / NTH
    for k in range(6):
        c.axvspan(30 + 60 * k - 7, 30 + 60 * k + 7, color="#9a9855", alpha=0.16, lw=0)
    for q, cor, lbl, am, fa in ((ql, COR_LISO, "sem obstaculo", al, fl),
                                (qr, COR_RUG, "com obstaculos (R1)", ar, fr)):
        rel = q["R"] / q["R"].mean()
        c.plot(th, rel, color=cor, lw=0.7, alpha=0.28)
        # componente m=6 isolada: e exatamente o que a metrica mede (amplitude + fase)
        c.plot(th, 1.0 + am * np.cos(np.radians(6 * (th - fa))), color=cor, lw=2.4,
               label=f"{lbl} — $a_6$={am:.3f}, fase={fa:.0f}°")
    c.axhline(1.0, color="#999", lw=0.7, ls=":")
    c.set_xlim(0, 360); c.set_xticks(np.arange(0, 361, 60))
    c.set_xlabel(r"$\theta$ (graus)"); c.set_ylabel(r"$R(\theta)\,/\,\bar{R}$")
    c.set_title("(d) contorno desenrolado no instante do pico\n"
                "claro = $R(\\theta)$ cru;  grosso = componente $m$=6 medida", fontsize=10, pad=6)
    c.legend(fontsize=7.5, frameon=False, loc="lower left", ncol=1)
    c.grid(alpha=0.25, lw=0.5)

    fig.suptitle("A mesma expansao, com organizacao diferente", fontsize=12.5, y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    fig.savefig(a.out, dpi=190)
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
