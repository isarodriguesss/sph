"""Painel (a) da figura de tese para varios runs lado a lado, mesma escala e mesmas cores.

    python plots/cmp_colonia.py runs/P2_fillerdonor runs/P2S_wakeseg [--t 50] [--L 5]
                                [--nomes P2 P2S] [--out plots/cmp_colonia_P2_P2S]
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fig_tese as FT


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t", type=float, default=50.0)
    ap.add_argument("--L", type=float, default=5.0)
    ap.add_argument("--n", type=int, default=700)
    ap.add_argument("--fecha", type=float, default=0.0)
    ap.add_argument("--iso", type=float, default=2.0, help="ignora agar isolado (dx); 0 desliga")
    ap.add_argument("--nomes", nargs="+", default=None)
    ap.add_argument("--out", default="plots/cmp_colonia")
    a = ap.parse_args()
    nomes = a.nomes or [r.rstrip("/").split("/")[-1] for r in a.runs]

    fig, ax = plt.subplots(1, len(a.runs), figsize=(3.2 * len(a.runs), 3.3),
                           constrained_layout=True)
    ext = [-a.L, a.L, -a.L, a.L]
    for k, (run, nome) in enumerate(zip(a.runs, nomes)):
        t, f = min(FT.frames(run), key=lambda p: abs(p[0] - a.t))
        g, img, borda, _ = FT.campo(FT.carrega(f), a.L, a.n, False, a.fecha, agar_iso=a.iso)
        rmax, rmin = FT.raios(g, borda)
        ax[k].imshow(img, origin="lower", extent=ext, cmap=FT.CMAP_COL, vmin=0, vmax=1,
                     interpolation="bilinear")
        FT.eixos(ax[k], a.L, rot=(k == 0))
        ax[k].set_title(rf"({chr(97 + k)}) {nome}, $t={t:.0f}$", loc="left")
        ax[k].text(0.03, 0.03, rf"$R_{{max}}={rmax:.2f}$  $R_{{min}}={rmin:.2f}$",
                   transform=ax[k].transAxes, fontsize=7.5, color="#2b2b2b")
    for ext_ in ("png", "pdf"):
        fig.savefig(f"{a.out}.{ext_}", dpi=300, bbox_inches="tight", pad_inches=0.05)
    print(f"{a.out}.png / .pdf")


if __name__ == "__main__":
    main()
