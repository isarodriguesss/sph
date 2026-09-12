"""Ranking de rodadas com as mesmas reguas da figura de tese e do C5.

    python tools/rank_runs.py [--raio 2.4] [--painel plots/rank_t50.png] RUN...

Por rodada, na janela t in [35, min(50, t_fim)] (C5 e criterio de trajetoria; um instante
varia +-7 pontos, §2.2):
  soltos    pedacos da colonia RENDERIZADA (plots/fig_tese.py) fora do corpo principal
            com area > 20 dx^2 — o defeito "braco deslocado do nucleo"
  fora%     area da colonia fora do corpo principal
  C5a, C5b  alcance e fracao do componente do centro, ligacao 1.05 dx; C5b+L com limbo
No mesmo R99 (`--raio`; forma comparada em raio igual, nunca em t igual):
  amp, dedos, limpo   de tools/compara_osm.py
No ultimo frame ate t=50: Rmax/Rmin da borda renderizada (Rmin = fundo das baias), vivas.
Do log ate t=50: dR/dt em t[25,50], a_pressure mediana e % de picos > 4.
"""
import argparse
import csv
import os
import sys

import numpy as np
from scipy import ndimage as ndi

sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "plots"))
import compara_osm as CO  # noqa: E402
import fig_tese as FT  # noqa: E402

T_FIM = 50.5


def soltos(f, n=500, L=5.0):
    d = FT.carrega(f)
    g, img, borda, _ = FT.campo(d, L, n, False, 3.5)
    lab, _ = ndi.label(img > 0.5)
    cel = g[1] - g[0]
    c = len(g) // 2
    area = np.bincount(lab.ravel())[1:] * cel * cel
    nuc = lab[c, c]
    fora = area.sum() - (area[nuc - 1] if nuc > 0 else 0.0)
    big = int(np.sum((area > 20 * FT.DX ** 2) & (np.arange(1, len(area) + 1) != nuc)))
    return big, 100 * fora / max(area.sum(), 1e-12), g, borda


def log_stats(run):
    rows = [r for r in csv.DictReader(open(os.path.join(run, "log.csv")))
            if float(r["t"]) <= T_FIM]
    ap = np.array([float(r["a_pressure"]) for r in rows])
    return float(np.median(ap)), 100 * float(np.mean(ap > 4))


def avalia(run, raio):
    fr = [p for p in FT.frames(run) if p[0] <= T_FIM]
    if not fr or fr[-1][0] < 45:
        return None
    S = {s["t"]: s for s in CO.serie(run)}
    jan = [(t, f) for t, f in fr if t >= 35]
    so, fo, c5a, c5b, c5l = [], [], [], [], []
    for t, f in jan:
        b, pf, _, _ = soltos(f)
        so.append(b)
        fo.append(pf)
        s = S.get(t)
        if s:
            c5a.append(s["c5a"])
            c5b.append(s["c5b"])
            c5l.append(s["c5bl"])
    t_fim, f_fim = fr[-1]
    _, _, g, borda = soltos(f_fim)
    Rmax, Rmin = FT.raios(g, borda)
    s_fim = S[t_fim]
    ts = np.array(sorted(S))
    R = np.array([S[t]["R99"] for t in ts])
    w = (ts >= 25) & (ts <= T_FIM)
    vel = np.polyfit(ts[w], R[w], 1)[0] if w.sum() >= 3 else np.nan
    forma = None
    if R.max() >= raio:
        i = int(np.argmax(R >= raio))
        viz = [S[t] for t in ts[max(0, i - 1):i + 2]]
        forma = tuple(float(np.median([q[k] for q in viz])) for k in ("amp", "ded", "limpo"))
    apm, apk = log_stats(run)
    return dict(t=t_fim, so=np.mean(so), fo=np.mean(fo), c5a=np.mean(c5a), c5b=np.mean(c5b),
                c5l=np.mean(c5l), R99=s_fim["R99"], Rmax=Rmax, Rmin=Rmin, viva=s_fim["viva"],
                vel=vel, forma=forma, apm=apm, apk=apk)


def painel(runs, out, ncol=6):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nr = int(np.ceil(len(runs) / ncol))
    fig, ax = plt.subplots(nr, ncol, figsize=(2.4 * ncol, 2.55 * nr), constrained_layout=True)
    ax = np.atleast_2d(ax)
    for a in ax.ravel():
        a.axis("off")
    for a, run in zip(ax.ravel(), runs):
        fr = [p for p in FT.frames(run) if p[0] <= T_FIM]
        t, f = fr[-1]
        g, img, _, _ = FT.campo(FT.carrega(f), 5.0, 400, False, 3.5)
        a.imshow(img, origin="lower", extent=[-5, 5, -5, 5], cmap=FT.CMAP_COL, vmin=0, vmax=1)
        a.set_title(f"{os.path.basename(run.rstrip('/'))[:22]}  t={t:.0f}", fontsize=7)
    fig.savefig(out, dpi=130)
    print("->", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--raio", type=float, default=2.4)
    ap.add_argument("--painel", default=None)
    a = ap.parse_args()
    print(f"{'run':24s} {'t':>4s} {'soltos':>6s} {'fora%':>5s} {'C5a':>5s} {'C5b':>5s} "
          f"{'C5b+L':>5s} | R99={a.raio}: {'amp':>5s} {'dedos':>5s} {'limpo':>5s} | "
          f"{'R99':>4s} {'Rmax':>4s} {'Rmin':>4s} {'vivas':>5s} {'dR/dt':>6s} {'ap_med':>6s} "
          f"{'ap>4':>4s}")
    for run in a.runs:
        r = avalia(run, a.raio)
        nome = os.path.basename(run.rstrip("/"))[:24]
        if r is None:
            print(f"{nome:24s}  (nao chega a t=45)")
            continue
        fm = r["forma"]
        fs = f"{fm[0]:5.3f} {fm[1]:5.0f} {100 * fm[2]:4.0f}%" if fm else "  (nao alcanca)   "
        print(f"{nome:24s} {r['t']:4.0f} {r['so']:6.1f} {r['fo']:5.1f} {r['c5a']:5.2f} "
              f"{r['c5b']:5.1f} {r['c5l']:5.1f} |           {fs} | {r['R99']:4.2f} "
              f"{r['Rmax']:4.2f} {r['Rmin']:4.2f} {r['viva']:5d} {r['vel']:6.4f} "
              f"{r['apm']:6.2f} {r['apk']:3.0f}%")
    if a.painel:
        painel(a.runs, a.painel)


if __name__ == "__main__":
    main()
