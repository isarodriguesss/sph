"""Figura de expansao da colonia — 4 vistas complementares.

Uso:  python tools/plot_expansion.py runs/C4_t100 [--window 5 55] [--out FIG.png]

(a) R(t) em log-log com ajustes de lei de potencia e as inclinacoes de referencia
    da literatura sobrepostas ([T2] alpha=1 nutrient-rich, [T3] alpha=0.45).
(b) R(t) linear com o nutriente `c_n` no eixo secundario — mostra POR QUE o regime
    muda: a transicao de alpha coincide com c_n cruzando o piso do gate.
(c) Kymograph: raio da frente por setor angular ao longo do tempo. Cada faixa
    clara e UM dendrito; permite ver nascimento, competicao e estagnacao.
(d) Contorno da frente em 12 instantes, cor = tempo (colorbar). Mostra a FORMA
    da expansao: dedos que avancam vs baias que estagnam.
(e) Numero de dendritos ao longo do tempo, com as faixas de referencia.
(f) Velocidade da frente dR/dt — o patamar e o steady-state de [T2].
"""

import sys
import os
import re
import csv
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

NBINS = 180


def load(run):
    rows = list(csv.DictReader(open(os.path.join(run, "log.csv"))))
    it2t = {int(r["iteration"]): float(r["t"]) for r in rows}
    its = np.array(sorted(it2t))
    ts = np.array([it2t[i] for i in its])
    out = []
    for fn in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        it = int(re.search(r"main_(\d+)", fn).group(1))
        out.append((float(np.interp(it, its, ts)), fn))
    return out, rows


def front(fn, nbins=NBINS, thr=0.1):
    f = h5py.File(fn, "r")
    a = f["particles"]["fluid"]["arrays"]
    x, y = np.array(a["x"]), np.array(a["y"])
    rb = np.array(a["rho_b_grown"])
    m = rb > thr
    if m.sum() < 10:
        return np.zeros(nbins), 0.0
    th = np.arctan2(y[m], x[m])
    r = np.hypot(x[m], y[m])
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    idx = np.clip(np.digitize(th, edges) - 1, 0, nbins - 1)
    prof = np.zeros(nbins)
    for b in range(nbins):
        s = r[idx == b]
        if s.size:
            prof[b] = np.percentile(s, 95)
    # preenche setores vazios por interpolacao circular
    bad = prof == 0
    if bad.any() and (~bad).any():
        good = np.where(~bad)[0]
        prof[bad] = np.interp(np.where(bad)[0], good, prof[good], period=nbins)
    return prof, float(np.percentile(r, 99))


def fit(T, R, lo, hi):
    m = (T >= lo) & (T <= hi) & (T > 0)
    if m.sum() < 4:
        return None
    p = np.polyfit(np.log(T[m]), np.log(R[m]), 1)
    r2 = 1 - np.var(np.log(R[m]) - np.polyval(p, np.log(T[m]))) / np.var(np.log(R[m]))
    return p, r2, m


def main(argv):
    win = [5.0, 55.0]
    out = "runs/expansao.png"
    runs = []
    i = 0
    while i < len(argv):
        if argv[i] == "--window":
            win = [float(argv[i + 1]), float(argv[i + 2])]
            i += 3
        elif argv[i] == "--out":
            out = argv[i + 1]
            i += 2
        else:
            runs.append(argv[i])
            i += 1
    run = runs[0] if runs else "runs/C4_t100"

    series, rows = load(run)
    T, R, P = [], [], []
    for t, fn in series:
        if t <= 0:
            continue
        prof, rr = front(fn)
        if rr == 0:
            continue
        T.append(t)
        R.append(rr)
        P.append(prof)
    T, R, P = np.array(T), np.array(R), np.array(P)

    tl = np.array([float(r["t"]) for r in rows])
    cn = np.array([float(r["c_n_bio_arms"]) for r in rows])

    fig = plt.figure(figsize=(21, 12))
    gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.26)

    # ---- (a) log-log com referencias da literatura ----
    ax = fig.add_subplot(gs[0, 0])
    ax.loglog(T, R, "o", ms=5, color="#222", label="modelo (SPH)", zorder=3)
    f1 = fit(T, R, 0.6 * win[1], win[1])
    f2 = fit(T, R, win[1], T.max())
    if f1:
        p, r2, m = f1
        ax.loglog(T[m], np.exp(np.polyval(p, np.log(T[m]))), "-", lw=2.5,
                  color="#1b7837",
                  label=f"nutrient-rich: α={p[0]:.2f} (R²={r2:.3f})", zorder=4)
    if f2:
        p, r2, m = f2
        ax.loglog(T[m], np.exp(np.polyval(p, np.log(T[m]))), "-", lw=2.5,
                  color="#b2182b",
                  label=f"pós-esgotamento: α={p[0]:.2f} (R²={r2:.3f})", zorder=4)
    t0 = np.array([T.min(), T.max()])
    for al, c, lab in [(1.0, "#1b7837", "[T2] Srinivasan α=1 (swarming steady-state)"),
                       (0.45, "#b2182b", "[T3] Giverso α=0.45 (difusão-limitado)")]:
        ax.loglog(t0, R[0] * (t0 / T[0]) ** al, "--", lw=1.4, color=c, alpha=0.65,
                  label=lab)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("R da colônia (p99)")
    ax.set_title("(a) Lei de expansão vs literatura", fontweight="bold")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, which="both")

    # ---- (b) R(t) linear + nutriente ----
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(T, R, "o-", ms=4, color="#222", label="R(t)")
    ax.axvspan(win[0], win[1], color="#1b7837", alpha=0.08)
    ax.axvline(win[1], color="#b2182b", ls="--", lw=1.5)
    ax.text(win[1], R.max() * 0.35, " nutriente esgota\n (c_n < 0.4)",
            color="#b2182b", fontsize=9, va="center")
    ax.set_xlabel("t (s)")
    ax.set_ylabel("R da colônia")
    ax2 = ax.twinx()
    ax2.plot(tl, cn, color="#2166ac", lw=2, alpha=0.8)
    ax2.axhline(0.4, color="#2166ac", ls=":", lw=1.2)
    ax2.set_ylabel("c_n nos braços", color="#2166ac")
    ax2.tick_params(axis="y", colors="#2166ac")
    ax.set_title("(b) A transição de regime é o nutriente", fontweight="bold")
    ax.grid(alpha=0.25)

    # ---- (c) kymograph ----
    ax = fig.add_subplot(gs[0, 2])
    im = ax.pcolormesh(np.degrees(np.linspace(-180, 180, NBINS) * np.pi / 180),
                       T, P, shading="auto", cmap="magma")
    ax.axhline(win[1], color="w", ls="--", lw=1.2, alpha=0.8)
    ax.set_xlabel("ângulo (graus)")
    ax.set_ylabel("t (s)")
    ax.set_title("(c) Kymograph — cada faixa clara é um dendrito",
                 fontweight="bold")
    plt.colorbar(im, ax=ax, label="raio da frente")

    # ---- (d) contornos da frente, cor = tempo ----
    ax = fig.add_subplot(gs[1, 0])
    th = np.linspace(-np.pi, np.pi, NBINS)
    tt = np.append(th, th[0])
    norm = matplotlib.colors.Normalize(vmin=T.min(), vmax=T.max())
    cmap = plt.get_cmap("viridis")
    for k in np.linspace(0, len(T) - 1, 12).astype(int):
        pr = np.append(P[k], P[k][0])
        ax.plot(pr * np.cos(tt), pr * np.sin(tt), lw=1.6,
                color=cmap(norm(T[k])), alpha=0.9)
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("(d) Contorno da frente — cor = tempo", fontweight="bold")
    plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                 label="t (s)")

    # ---- (e) numero de dendritos ----
    nf = []
    for pr in P:
        f_ = pr - pr.mean()
        n_ = len(f_)
        nf.append(sum(1 for i in range(n_)
                      if f_[i] > 0 and f_[i] >= f_[i - 1] and f_[i] >= f_[(i + 1) % n_]))
    nf = np.array(nf)
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(T, nf, "o-", ms=4, lw=1.6, color="#222")
    ax.axhspan(15, 20, color="#1b7837", alpha=0.18)
    ax.text(T.max() * 0.45, 17.5, "PA14: 15–20 dendritos (reference.jpg)",
            fontsize=9, color="#1b7837")
    ax.axhspan(7, 9, color="#762a83", alpha=0.15)
    ax.text(T.max() * 0.45, 8, "[T1] Trinschek painel (b): 7–9", fontsize=9,
            color="#762a83")
    ax.axvline(win[1], color="#b2182b", ls="--", lw=1.2)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("nº de dendritos")
    ax.set_title("(e) Seleção competitiva de dendritos", fontweight="bold")
    ax.grid(alpha=0.25)

    # ---- (f) velocidade da frente ----
    ax = fig.add_subplot(gs[1, 2])
    v = np.gradient(R, T)
    ax.plot(T, v, "o-", ms=4, lw=1.6, color="#222", label="dR/dt (modelo)")
    if f1:
        ax.axhline(np.median(v[(T >= 0.6 * win[1]) & (T <= win[1])]), color="#1b7837",
                   ls="--", lw=1.6,
                   label="[T2] steady-state: V constante")
    ax.axvline(win[1], color="#b2182b", ls="--", lw=1.2)
    ax.text(win[1], v.max() * 0.85, " nutriente esgota", color="#b2182b", fontsize=9)
    ax.set_xlabel("t (s)")
    ax.set_ylabel("dR/dt")
    ax.set_title("(f) Velocidade da frente", fontweight="bold")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)

    fig.suptitle(
        f"Expansão da colônia — {os.path.basename(run)}   "
        f"(janela útil t={win[0]:.0f}–{win[1]:.0f}s)",
        fontsize=14, fontweight="bold")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    plt.savefig(out, dpi=110, bbox_inches="tight")
    print(f"figura salva em {out}")


if __name__ == "__main__":
    main(sys.argv[1:])
