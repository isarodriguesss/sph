"""De onde vem o material que enche as baias depois de t=50 (o "disco tardio" do P2).

    python tools/diag_disco.py runs/P2_t100 runs/E11_t100 [--t0 50] [--t1 80]

Regiao F = baia em t0 (fora da borda renderizada de t0) que esta DENTRO do disco de raio
R_min(t1). Cada particula de colonia (rho_b >= 0.1 ou filler) em F no instante t1 e
classificada pelo estado em t0 — indice = identidade, porque nenhuma particula e removida e as
novas sao anexadas ao fim do array (premissa verificada no campo distante):
  nova-wake / nova-insert   criada depois de t0
  convertida-agar/limbo     existia em t0 fora da colonia; cruzou rho_b = 0.1
  colonia-que-entrou        ja era colonia em t0, fora de F; foi advectada para dentro
  colonia-ja-em-F           ja era colonia em t0 dentro de F (so borda: F e fora da colonia)
Para as convertidas, o peso de doador `sum V rho_b W` no instante medio, separado em filler e
vivas, diz se foi o filler doador (COL_FILLER_DONOR, licao #88) que as recrutou.
"""

import argparse
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, "plots")
import fig_tese as FT  # noqa: E402

DX = 14.0 / 260
L, N = 6.5, 450


def le(f):
    with h5py.File(f, "r") as h:
        a = h["particles"]["fluid"]["arrays"]
        d = {k: np.asarray(a[k], float) for k in
             ("x", "y", "m", "rho", "rho_b_grown", "is_filler", "is_wake", "cs")}
        d["t"] = float(h["solver_data"].attrs["t"])
    return d


def perto(fr, t):
    return min(fr, key=lambda p: abs(p[0] - t))[1]


def dentro(g, mask, x, y):
    cel = g[1] - g[0]
    i = np.clip(((y - g[0]) / cel).round().astype(int), 0, len(g) - 1)
    j = np.clip(((x - g[0]) / cel).round().astype(int), 0, len(g) - 1)
    return mask[i, j]


def doadores(d, idx):
    """Peso de doador sum V rho_b W sobre vizinhos filler e vivos (kernel cubico, h=1.8 dx)."""
    h = 1.8 * DX
    T = cKDTree(np.c_[d["x"], d["y"]])
    viz = T.query_ball_point(np.c_[d["x"][idx], d["y"][idx]], 2 * h)
    wf, wv = np.zeros(len(idx)), np.zeros(len(idx))
    V = d["m"] / d["rho"]
    for k, (i, v) in enumerate(zip(idx, viz)):
        v = np.asarray(v)
        v = v[v != i]
        r = np.hypot(d["x"][v] - d["x"][i], d["y"][v] - d["y"][i])
        w = FT.spline(r / h) * V[v] * d["rho_b_grown"][v]
        fil = d["is_filler"][v] > 0.5
        wf[k], wv[k] = w[fil].sum(), w[~fil].sum()
    return wf, wv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--t0", type=float, default=50.0)
    ap.add_argument("--t1", type=float, default=80.0)
    a = ap.parse_args()

    for run in a.runs:
        fr = FT.frames(run)
        d0, d1 = le(perto(fr, a.t0)), le(perto(fr, a.t1))
        dm = le(perto(fr, 0.5 * (a.t0 + a.t1)))
        n0 = len(d0["x"])

        # premissa: indice = identidade (campo distante nao se move)
        far = np.hypot(d0["x"], d0["y"]) > 6.0
        desloc = np.hypot(d1["x"][:n0][far] - d0["x"][far], d1["y"][:n0][far] - d0["y"][far])

        g, _, b0, _ = FT.campo(FT.carrega(perto(fr, a.t0)), L, N, False, 3.5)
        _, _, b1, _ = FT.campo(FT.carrega(perto(fr, a.t1)), L, N, False, 3.5)
        R0max, R0min = FT.raios(g, b0)
        R1max, R1min = FT.raios(g, b1)
        GX, GY = np.meshgrid(g, g)
        F = (np.hypot(GX, GY) < R1min) & ~b0
        cel = g[1] - g[0]

        col1 = (d1["rho_b_grown"] >= 0.1) | (d1["is_filler"] > 0.5)
        emF1 = dentro(g, F, d1["x"], d1["y"])
        alvo = np.where(col1 & emF1)[0]

        velho = alvo[alvo < n0]
        nova = alvo[alvo >= n0]
        col0 = (d0["rho_b_grown"] >= 0.1) | (d0["is_filler"] > 0.5)
        emF0 = dentro(g, F, d0["x"][velho], d0["y"][velho])
        conv = velho[~col0[velho]]
        agar0 = d0["rho_b_grown"][conv] <= 1e-12
        entrou = velho[col0[velho] & ~emF0]
        jaF = velho[col0[velho] & emF0]
        nw = int(np.sum(d1["is_wake"][nova] > 0.5))

        print(f"\n=== {run}   t0={d0['t']:.1f} -> t1={d1['t']:.1f}")
        print(f"premissa indice=identidade: desloc. max no campo distante {desloc.max():.1e} "
              f"(dx={DX:.3f})")
        print(f"R_min {R0min:.2f} -> {R1min:.2f}   R_max {R0max:.2f} -> {R1max:.2f}   "
              f"area de baia engolida F = {F.sum() * cel * cel:.2f} ({F.sum() * cel * cel / DX**2:.0f} dx^2)")
        tot = len(alvo)
        linhas = [
            ("nova — wake", nw),
            ("nova — insert", len(nova) - nw),
            ("convertida (era agar em t0)", int(agar0.sum())),
            ("convertida (era limbo em t0)", int((~agar0).sum())),
            ("colonia que entrou em F", len(entrou)),
            ("colonia ja em F em t0", len(jaF)),
        ]
        print(f"colonia em F no t1: {tot} particulas")
        for nome, n in linhas:
            print(f"  {nome:30s} {n:6d}  {100 * n / max(tot, 1):5.1f}%")

        if len(conv):
            # as convertidas estavam em F ja em t0, ou vieram de fora?
            conv_em_F0 = dentro(g, F, d0["x"][conv], d0["y"][conv])
            print(f"  convertidas que ja estavam em F em t0: {100 * conv_em_F0.mean():.0f}%")
            # quem doou: filler ou vivas, no instante medio (so as que ainda eram sub-quorum la)
            sub = conv[(dm["rho_b_grown"][conv] < 0.1) & (dm["is_filler"][conv] < 0.5)]
            if len(sub) > 20:
                s = np.random.default_rng(0).choice(sub, min(1500, len(sub)), replace=False)
                wf, wv = doadores(dm, s)
                frac = wf / np.maximum(wf + wv, 1e-30)
                print(f"  doador no instante medio (t={dm['t']:.0f}, n={len(s)}): fracao do peso "
                      f"vinda de FILLER p50 {np.median(frac):.2f}, media {frac.mean():.2f}; "
                      f"sem doador nenhum {100 * np.mean(wf + wv < 1e-12):.0f}%")


if __name__ == "__main__":
    main()
