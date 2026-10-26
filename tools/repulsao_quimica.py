# pass-l-aprovado: Pass L desbloqueado em 2026-09-17 (CLAUDE.md §12). Ferramenta de
# pos-processamento: le pilares.json e os HDF5, nao cria nem altera geometria de contorno.
"""Repulsao quimica: o acumulo de c_s a montante do pilar empurra o lider PARA LONGE dele?

Hipotese da CRITERIOS_RUGOSIDADE.md secao 18.4. O pilar enriquece c_s no lado voltado para o
inoculo (secao 18.2, fator 1.66 do dipolo isolante). Entao, entre um lider que chega e o pilar,
c_s cresce em direcao ao pilar e -grad(c_s) aponta para longe dele: as forcas que seguem
-grad(c_s) (Marangoni e flagelo) freariam o lider antes do contato mecanico.

    python tools/repulsao_quimica.py runs/rugosidade/R1_lambda28 runs/swarm/P2R23_conduz

DEFINICOES, fixadas antes de rodar:
  lider        (~filler) & (rho_b >= 0.1) & (r >= 0.7*R99) — a mesma de desfecho_contato.py
  pilar alvo   o pilar mais proximo do lider, so se estiver A FRENTE: cos(c - x, x) > 0
  d            distancia do lider a SUPERFICIE desse pilar, em dx
  e            vetor unitario do lider para o centro do pilar
  medida       cos_g = (-grad c_s) . e / |grad c_s|  — direcao do flagelo em relacao ao pilar
               a_p   = a_mar . e                     — Marangoni projetada no pilar
               O HDF5 so grava magnitudes (au_mar, au_flag); as direcoes sao RECONSTRUIDAS com as
               formulas do solver e validadas contra au_mar (correlacao impressa).
  PAREAMENTO   o mesmo calculo no LISO, para os lideres do liso, contra as MESMAS posicoes de
               pilar (que no liso nao existem). A diferenca R1 - liso por faixa de d e o efeito.
  faixas       d em [0,1), [1,2), [2,3), [3,5), [5,7), [7,10), [10,15) dx
  hipotese     CONFIRMADA se cos_g do R1 cai abaixo do liso nas faixas proximas (d < 3) e a
               diferenca some nas distantes (d >= 7). A forma por distancia e a prova (licao #68);
               uma diferenca constante em todas as faixas seria composicao, nao pilar.
  quadros      t >= 8 (antes disso nao ha lider perto de pilar)
"""

import argparse
import glob
import json
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

DX = 14.0 / 260.0
FAIXAS = [(0, 1), (1, 2), (2, 3), (3, 5), (5, 7), (7, 10), (10, 15)]
CHAVES = ("x", "y", "m", "rho", "h", "cs", "rho_b_grown", "is_filler", "au_mar")
BETA, GRAD_LO, GRAD_HI, DET_MIN = 5.0, 0.15, 0.5, 0.25


def dw(r, h):
    """dW/dr do CubicSpline 2D do PySPH."""
    s = 10.0 / (7.0 * np.pi * h * h)
    q = r / h
    d = np.where(q < 1.0, -3.0 * q + 2.25 * q * q, np.where(q < 2.0, -0.75 * (2.0 - q) ** 2, 0.0))
    return s * d / h


def gradientes(d, sel, T):
    """-grad(c_s) dos lideres com as formulas do solver: fontes = fluido (pilar fora);
    c_s sem filler como fonte (MarangoniForce/FlagellarForce); KGC com todo o fluido;
    gate da Marangoni por |grad rho_b| com todo o fluido. Devolve direcao do flagelo
    (sem KGC), a_mar reconstruida (com KGC e gate) e o |a_mar| para validacao."""
    x, y = d["x"], d["y"]
    H = float(np.median(d["h"]))
    viz = T.query_ball_point(np.c_[x[sel], y[sel]], 2.0 * H)
    ii = np.concatenate([np.full(len(v), k) for k, v in enumerate(viz)])
    jj = np.concatenate(viz).astype(int)
    di = sel[ii]
    k = jj != di
    ii, jj, di = ii[k], jj[k], di[k]
    xij, yij = x[di] - x[jj], y[di] - y[jj]
    r = np.hypot(xij, yij)
    g = dw(r, 0.5 * (d["h"][di] + d["h"][jj])) / np.maximum(r, 1e-12)
    DX, DY = g * xij, g * yij
    vol = d["m"][jj] / d["rho"][jj]
    n = len(sel)
    acc = lambda w: np.bincount(ii, weights=w, minlength=n)  # noqa: E731
    fonte = d["is_filler"][jj] < 0.5
    dcs = (d["cs"][jj] - d["cs"][di]) * fonte
    gx, gy = acc(vol * dcs * DX), acc(vol * dcs * DY)
    drb = d["rho_b_grown"][jj] - d["rho_b_grown"][di]
    grb = np.hypot(acc(vol * drb * DX), acc(vol * drb * DY))
    Mxx, Mxy = acc(-vol * xij * DX), acc(-vol * xij * DY)
    Myx, Myy = acc(-vol * yij * DX), acc(-vol * yij * DY)
    det = Mxx * Myy - Mxy * Myx
    ok = det >= DET_MIN
    sd = np.where(ok, det, 1.0)
    Lxx, Lxy = np.where(ok, Myy / sd, 1.0), np.where(ok, -Mxy / sd, 0.0)
    Lyx, Lyy = np.where(ok, -Myx / sd, 0.0), np.where(ok, Mxx / sd, 1.0)
    t = np.clip((grb - GRAD_LO) / (GRAD_HI - GRAD_LO), 0.0, 1.0)
    gate = np.where(grb >= GRAD_LO, np.where(t >= 1.0, 1.0, t * t * (3.0 - 2.0 * t)), 0.0)
    ax = -BETA * gate * (Lxx * gx + Lxy * gy)
    ay = -BETA * gate * (Lyx * gx + Lyy * gy)
    return np.c_[-gx, -gy], np.c_[ax, ay]


def coleta(run, C, rp):
    tree = cKDTree(C)
    linhas, val = [], []
    for fn in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        with h5py.File(fn, "r") as f:
            a = f["particles"]["fluid"]["arrays"]
            d = {k: np.asarray(a[k], dtype=float) for k in CHAVES}
            t = float(f["solver_data"].attrs["t"])
        if t < 8.0:
            continue
        r = np.hypot(d["x"], d["y"])
        fil = d["is_filler"] > 0.5
        r99 = np.percentile(r[(d["rho_b_grown"] >= 0.1) | fil], 99)
        L = np.where((~fil) & (d["rho_b_grown"] >= 0.1) & (r >= 0.7 * r99))[0]
        X = np.column_stack([d["x"][L], d["y"][L]])
        dist, k = tree.query(X)
        v = C[k] - X
        e = v / np.maximum(dist, 1e-12)[:, None]
        frente = np.einsum("ij,ij->i", v, X) > 0
        mg, am = gradientes(d, L, cKDTree(np.c_[d["x"], d["y"]]))
        val.append(np.c_[np.hypot(am[:, 0], am[:, 1]), d["au_mar"][L]])
        cos_g = np.einsum("ij,ij->i", mg, e) / np.maximum(np.hypot(mg[:, 0], mg[:, 1]), 1e-30)
        a_p = np.einsum("ij,ij->i", am, e)
        ds = (dist - rp) / DX
        for i in np.where(frente)[0]:
            linhas.append((t, ds[i], cos_g[i], a_p[i]))
    v = np.vstack(val)
    m = v[:, 1] > 0
    print(f"  validacao {run}: corr(|a_mar| reconstruida, au_mar gravada) = "
          f"{np.corrcoef(v[m, 0], v[m, 1])[0, 1]:.4f}, mediana da razao {np.median(v[m, 0] / v[m, 1]):.4f}")
    return np.array(linhas)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("rugoso")
    p.add_argument("liso")
    a = p.parse_args()

    meta = json.load(open(os.path.join(a.rugoso, "main_output", "pilares.json")))
    C, rp = np.asarray(meta["centros"]), float(meta["raio"])
    R = coleta(a.rugoso, C, rp)
    S = coleta(a.liso, C, rp)

    print(f"observacoes lider-quadro com pilar a frente: R1 {len(R)}, liso {len(S)}\n")
    print(f"{'d (dx)':>9} | {'n R1':>5} {'cos_g R1':>9} {'n liso':>7} {'cos_g liso':>10} "
          f"{'dif':>7} | {'a_p R1':>8} {'a_p liso':>9}")
    for lo, hi in FAIXAS:
        mr = (R[:, 1] >= lo) & (R[:, 1] < hi)
        ml = (S[:, 1] >= lo) & (S[:, 1] < hi)
        if mr.sum() == 0 or ml.sum() == 0:
            print(f"{lo:>4}-{hi:<4} | {int(mr.sum()):5d} {'-':>9} {int(ml.sum()):7d}")
            continue
        cr, cl = np.nanmedian(R[mr, 2]), np.nanmedian(S[ml, 2])
        er = 1.2533 * np.nanstd(R[mr, 2]) / np.sqrt(mr.sum())
        print(f"{lo:>4}-{hi:<4} | {int(mr.sum()):5d} {cr:+9.3f} {int(ml.sum()):7d} {cl:+10.3f} "
              f"{cr - cl:+7.3f} | {np.median(R[mr, 3]):+8.3f} {np.median(S[ml, 3]):+9.3f}"
              f"   (erro p50 R1 ~{er:.3f})")
    print("\ncos_g = cosseno entre -grad(c_s) e a direcao lider->pilar: +1 empurra PARA o pilar,"
          " -1 empurra PARA LONGE.")


if __name__ == "__main__":
    main()
