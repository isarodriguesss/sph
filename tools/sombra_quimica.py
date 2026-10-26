# pass-l-aprovado: Pass L desbloqueado em 2026-09-17 (CLAUDE.md §12). Ferramenta de
# pos-processamento: le pilares.json e os HDF5, nao cria nem altera geometria de contorno.
"""Sombra quimica em volta dos pilares: o fluxo zero de c_s e c_n deixa marca mensuravel?

Verificacao prevista em docs/PLANO_RUGOSIDADE.md secao 1 e nunca feita. Responde a pergunta 2
da serie (bloqueio quimico). Pos-processamento puro dos HDF5; nao toca o solver.

    python tools/sombra_quimica.py runs/rugosidade/R1_lambda28 runs/swarm/P2R23_conduz

DEFINICOES, fixadas antes de olhar o resultado:
  anel         fluido com distancia a SUPERFICIE do pilar em (0, W], W = 2 dx (plano secao 1)
  montante     anel com cos(x - c, c) <= -0.5 (lado voltado para o inoculo)
  jusante      anel com cos(x - c, c) >= +0.5 (lado voltado para fora)
  assimetria   S = ln( media(montante) / media(jusante) ), por pilar e quadro
  PAREAMENTO   S e calculado no rugoso E no liso, nas MESMAS posicoes de pilar, no quadro do
               liso de t mais proximo. O liso tem fluido dentro do disco; o anel exclui o disco
               nos dois. O efeito do pilar e dS = S_rugoso - S_liso: o gradiente radial natural
               (o lado de dentro sempre tem mais c_s) cancela.
  nivel        N = ln( media(anel rugoso) / media(anel liso) ): o pilar mudou o NIVEL local?
  classes      'dentro'  centro do pilar em r < 0.9 * R99 do rugoso
               'frente'  centro em R99 <= r < R99 + 1.03 (L_D do c_s no agar, plano secao 1) —
                         e onde o bloqueio quimico ANTES do toque apareceria (o halo)
  CONTROLE     'falsos' (-F): as mesmas posicoes giradas 30 graus em torno da origem, mantidas
               so as cujo anel fica a > 2 dx da superficie de qualquer pilar real. Mesmo calculo
               (disco excluido nos dois runs). Se o R1 so tiver um campo mais ingreme que o
               liso, dS dos falsos sai igual ao dos reais; o efeito do pilar e real - falso.
  janela       t in [35, 50] (CRITERIOS secao 4); a trajetoria inteira e impressa junto
  previsao     plano secao 1: diferenca montante/jusante de c_s < 30%, de c_n < 10%,
               ou seja |exp(dS) - 1| abaixo disso. Acima = canal quimico mensuravel.
"""

import argparse
import glob
import json
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

DX = 14.0 / 260.0
W = 2.0
L_D = 1.03
CAMPOS = ("cs", "c_n")


def le(fn):
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        d = {k: a[k][:] for k in ("x", "y", "rho_b_grown", "is_filler") + CAMPOS}
        d["t"] = float(f["solver_data"].attrs["t"])
    return d


def r99(d):
    r = np.hypot(d["x"], d["y"])
    return float(np.percentile(r[(d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)], 99))


def lados(d, C, rp):
    """Para cada pilar: indices de montante, jusante e anel inteiro."""
    P = np.column_stack([d["x"], d["y"]])
    tree = cKDTree(P)
    out = []
    for c in C:
        idx = np.asarray(tree.query_ball_point(c, rp + W * DX), dtype=int)
        v = P[idx] - c
        dist = np.hypot(v[:, 0], v[:, 1])
        ok = dist > rp
        idx, v, dist = idx[ok], v[ok], dist[ok]
        rc = np.hypot(*c)
        cosang = (v @ (c / rc)) / dist if rc > 0 else np.zeros(len(idx))
        out.append((idx[cosang <= -0.5], idx[cosang >= 0.5], idx))
    return out


def medias(d, sel, campo):
    m, j, a = sel
    if len(m) < 3 or len(j) < 3:
        return None
    f = d[campo]
    return f[m].mean(), f[j].mean(), f[a].mean()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("rugoso")
    p.add_argument("liso")
    a = p.parse_args()

    meta = json.load(open(os.path.join(a.rugoso, "main_output", "pilares.json")))
    C, rp = np.asarray(meta["centros"]), float(meta["raio"])
    ang = np.radians(30.0)
    F = C @ np.array([[np.cos(ang), np.sin(ang)], [-np.sin(ang), np.cos(ang)]])
    F = F[cKDTree(C).query(F)[0] - 2.0 * rp > (W + 2.0) * DX]
    rc = np.hypot(C[:, 0], C[:, 1])
    rf = np.hypot(F[:, 0], F[:, 1])
    fr = sorted(glob.glob(os.path.join(a.rugoso, "main_output", "main_*.hdf5")))
    fl = sorted(glob.glob(os.path.join(a.liso, "main_output", "main_*.hdf5")))
    tl = np.array([le(f)["t"] for f in fl])

    linhas = []
    SL = {}
    print(f"{len(C)} pilares, {len(F)} falsos, R_p = {rp / DX:.1f} dx, anel {W:.0f} dx\n")
    print(f"{'t_r':>6} {'t_l':>6} {'R99':>5} {'classe':>7} {'n':>3} | "
          + " | ".join(f"{c:>4}: dS p50  exp-1   N p50" for c in CAMPOS))
    for f in fr[1:]:
        dr = le(f)
        dl = le(fl[int(np.argmin(np.abs(tl - dr["t"])))])
        R = r99(dr)
        cls = {"dentro": rc < 0.9 * R, "frente": (rc >= R) & (rc < R + L_D)}
        cls["dentro-F"] = rf < 0.9 * R
        cls["frente-F"] = (rf >= R) & (rf < R + L_D)
        sel = {"": (lados(dr, C, rp), lados(dl, C, rp)), "-F": (lados(dr, F, rp), lados(dl, F, rp))}
        for nome, mask in cls.items():
            ar, al = sel["-F" if nome.endswith("-F") else ""]
            res = {}
            for campo in CAMPOS:
                dS, N = [], []
                for k in np.where(mask)[0]:
                    mr, ml = medias(dr, ar[k], campo), medias(dl, al[k], campo)
                    if mr is None or ml is None or min(mr + ml) <= 0:
                        continue
                    dS.append(np.log(mr[0] / mr[1]) - np.log(ml[0] / ml[1]))
                    SL.setdefault((nome, campo), []).append(np.log(ml[0] / ml[1]))
                    N.append(np.log(mr[2] / ml[2]))
                res[campo] = (np.array(dS), np.array(N))
            n = len(res["cs"][0])
            if n == 0:
                continue
            linhas.append((dr["t"], nome, res))
            txt = " | ".join(
                f"{np.median(res[c][0]):+8.3f} {np.expm1(np.median(res[c][0])):+6.1%} "
                f"{np.median(res[c][1]):+6.3f}"
                for c in CAMPOS
            )
            print(f"{dr['t']:6.1f} {dl['t']:6.1f} {R:5.2f} {nome:>8} {n:3d} | {txt}")

    print("\nS do LISO (escala da assimetria natural), mediana sobre a trajetoria:")
    for (nome, campo), v in sorted(SL.items()):
        print(f"  {nome:>8} {campo:>4}: S_liso = {np.median(v):+.3f}")
    print("\nJANELA t in [35, 50] — media +- desvio das medianas por quadro")
    for nome in ("dentro", "dentro-F", "frente", "frente-F"):
        for campo in CAMPOS:
            ds = [np.median(r[campo][0]) for t, c, r in linhas if c == nome and 35 <= t <= 50]
            nn = [np.median(r[campo][1]) for t, c, r in linhas if c == nome and 35 <= t <= 50]
            if not ds:
                continue
            ds, nn = np.array(ds), np.array(nn)
            print(f"  {nome:>7} {campo:>4}: dS = {ds.mean():+.3f} +- {ds.std():.3f} "
                  f"(montante/jusante x{np.exp(ds.mean()):.3f}), "
                  f"N = {nn.mean():+.3f} +- {nn.std():.3f} (nivel x{np.exp(nn.mean()):.3f}), "
                  f"{len(ds)} quadros")


if __name__ == "__main__":
    main()
