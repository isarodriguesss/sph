"""Desfecho do contato lider-pilar: PARADA / DESVIO / DIVISAO (CRITERIOS_RUGOSIDADE.md §4).

E a metrica MECANISTICA pre-registrada: sem ela um resultado quase-nulo nao e publicavel
(§7), porque "nao atrapalha" exige demonstrar que houve contato E o que ele fez.

Mecanismo que a motiva (licoes #98, #104): cada dedo e o rastro de UM lider; lider que
PARA deixa de ser lider (o limiar `r >= 0.7*r99` sobe com a colonia) e o dedo morre com
ele. Logo "perda de dedo" e extincao de lider, nao afinamento.

O CONTROLE e interno e e o ponto: lideres param tambem SEM pilar (licao #104 mediu o
P2R16 perdendo 33 de 53 lideres entre t=10 e 17). Entao a medida e
P(parada | contato) contra P(parada | livre) DENTRO da mesma rodada, nunca a taxa
absoluta.

Identidade = indice do array (verificado no R1: n cresce monotonicamente, insercoes vao
para o fim, deslocamento maximo de indice comum 5-7 dx = movimento real).

    python tools/desfecho_contato.py runs/rugosidade/R1_lambda28

DEFINICOES, fixadas antes de olhar o resultado:
  lider        (~filler) & (rho_b >= 0.1) & (r >= 0.7*r99), igual a `_alarga_rastro`
  d_surf       distancia ao centro de pilar mais proximo menos R_p
  em contato   d_surf < LIM dx. **A VARREDURA DO LIMIAR E OBRIGATORIA, nao escolha livre**
               (licao #68): o efeito e de curto alcance, entao a razao tem de CRESCER quando o
               limiar aperta e convergir a 1 quando ele abre. Se nao escalar com a distancia ao
               pilar, nao e contato — e composicao. Medido no R1: 10.8x a 1 dx, 1.1x a 10 dx,
               com o controle (livre) achatado em 7.7-10.0% em todos os limiares.
  PARADA       avanco radial no intervalo < 0.25 x mediana do avanco dos lideres
  DESVIO       deslocamento azimutal com SINAL, r*dtheta (serve tambem a quiralidade)
  DIVISAO      nivel de GRUPO: grupos de lideres ligados a 6 dx (RASTRO_PONTA_LINK)
"""

import argparse
import glob
import json
import os

import h5py
import numpy as np
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

DX = 14.0 / 260.0
R_MIN_FRAC = 0.7
LINK = 6.0
PARADA_FRAC = 0.25


def frames(run):
    for fn in sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5"))):
        with h5py.File(fn, "r") as f:
            a = f["particles"]["fluid"]["arrays"]
            d = {k: a[k][:] for k in ("x", "y", "rho_b_grown", "is_filler")}
            d["t"] = float(f["solver_data"].attrs["t"])
        yield d


def grupos(d, tree, rp):
    L, _, _ = lideres(d)
    if len(L) == 0:
        return L, np.array([], int), 0
    t = cKDTree(np.column_stack([d["x"][L], d["y"][L]]))
    ng, g = connected_components(t.sparse_distance_matrix(t, LINK * DX), directed=False)
    return L, g, ng


def lideres(d):
    rb, fil = d["rho_b_grown"], d["is_filler"] > 0.5
    r = np.hypot(d["x"], d["y"])
    r99 = float(np.percentile(r[(rb >= 0.1) | fil], 99))
    return np.where((~fil) & (rb >= 0.1) & (r >= R_MIN_FRAC * r99))[0], r, r99


def main():
    p = argparse.ArgumentParser()
    p.add_argument("run")
    p.add_argument("--lim", type=float, default=1.0, help="limiar de contato, em dx")
    a = p.parse_args()

    meta = json.load(open(os.path.join(a.run, "main_output", "pilares.json")))
    C = np.asarray(meta["centros"])
    rp = float(meta["raio"])
    tree = cKDTree(C)
    F = list(frames(a.run))

    reg = []   # (t, idx, contato, dr, arco, segue_lider, dsurf)
    print(f"pilares: {len(C)} centros, R_p = {rp / DX:.1f} dx; limiar de contato "
          f"{a.lim:.1f} dx\n")
    print(f"{'t':>6} {'lideres':>8} {'em contato':>11} {'grupos':>7} {'dr_med(dx)':>11}")
    for k in range(len(F) - 1):
        d0, d1 = F[k], F[k + 1]
        L, r0, _ = lideres(d0)
        if len(L) == 0:
            continue
        L1, r1, _ = lideres(d1)
        set1 = set(L1.tolist())
        dsurf = tree.query(np.column_stack([d0["x"][L], d0["y"][L]]))[0] - rp
        cont = dsurf < a.lim * DX
        # avanco radial e arco azimutal com sinal
        th0 = np.arctan2(d0["y"][L], d0["x"][L])
        th1 = np.arctan2(d1["y"][L], d1["x"][L])
        dth = (th1 - th0 + np.pi) % (2 * np.pi) - np.pi
        dr = r1[L] - r0[L]
        arco = r0[L] * dth
        segue = np.array([int(i) in set1 for i in L])
        # grupos de lideres
        tl = cKDTree(np.column_stack([d0["x"][L], d0["y"][L]]))
        ng, _ = connected_components(tl.sparse_distance_matrix(tl, LINK * DX),
                                     directed=False)
        print(f"{d0['t']:6.1f} {len(L):8d} {int(cont.sum()):5d} ({100*cont.mean():4.0f}%)"
              f" {ng:7d} {np.median(dr)/DX:11.2f}")
        for i in range(len(L)):
            reg.append((d0["t"], int(L[i]), bool(cont[i]), dr[i], arco[i],
                        bool(segue[i]), dsurf[i]))

    if not reg:
        raise SystemExit("sem lideres")
    T = np.array([x[0] for x in reg]); CT = np.array([x[2] for x in reg])
    ID = np.array([x[1] for x in reg]); DS = np.array([x[6] for x in reg])
    DR = np.array([x[3] for x in reg]); AR = np.array([x[4] for x in reg])
    SG = np.array([x[5] for x in reg])

    # limiar de parada por intervalo (auto-normalizado)
    par = np.zeros(len(reg), bool)
    for t in np.unique(T):
        m = T == t
        par[m] = DR[m] < PARADA_FRAC * np.median(DR[m])

    def taxa(mask, sub):
        n = int(mask.sum())
        if n == 0:
            return "   n=0"
        k = int((mask & sub).sum())
        p = k / n
        return f"{100*p:5.1f}% +-{100*np.sqrt(p*(1-p)/n):4.1f}  (n={n})"

    print(f"\n=== DESFECHO, {len(reg)} observacoes lider-quadro ===")
    print(f"{'':16} {'EM CONTATO':>26} {'LIVRE':>26}")
    for lbl, sub in (("PARADA", par), ("deixa de ser lider", ~SG)):
        print(f"{lbl:16} {taxa(CT, sub):>26} {taxa(~CT, sub):>26}")
    print(f"\n{'avanco radial dr (dx)':24} contato p50={np.median(DR[CT])/DX:6.2f}"
          f"   livre p50={np.median(DR[~CT])/DX:6.2f}")
    print(f"{'|arco| azimutal (dx)':24} contato p50={np.median(np.abs(AR[CT]))/DX:6.2f}"
          f"   livre p50={np.median(np.abs(AR[~CT]))/DX:6.2f}")
    print(f"{'arco COM SINAL (dx)':24} contato med={np.mean(AR[CT])/DX:+6.2f}"
          f"   livre med={np.mean(AR[~CT])/DX:+6.2f}   (quiralidade: != 0 ?)")
    s = AR[CT] / DX
    print(f"{'  erro do sinal':24} +-{np.std(s)/np.sqrt(len(s)):.2f} (contato)")

    print(f"\n=== VARREDURA DO LIMIAR (obrigatoria — licao #68) ===")
    print(f"{'lim(dx)':>8} {'n_cont':>7} {'distintas':>10} {'P(par|cont)':>14}"
          f" {'P(par|livre)':>14} {'razao':>7} {'dr cont/livre':>14}")
    for lim in (1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0):
        c = DS < lim * DX
        if c.sum() < 3:
            print(f"{lim:8.1f} {int(c.sum()):7d}   (poucos)")
            continue
        pc, pl = par[c].mean(), par[~c].mean()
        nc = int(c.sum())
        ec = 100 * np.sqrt(pc * (1 - pc) / nc)
        print(f"{lim:8.1f} {nc:7d} {len(np.unique(ID[c])):10d} {100*pc:8.1f}%+-{ec:4.1f}"
              f" {100*pl:13.1f}% {pc/pl if pl > 0 else np.nan:7.1f}x"
              f" {np.median(DR[c])/np.median(DR[~c]):14.2f}")

    print(f"\n=== DIVISAO / FUSAO / EXTINCAO de grupo (ligacao {LINK:.0f} dx) ===")
    tot = dict(div=0, fus=0, morre=0, nasce=0, div_pil=0, morre_pil=0)
    for k in range(len(F) - 1):
        L0, g0, n0 = grupos(F[k], tree, rp)
        L1, g1, n1 = grupos(F[k + 1], tree, rp)
        if not len(L0) or not len(L1):
            continue
        m1 = {int(i): gg for i, gg in zip(L1, g1)}
        ds0 = tree.query(np.column_stack([F[k]["x"][L0], F[k]["y"][L0]]))[0] - rp
        for gg in range(n0):
            mem = [int(i) for i in L0[g0 == gg]]
            dest = {m1[i] for i in mem if i in m1}
            perto = ds0[g0 == gg].min() / DX < 3.0
            if not dest:
                tot["morre"] += 1
                tot["morre_pil"] += perto
            elif len(dest) > 1:
                tot["div"] += 1
                tot["div_pil"] += perto
        org = {}
        for i, gg in zip(L0, g0):
            if int(i) in m1:
                org.setdefault(m1[int(i)], set()).add(gg)
        tot["fus"] += sum(1 for v in org.values() if len(v) > 1)
        tot["nasce"] += n1 - len(org)
    print(f"  divisoes={tot['div']} (com pilar a <3 dx: {tot['div_pil']})   fusoes={tot['fus']}")
    print(f"  grupos extintos={tot['morre']} (com pilar a <3 dx: {tot['morre_pil']})"
          f"   grupos novos={tot['nasce']}")


if __name__ == "__main__":
    main()
