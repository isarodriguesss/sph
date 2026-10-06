"""Caracteriza o MEIO antes de gastar rodada: regular contra desordenado, a phi fixo.

Sugestao do Prof. Cesar — descrever o obstaculo como GRAFO com regioes proibidas, com a
condicao de contorno descrita estocasticamente. O grafo e a triangulacao de Delaunay dos
centros (no = pilar, aresta = garganta, peso = distancia de superficie a superficie), que
e a abstracao padrao de meio poroso (pore network) e vale igual para rede e para campo
sorteado — e por isso resolve o problema que fazia a dose nao ser comparavel entre os dois.

`PILAR_JITTER` desloca cada centro por uma gaussiana de JITTER*Lambda, com rejeicao de
sobreposicao / dominio / zona de exclusao, entao a CONTAGEM e `phi` ficam IDENTICOS ao caso
regular: a unica coisa que muda e a distribuicao de gargantas (delta -> larga). E uma dose
de DESORDEM a phi fixo.

O que se mede aqui, tudo sem solver:
  1. distribuicao de gargantas (a dose no caso desordenado)
  2. percolacao do subgrafo de gargantas > largura do braco — na rede regular e uma CHAVE
     (todas iguais), so com desordem vira LIMIAR
  3. contatos por braco (varredura angular deterministica, 20 000 raios)
  4. ANISOTROPIA DO PROPRIO MEIO: espectro azimutal de d(theta) = distancia radial ate o
     primeiro pilar. E a causa geometrica do travamento de fase m=6 que o R1 mostrou
     (licao em CRITERIOS_RUGOSIDADE.md §12), e e ela que tem de DESAPARECER com desordem
     para o travamento ser atribuivel a rede.

    python tools/compara_meio.py --lam 28 --jitter 0 0.05 0.1 0.2 0.3
"""

import argparse

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import Delaunay, cKDTree

DX = 14.0 / 260.0
DOM = 7.0
R_EXCL = 1.2
BRACO = 5.4 * DX
NTH = 7200


def campo(lam_dx, a_dx, jitter, seed, rot=0.0):
    """Replica o gerador de main.py: rede, clip, exclusao, depois jitter a contagem fixa."""
    lam, rp = lam_dx * DX, 0.5 * a_dx * DX
    ang = np.radians(rot)
    n = int(np.ceil(2.0 * DOM / lam)) + 2
    i, j = np.meshgrid(np.arange(-n, n + 1), np.arange(-n, n + 1))
    cx = lam * (i + 0.5 * j).ravel()
    cy = lam * (np.sqrt(3.0) / 2.0) * j.ravel()
    cx, cy = (cx * np.cos(ang) - cy * np.sin(ang), cx * np.sin(ang) + cy * np.cos(ang))
    b = 2.0 * DX
    ok = ((cx > -DOM + rp + b) & (cx < DOM - rp - b)
          & (cy > -DOM + rp + b) & (cy < DOM - rp - b)
          & (np.hypot(cx, cy) >= R_EXCL))
    base = np.column_stack([cx[ok], cy[ok]])
    if jitter <= 0.0:
        return base, rp
    rng = np.random.default_rng(seed)
    amp = jitter * lam

    def valido(p):
        v = ((p[:, 0] > -DOM + rp + b) & (p[:, 0] < DOM - rp - b)
             & (p[:, 1] > -DOM + rp + b) & (p[:, 1] < DOM - rp - b)
             & (np.hypot(p[:, 0], p[:, 1]) >= R_EXCL))
        par = cKDTree(p).query_pairs(2.0 * rp, output_type="ndarray")
        if len(par):
            v[np.unique(par[:, 1])] = False
        return v

    pos = base + amp * rng.standard_normal(base.shape)
    for _ in range(200):
        ruim = ~valido(pos)
        if not ruim.any():
            break
        pos[ruim] = base[ruim] + amp * rng.standard_normal((int(ruim.sum()), 2))
    else:
        pos[~valido(pos)] = base[~valido(pos)]
    return pos, rp


def grafo(C, rp):
    ar = set()
    for s in Delaunay(C).simplices:
        for a, b in ((0, 1), (1, 2), (2, 0)):
            ar.add((min(s[a], s[b]), max(s[a], s[b])))
    e = np.array(sorted(ar))
    g = np.linalg.norm(C[e[:, 0]] - C[e[:, 1]], axis=1) - 2.0 * rp
    return e, g


def percola(C, e, g, lim):
    m = g > lim
    if m.sum() == 0:
        return 0, 0.0
    M = coo_matrix((np.ones(int(m.sum())), (e[m, 0], e[m, 1])), shape=(len(C), len(C)))
    nc, lab = connected_components(M, directed=False)
    return nc, np.bincount(lab).max() / len(C)


def raios(C, rp, rmax, captura):
    """d(theta) ate o primeiro pilar e contatos por raio, varredura deterministica."""
    th = np.arange(NTH) * 2.0 * np.pi / NTH
    u = np.column_stack([np.cos(th), np.sin(th)])
    proj = C @ u.T
    perp = np.abs(C[:, 0:1] * u[:, 1] - C[:, 1:2] * u[:, 0])
    hit = (perp < captura) & (proj > R_EXCL) & (proj < rmax)
    nc = hit.sum(0)
    ent = np.where(hit, proj - np.sqrt(np.maximum(captura**2 - perp**2, 0.0)), np.inf)
    d1 = ent.min(0)
    d1[~np.isfinite(d1)] = rmax
    return d1, nc


def modo(f, m):
    F = np.fft.rfft(f)
    return 2.0 * np.abs(F[m]) / np.abs(F[0]), np.degrees(-np.angle(F[m]) / m) % (360.0 / m)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--lam", type=float, default=28.0)
    p.add_argument("--a", type=float, default=10.0)
    p.add_argument("--jitter", type=float, nargs="+", default=[0.0, 0.05, 0.1, 0.2, 0.3])
    p.add_argument("--seed", type=int, default=20261006)
    p.add_argument("--rmax", type=float, default=3.852)
    a = p.parse_args()

    print(f"Lambda = {a.lam:.0f} dx, A = {a.a:.0f} dx, r_max = {a.rmax:.3f}, "
          f"braco = 5.4 dx, semente {a.seed}\n")
    print(f"{'jitter':>7} {'N':>4} {'phi':>6} | {'garganta p10/p50/p90 (dx)':>25} {'desvio':>7}"
          f" | {'<braco':>6} {'perc%':>6} | {'contatos':>8} {'>=1':>5}"
          f" | {'a6':>6} {'fase':>6} {'a12':>6}")
    for jt in a.jitter:
        C, rp = campo(a.lam, a.a, jt, a.seed)
        phi = len(C) * np.pi * rp**2 / (2 * DOM) ** 2
        e, g = grafo(C, rp)
        estreito = 100.0 * (g < BRACO).mean()
        nc, frac = percola(C, e, g, BRACO)
        d1, ncont = raios(C, rp, a.rmax, rp + BRACO / 2)
        a6, f6 = modo(d1, 6)
        a12, _ = modo(d1, 12)
        print(f"{jt:7.2f} {len(C):4d} {100*phi:5.1f}% | "
              f"{np.percentile(g,10)/DX:7.2f} {np.median(g)/DX:7.2f} {np.percentile(g,90)/DX:7.2f}"
              f" {np.std(g)/DX:7.2f} | {estreito:5.1f}% {100*frac:5.0f}% | "
              f"{ncont.mean():8.2f} {100*(ncont>=1).mean():4.0f}% | "
              f"{a6:6.3f} {f6:6.1f} {a12:6.3f}")
    print("\n`<braco` = fracao das gargantas mais estreitas que o braco (intransponiveis sem")
    print("deformar);  `perc%` = fracao dos nos no maior componente do subgrafo de gargantas")
    print("largas;  `a6/fase` = anisotropia do MEIO (azimute do primeiro pilar).")


if __name__ == "__main__":
    main()
