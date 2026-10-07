"""Morfometria das fotos de referencia, com recorte e limiar DECLARADOS por imagem.

Por que existe: a medida de nucleo depende da DEFINICAO, e as duas em uso discordam por
3-5x. O criterio de ocupacao (`maior raio que cobre >=98% das direcoes`) e fragil em foto
de placa — um entalhe escuro em uma direcao o zera, e no painel A da reference.jpg ele
varre de 0.000 a 0.653 conforme o limiar. O circulo inscrito (EDT no centroide) nao exige
cobrir todas as direcoes e e estavel. Aqui as duas sao medidas lado a lado, na mesma
mascara, para que a comparacao com os runs tenha base.

Cada fonte declara o seu recorte e o seu metodo de limiar — nao ha segmentacao automatica
que sirva para foto de placa com aro brilhante, anotacao em vermelho e iluminacao desigual.
`--mascaras` salva o painel de conferencia: **olhar as mascaras antes de citar qualquer
numero** (licoes #68, #101 — calibrar a regua antes de ranquear por ela).

    python tools/lit_fotos.py
    python tools/lit_fotos.py --mascaras plots/lit_fotos_mascaras.png
"""

import argparse
import os

import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.signal import find_peaks

NTH = 720

# nome, arquivo, (y0,y1,x0,x1), metodo, parametro
FONTES = [
    # nome, arquivo, (y0,y1,x0,x1), metodo, parametro. Limiar ajustado POR IMAGEM e
    # conferido na mascara (--mascaras): foto de placa tem aro brilhante, anotacao e
    # iluminacao desigual, e nao ha automatico que sirva para todas.
    ("PA14 ref.jpg A",   "assets/reference.jpg",                            (4, 152, 3, 158),      "vale", 0),
    ("PA14 ref.jpg C",   "assets/reference.jpg",                            (4, 152, 322, 476),    "vale", 0),
    ("PA14 BHI",         "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (15, 315, 10, 350),   "pct", 92),
    ("PA14 LB",          "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (15, 315, 360, 700),  "pct", 92),
    ("PA14 M8",          "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (355, 650, 35, 330),  "pct", 72),
    ("PA14 M9",          "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (355, 650, 385, 680), "pct", 62),
    ("PA14 BM2",         "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (685, 980, 35, 330),  "pct", 68),
    ("PA14 PGM",         "assets/Pseudomonas Aeruginosa Swarm Pattern.webp", (700, 980, 385, 680), "pct", 62),
    ("PA14 M9+CAA",      "assets/Motilidade de Pseudomonas Aeruginosa.jpg",  (120, 900, 140, 920), "pct", 62),
    ("PA14 placa 16",    "assets/Swarm Pseudomonas.jpeg",                    (55, 400, 50, 400),   "pct", 55),
    ("PA14 WT 0% 24h",   "assets/gr1 lrg 1660x2544.jpg",                     (130, 360, 1050, 1290), "pct", 72),
    ("PA14 WT .02% 24h", "assets/gr1 lrg 1660x2544.jpg",                     (690, 925, 1050, 1290), "pct", 68),
]



def vale(v):
    h, be = np.histogram(v, 64, (0, 256))
    h = ndi.gaussian_filter1d(h.astype(float), 1.5)
    c = (be[:-1] + be[1:]) / 2
    pico = int(np.argmax(h[:32]))
    mins = [i for i in np.where((h[1:-1] < h[:-2]) & (h[1:-1] < h[2:]))[0] + 1 if i > pico]
    return c[mins[0]] if mins else c[pico] + 40


def mascara(p, metodo, par):
    lim = vale(p) if metodo == "vale" else np.percentile(p, par)
    m = ndi.binary_opening(p > lim, np.ones((3, 3)))
    lab, n = ndi.label(m)
    if n == 0:
        return None
    ar = ndi.sum(m, lab, np.arange(1, n + 1))
    return ndi.binary_fill_holes(lab == 1 + np.argmax(ar))


def mede(mm):
    cy, cx = ndi.center_of_mass(mm)
    edt = ndi.distance_transform_edt(mm)
    r = np.arange(0.0, 0.995 * min(mm.shape) / 2, 0.5)
    th = np.linspace(-np.pi, np.pi, NTH, endpoint=False)
    X = cx + r[None, :] * np.cos(th[:, None])
    Y = cy + r[None, :] * np.sin(th[:, None])
    d = ndi.map_coordinates(mm.astype(float), [Y, X], order=0) > 0.5
    if not d.any():
        return None
    Rt = np.where(d.any(1), r[d.shape[1] - 1 - np.argmax(d[:, ::-1], 1)], 0.0)
    R1 = Rt.max()
    o = d.mean(0)
    i = np.argmax(o < 0.98)
    n_occ = (r[i] if i > 0 else 0.0) / R1
    n_ins = edt[int(round(cy)), int(round(cx))] / R1
    s = ndi.uniform_filter1d(np.r_[Rt, Rt, Rt], 5)[NTH:2 * NTH]
    pk, _ = find_peaks(np.r_[s, s], prominence=0.06 * s.max(), distance=NTH / 36)
    nb = int(np.sum((pk >= NTH // 2) & (pk < 3 * NTH // 2)))
    gy, gx = np.mgrid[:mm.shape[0], :mm.shape[1]]
    rr = np.hypot(gx - cx, gy - cy)
    braco = mm & (rr > 0.5 * R1) & (rr < 0.95 * R1)
    larg = 2 * np.median(edt[braco]) if braco.any() else np.nan
    return dict(n_occ=n_occ, n_ins=n_ins, R1=R1, n=nb,
                ocup=mm.sum() / (np.pi * R1 ** 2),
                larg_R=larg / R1, AR=(R1 - edt[int(round(cy)), int(round(cx))]) / larg
                if larg and larg == larg else np.nan,
                passo=2 * np.pi * 0.6 / nb if nb else np.nan, mask=mm)


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--mascaras", default=None)
    a = a.parse_args()
    linhas, masks = [], []
    for nome, arq, (y0, y1, x0, x1), met, par in FONTES:
        if not os.path.exists(arq):
            print(f"{nome:20} arquivo ausente"); continue
        im = np.asarray(Image.open(arq).convert("L"), float)
        p = im[y0:y1, x0:x1]
        mm = mascara(p, met, par)
        if mm is None:
            print(f"{nome:20} sem mascara"); continue
        r = mede(mm)
        if r is None:
            print(f"{nome:20} sem medida"); continue
        linhas.append((nome, r)); masks.append((nome, p, mm))
    print(f"\n{'fonte':>18} {'nuc occ':>8} {'nuc ins':>8} {'dedos':>6} {'AR':>6}"
          f" {'larg/R':>7} {'passo/R':>8} {'area/disco':>11}")
    for nome, r in linhas:
        print(f"{nome:>18} {r['n_occ']:8.3f} {r['n_ins']:8.3f} {r['n']:6d} {r['AR']:6.1f}"
              f" {r['larg_R']:7.3f} {r['passo']:8.3f} {r['ocup']:11.3f}")
    if a.mascaras and masks:
        from matplotlib import pyplot as plt
        k = len(masks); c = 4; l = (k + c - 1) // c
        fig, ax = plt.subplots(l, c, figsize=(3.1 * c, 3.2 * l))
        for e, (nome, p, mm) in zip(np.ravel(ax), masks):
            e.imshow(p, cmap="gray"); e.contour(mm, [0.5], colors="#d23a3a", linewidths=0.9)
            e.set_title(nome, fontsize=8); e.set_xticks([]); e.set_yticks([])
        for e in np.ravel(ax)[len(masks):]:
            e.axis("off")
        fig.suptitle("mascaras — conferir ANTES de citar os numeros", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.97]); fig.savefig(a.mascaras, dpi=150)
        print(f"\n-> {a.mascaras}")


if __name__ == "__main__":
    main()
