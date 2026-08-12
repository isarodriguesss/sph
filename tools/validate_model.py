"""Bateria de validacao do modelo com ancoras na literatura.

Uso:  python tools/validate_model.py runs/C4_t100 [--window 5 55]

Cada teste cita a referencia e o valor esperado. Nao ha "passa/falha" cego:
onde o modelo diverge, o relatorio diz por quanto e em que regime.

REFERENCIAS (ver CLAUDE.md §3.0)
  [T1] Trinschek, John & Thiele 2018, Soft Matter 14, 4464 — halo de surfactante
       alem da biomassa; painel (b) Fingering com 7-9 dedos
  [T2] Srinivasan, Kaplan & Mahadevan 2019, eLife 8, e42697 — swarming em regime
       nutrient-rich e steady-state de velocidade CONSTANTE (R ~ t, alpha=1)
  [T3] Giverso, Verani & Ciarletta 2016, Biomech Model Mechanobiol 15, 643 —
       dedos difusao-limitados crescem com t^0.45 (~sqrt(t))
  [T6] Liu & Liu 2003 §3.3 — ~20 vizinhos p/ gradiente, ~35 p/ Laplaciano
  [T7] Violeau 2012 §3.4/3.6 — particao da unidade; sigma_a<0.85 ~ 15% de erro
  reference.jpg (Michiels et al., PA14) — 15-20 dendritos radiais, AR >= 1:5
"""

import sys
import os
import re
import csv
import glob
import numpy as np
import h5py
from scipy.spatial import cKDTree

DX = 14.0 / 260
H_FACTOR = 1.8


def load_series(run):
    """Retorna [(t, arrays)] ordenado, com t vindo do log.csv."""
    logp = os.path.join(run, "log.csv")
    rows = list(csv.DictReader(open(logp)))
    it2t = {int(r["iteration"]): float(r["t"]) for r in rows}
    its = np.array(sorted(it2t))
    ts = np.array([it2t[i] for i in its])
    out = []
    for fn in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        it = int(re.search(r"main_(\d+)", fn).group(1))
        out.append((float(np.interp(it, its, ts)), fn))
    return out, rows


def arrays(fn):
    f = h5py.File(fn, "r")
    a = f["particles"]["fluid"]["arrays"]
    return {k: np.array(a[k]) for k in a if np.array(a[k]).size}


def front_radius(x, y, rb, nbins=72, thr=0.1):
    """Raio da frente por setor angular — base do kymograph e da contagem de dedos."""
    m = rb > thr
    if m.sum() < 10:
        return np.zeros(nbins)
    th = np.arctan2(y[m], x[m])
    r = np.hypot(x[m], y[m])
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    idx = np.clip(np.digitize(th, edges) - 1, 0, nbins - 1)
    out = np.zeros(nbins)
    for b in range(nbins):
        s = r[idx == b]
        if s.size:
            out[b] = np.percentile(s, 95)
    return out


def count_fingers(fr):
    """Picos locais no perfil angular da frente (circular)."""
    f = fr - fr.mean()
    n = len(f)
    return int(
        sum(
            1
            for i in range(n)
            if f[i] > 0 and f[i] >= f[i - 1] and f[i] >= f[(i + 1) % n]
        )
    )


def expansion_law(series):
    T, R = [], []
    for t, fn in series:
        if t <= 0:
            continue
        d = arrays(fn)
        m = d["rho_b_grown"] > 0.1
        if m.sum() < 20:
            continue
        T.append(t)
        R.append(np.percentile(np.hypot(d["x"], d["y"])[m], 99))
    return np.array(T), np.array(R)


def fit_alpha(T, R, lo, hi):
    m = (T >= lo) & (T <= hi)
    if m.sum() < 4:
        return np.nan, np.nan, 0
    p = np.polyfit(np.log(T[m]), np.log(R[m]), 1)
    r2 = 1 - np.var(np.log(R[m]) - np.polyval(p, np.log(T[m]))) / np.var(np.log(R[m]))
    return p[0], r2, int(m.sum())


def report(run, win):
    series, rows = load_series(run)
    T, R = expansion_law(series)
    # frame mais proximo do fim da janela limpa
    t_ref = win[1]
    tf, fn = min(series, key=lambda s: abs(s[0] - t_ref))
    d = arrays(fn)
    x, y, rb = d["x"], d["y"], d["rho_b_grown"]
    r = np.hypot(x, y)
    Rc = float(np.percentile(r[rb > 0.1], 99))
    h0 = float(d["h"][0])

    W = 78
    print("=" * W)
    print(
        f"VALIDACAO DO MODELO — {os.path.basename(run)}   (frame de referencia t={tf:.1f}s)"
    )
    print("=" * W)

    # ---------- 1. LEI DE EXPANSAO ----------
    print("\n1. LEI DE EXPANSAO  R(t) ~ t^alpha")
    print(
        "   [T2] Srinivasan 2019: swarming nutrient-rich -> velocidade CONSTANTE, alpha=1"
    )
    print("   [T3] Giverso 2016:    dedos difusao-limitados -> alpha ~ 0.45")
    a1, r1, n1 = fit_alpha(T, R, win[0], win[1])
    a2, r2_, n2 = fit_alpha(T, R, 0.6 * win[1], win[1])
    a3, r3, n3 = fit_alpha(T, R, win[1], T.max())
    print(
        f"     janela util   t={win[0]:.0f}-{win[1]:.0f}   alpha={a1:.3f}  R2={r1:.4f}  n={n1}"
    )
    print(
        f"     estabelecido  t={0.6 * win[1]:.0f}-{win[1]:.0f}   alpha={a2:.3f}  R2={r2_:.4f}  n={n2}"
    )
    if n3 >= 4:
        print(
            f"     pos-esgotamento t={win[1]:.0f}-{T.max():.0f} alpha={a3:.3f}  R2={r3:.4f}  n={n3}"
        )
    print(
        f"   -> regime nutrient-rich reproduz [T2] (alpha~1): {'SIM' if abs(a2 - 1) < 0.15 else 'nao'}"
    )
    if n3 >= 4:
        print(
            f"   -> apos esgotar nutriente migra p/ [T3] (alpha~0.45): "
            f"{'SIM' if abs(a3 - 0.45) < 0.15 else 'nao'}"
        )

    # ---------- 2. MORFOLOGIA ----------
    print("\n2. MORFOLOGIA DENDRITICA")
    print("   reference.jpg (PA14, Michiels): 15-20 dendritos radiais, AR >= 1:5")
    print("   [T1] painel (b) Fingering: 7-9 dedos finos de nucleo compacto")
    fr = front_radius(x, y, rb)
    nf = count_fingers(fr)
    r_nuc = float(r[rb >= 0.8].max()) if (rb >= 0.8).any() else 0.0
    # largura: extensao angular ocupada por cristas, em 3 raios
    ws = []
    for f_ in (0.5, 0.65, 0.8):
        rad = f_ * Rc
        sel = (np.abs(r - rad) < 1.5 * DX) & (d["rho"] > 1.05)
        if sel.sum() < 10:
            continue
        th = np.sort(np.arctan2(y[sel], x[sel]))
        gaps = np.diff(th)
        big = gaps[gaps > np.deg2rad(4)]
        if len(big):
            ws.append((2 * np.pi - big.sum()) * rad / len(big))
    w = float(np.mean(ws)) if ws else np.nan
    ar = (Rc - r_nuc) / w if w == w and w > 0 else np.nan
    print(f"     dendritos contados: {nf}          (alvo PA14: 15-20 | [T1]: 7-9)")
    print(f"     R colonia={Rc:.2f}  nucleo pinado={r_nuc:.2f}  largura braco={w:.3f}")
    print(f"     AR = {ar:.1f}                      (criterio §2.2: >= 5)")

    # ---------- 3. HALO DE SURFACTANTE ----------
    print("\n3. HALO DE SURFACTANTE")
    print("   [T1] painel (b): campo de surfactante EXTENDE alem da biomassa")
    cs = d["cs"]
    isf = d.get("is_filler")
    real = (
        (isf < 0.5)
        if isf is not None and isf.size == cs.size
        else np.ones_like(cs, bool)
    )
    r_bio = float(np.percentile(r[(rb > 0.05) & real], 99))
    sel_cs = (cs > 0.02) & real
    r_cs = float(np.percentile(r[sel_cs], 99)) if sel_cs.any() else 0.0
    print(f"     raio biomassa (rho_b>0.05, p99) = {r_bio:.2f}")
    print(f"     raio surfactante (cs>0.02, p99) = {r_cs:.2f}")
    print(
        f"     halo = {r_cs - r_bio:+.2f}  ({r_cs / max(r_bio, 1e-9):.2f}x)"
        f"   -> {'PRESENTE' if r_cs > r_bio else 'AUSENTE'}"
    )

    # ---------- 4. CONSISTENCIA SPH ----------
    print("\n4. CONSISTENCIA SPH")
    print("   [T6] Liu §3.3: ~20 vizinhos p/ gradiente confiavel, ~35 p/ Laplaciano")
    tree = cKDTree(np.column_stack([x, y]))
    rng = np.random.default_rng(0)
    lim = 0.93 * float(np.max(np.abs(np.concatenate([x, y]))))
    interior = (np.abs(x) < lim) & (np.abs(y) < lim)
    print(
        f"     esperado numa rede uniforme: pi*(2h)^2/dx^2 = "
        f"{np.pi * (2 * h0) ** 2 / DX**2:.0f}"
    )
    for nm, msk in [
        ("nucleo (rho_b>0.8)", rb > 0.8),
        ("bracos (rho_b 0.1-0.5)", (rb >= 0.1) & (rb < 0.5)),
        ("agar interior", (rb < 0.05) & interior),
    ]:
        idx = np.where(msk)[0]
        if idx.size == 0:
            continue
        if idx.size > 300:
            idx = rng.choice(idx, 300, replace=False)
        nb = np.array(
            [len(tree.query_ball_point([x[i], y[i]], 2 * h0)) - 1 for i in idx]
        )
        print(
            f"     {nm:24s} media={nb.mean():5.1f}  >=20:{(nb >= 20).mean():4.0%}"
            f"  >=35:{(nb >= 35).mean():4.0%}"
        )

    print("\n   [T7] Violeau §3.4/3.6: particao da unidade; sigma_a<0.85 ~ 15% de erro")
    sa = d["sigma_a"]
    col = rb > 0.1
    print(
        f"     colonia: media={sa[col].mean():.4f}  frac<0.85={(sa[col] < 0.85).mean():.1%}"
        f"  p05={np.percentile(sa[col], 5):.3f}"
    )

    # ---------- 5. COBERTURA ESPACIAL ----------
    print("\n5. COBERTURA ESPACIAL (§2.5 C1)")
    n = 700
    gg = np.linspace(-Rc, Rc, n)
    GX, GY = np.meshgrid(gg, gg)
    xlim = float(np.max(np.abs(x)))
    ins = (GX**2 + GY**2 <= Rc**2) & (np.abs(GX) <= xlim) & (np.abs(GY) <= xlim)
    dd, _ = tree.query(np.column_stack([GX[ins], GY[ins]]))
    for thr in (0.7, 1.0, 1.5):
        fr_ = (dd > thr * DX).mean()
        print(
            f"     area sem vizinho a <{thr:.1f}dx: {fr_:7.3%}"
            f"   ({fr_ * np.pi * Rc**2 / DX**2:6.0f} dx2)"
        )

    # ---------- 6. CONSERVACAO ----------
    print("\n6. CONSERVACAO E ESTABILIDADE")
    m0 = float(rows[0]["mass_total"])
    mf = float(rows[-1]["mass_total"])
    ap = np.array([float(q["a_pressure"]) for q in rows[2:]])
    print(
        f"     massa {m0:.1f} -> {mf:.1f}  ({(mf / m0 - 1) * 100:+.1f}%)"
        f"   [insercao adiciona celulas; deve ser LIMITADO]"
    )
    print(
        f"     a_pressure mediana={np.median(ap):.2f}  picos>4 em {(ap > 4).mean():.0%}"
        f" das amostras   (orcamento §8: <3-4)"
    )
    print("=" * W)
    return T, R, series


def main(argv):
    win = [5.0, 55.0]
    runs = []
    i = 0
    while i < len(argv):
        if argv[i] == "--window":
            win = [float(argv[i + 1]), float(argv[i + 2])]
            i += 3
        else:
            runs.append(argv[i])
            i += 1
    for run in runs or ["runs/C4_t100"]:
        report(run, win)


if __name__ == "__main__":
    main(sys.argv[1:])
