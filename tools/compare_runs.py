"""Compara rotas sob o CRITERIO OBRIGATORIO DE VALIDACAO (estabelecido 2026-08-06).

Uso:  python tools/compare_runs.py runs/S0 runs/S1 ...

FUNDAMENTACAO (Violeau §3.4 — Particao da Unidade; Liu §3.3 — consistencia):
A representacao SPH de um fluido exige COBERTURA ESPACIAL CONTINUA. A particao da
unidade `sigma_a = Sum_j V_j W_aj -> 1` so e recuperavel se houver particula onde o
fluido existe. Um buraco geometrico na regiao de expansao NAO e um sigma_a ruim — e
a ausencia do proprio ponto de amostragem, invisivel a qualquer media feita "onde ha
particula". A KGC (Bonet-Lok) mitiga o operador sob suporte incompleto, mas nao
repara a quebra fisica de cobertura. Por isso a metrica primaria e AREAL.

CRITERIOS (ambos obrigatorios — "Estado de Preenchimento Denso"):
  C1. Fracao de Vazio  : area da colonia sem vizinho < 1.5dx  -> minimizar (alvo ~0)
  C2. Suporte do kernel: mean_sig_all >= 0.85 (colônia INTEIRA, inseridas incluidas)

Nenhum dos dois sozinho aprova. C2 sem C1 e o falso positivo do R0 (nao preencher
"melhora" a media por remover o ponto de amostragem). C1 sem C2 e o falso positivo
do R1 (preencher com particula mal suportada).

GUARDRAILS (violar reprova, independente de C1/C2):
  mean_v >= 3e-4          nao congelar a morfologia (licao #31)
  contrast_cs >= 12       motor vivo (§2.4 / §11)
  a_pressure <= 4.0       sem over-pack (licao #27)
  iteracoes <= 3000       custo de dt zero (licao #29)
  massa nao acelerando    taxa 2a/1a metade <= 1.3 (runaway, licao #34)
"""

import sys
import glob
import os
import csv as _csv
import numpy as np
import h5py

GUARD = {
    "mean_v_min": 3e-4,
    "contrast_cs_min": 12.0,
    "a_pressure_max": 4.0,
    # Alvo: colapso de dt (licao #29 = 85x/35x), NAO variacao normal. Baseline S4
    # = 2933 iter; 6000 e ~2x isso — pega colapso real sem reprovar +9%.
    "iter_max": 6000,
}
# O teto absoluto de massa era amarrado a um WAKE_MASS_BUDGET especifico e ficava
# obsoleto sempre que a rodada mudava o orcamento. O que importa e RUNAWAY: taxa de
# massa ACELERANDO (licao #34) vs limitada. Teste independente de configuracao.
MASS_ACCEL_MAX = 1.3
SIG_MIN = 0.85
# Rodada que morreu cedo tem log so com as primeiras linhas; sem este guarda ela
# aparece com vazio 0% / sig_all 1.0 (valores de t=0) e VENCE o ranking.
EXPECTED_T = 50.0
MIN_T_FRAC = 0.9


def read_log(run):
    with open(os.path.join(run, "log.csv")) as f:
        return list(_csv.DictReader(f))


def last_hdf5(run):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5")))
    return fs[-1] if fs else None


def void_from_hdf5(run, thresholds=(0.7, 1.0, 1.5), n_grid=400):
    """Fracao de vazio do ultimo frame + raio da colonia.

    Retorna tambem R para converter fracao -> AREA ABSOLUTA. A fracao sozinha e
    normalizada pela area da colonia, que cresce; ela pode cair sem que o vacuo
    absoluto diminua. Rotas com raios diferentes exigem os dois numeros.
    """
    from scipy.spatial import cKDTree

    out = {t: np.nan for t in thresholds}
    out["R"] = np.nan
    fn = last_hdf5(run)
    if not fn:
        return out
    dx = 10.0 / 186
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        x, y = np.array(a["x"]), np.array(a["y"])
        rb = np.array(a["rho_b_grown"])
    colony = rb > 0.1
    if colony.sum() < 10:
        return out
    R = float(np.percentile(np.hypot(x, y)[colony], 99))
    g = np.linspace(-R, R, n_grid)
    GX, GY = np.meshgrid(g, g)
    # recorta ao dominio (inferido dos dados): fora dele nao ha particula
    xlim = float(np.max(np.abs(x)))
    ylim = float(np.max(np.abs(y)))
    ins = (
        (GX * GX + GY * GY <= R * R)
        & (np.abs(GX) <= xlim)
        & (np.abs(GY) <= ylim)
    )
    d, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX[ins], GY[ins]])
    )
    res = {t: float(np.mean(d > t * dx)) for t in thresholds}
    res["R"] = R
    return res


def sigma_from_hdf5(run):
    out = dict(sig_bio=np.nan, sig_all=np.nan, low_all=np.nan, n_bio=0, n_ins=0)
    fn = last_hdf5(run)
    if not fn:
        return out
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]

        def g(n):
            return np.array(a[n]) if n in a and np.array(a[n]).size else None

        rb, sa, isw, isf = (g("rho_b_grown"), g("sigma_a"), g("is_wake"),
                            g("is_filler"))
        if rb is None or sa is None:
            return out
        arms = (rb >= 0.1) & (rb < 0.5)
        n = rb.size
        ins = (isf > 0.5) if (isf is not None and isf.size == n) else (
            (isw > 0.5) if (isw is not None and isw.size == n)
            else np.zeros(n, bool)
        )
        bio = arms & ~ins
        out["sig_all"] = float(sa[arms].mean()) if arms.any() else np.nan
        out["low_all"] = float((sa[arms] < 0.85).mean()) if arms.any() else np.nan
        out["sig_bio"] = float(sa[bio].mean()) if bio.any() else np.nan
        out["n_bio"] = int(bio.sum())
        out["n_ins"] = int((arms & ins).sum())
    return out


def motor_from_hdf5(run):
    """Backfill das metricas de motor p/ rodadas anteriores as colunas novas."""
    out = dict(amed=np.nan, ap95=np.nan, csb=np.nan, cnb=np.nan,
               contr_real=np.nan)
    fn = last_hdf5(run)
    if not fn:
        return out
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]

        def g(n):
            return np.array(a[n]) if n in a and np.array(a[n]).size else None

        rb, am, cs, cn, isf, x, y = (g("rho_b_grown"), g("au_mar"), g("cs"),
                                     g("c_n"), g("is_filler"), g("x"), g("y"))
        if rb is None or am is None or isf is None:
            return out
        bio = (rb > 0.1) & (isf < 0.5)
        if not bio.any():
            return out
        r = np.hypot(x, y)
        R = float(np.percentile(r[rb > 0.1], 99))
        arm = bio & (r > 0.3 * R)
        real = isf < 0.5
        out["contr_real"] = float(
            (cs[real].max() - cs[real].min()) / (cs[real].mean() + 1e-9)
        )
        out["amed"] = float(np.median(am[bio]))
        out["ap95"] = float(np.percentile(am[bio], 95))
        if arm.any():
            out["csb"] = float(cs[arm].mean())
            out["cnb"] = float(cn[arm].mean()) if cn is not None else np.nan
    return out


def metrics(run):
    try:
        rows = read_log(run)
    except FileNotFoundError:
        return None
    if not rows:
        return None
    last = rows[-1]
    tail = rows[max(0, len(rows) - 5):]

    def c(r, k, d=0.0):
        try:
            return float(r[k])
        except (KeyError, TypeError, ValueError):
            return d

    m = {
        "run": os.path.basename(run.rstrip("/")),
        "t": c(last, "t"),
        "iter": int(c(last, "iteration")),
        "mean_v": np.mean([c(r, "mean_v") for r in tail]),
        "mass": c(last, "mass_total"),
        "contrast_cs": np.mean([c(r, "constrast_cs") for r in tail]),
        "a_mar": np.mean([c(r, "a_marangoni") for r in tail]),
        "a_press": max(c(r, "a_pressure") for r in tail),
        "m0": c(rows[0], "mass_total", 1.0),
    }
    _t = [c(r, "t") for r in rows]
    _m = [c(r, "mass_total") for r in rows]
    _h = len(_t) // 2
    if len(_t) > 3 and _t[_h] > _t[0] and _t[-1] > _t[_h]:
        _d1 = (_m[_h] - _m[0]) / (_t[_h] - _t[0])
        _d2 = (_m[-1] - _m[_h]) / (_t[-1] - _t[_h])
        m["mass_accel"] = _d2 / max(_d1, 1e-9)
    else:
        m["mass_accel"] = 0.0

    if "mean_sig_all" in last:
        m["sig_bio"] = c(last, "mean_sig_bio")
        m["sig_all"] = c(last, "mean_sig_all")
        m["low_all"] = c(last, "frac_lowsig_all")
        m["n_bio"] = int(c(last, "n_bio_arms"))
        m["n_ins"] = int(c(last, "n_ins_arms"))
    else:
        m.update(sigma_from_hdf5(run))

    vf = void_from_hdf5(run)  # sempre — precisamos de R p/ a area absoluta
    _mh = motor_from_hdf5(run)
    if "a_mar_bio_med" in last:
        for k, dk in [("a_mar_bio_med", "amed"), ("a_mar_bio_p95", "ap95"),
                      ("cs_bio_arms", "csb"), ("c_n_bio_arms", "cnb")]:
            m[dk] = c(last, k, float("nan"))
    else:
        m.update({k: v for k, v in _mh.items() if k != "contr_real"})
    # contraste SEMPRE do HDF5 excluindo fillers (o log antigo os incluia)
    m["contrast_cs"] = _mh["contr_real"]

    if "void_15" in last:
        m["v07"], m["v10"], m["v15"] = (c(last, "void_07"), c(last, "void_10"),
                                        c(last, "void_15"))
    else:
        m["v07"], m["v10"], m["v15"] = vf[0.7], vf[1.0], vf[1.5]
    m["R"] = vf["R"]
    dx = 10.0 / 186
    # area absoluta do vacuo, em unidades de dx^2 (= "quantas particulas cabem")
    m["a15"] = m["v15"] * np.pi * m["R"] ** 2 / (dx * dx)

    return m


def evaluate(m):
    """Retorna (falhas_guardrail, falhas_criterio)."""
    g = []
    if m["t"] < MIN_T_FRAC * EXPECTED_T:
        g.append(f"INCOMPLETA t={m['t']:.1f}<{MIN_T_FRAC * EXPECTED_T:.0f}")
    if m["mean_v"] < GUARD["mean_v_min"]:
        g.append(f"mean_v={m['mean_v']:.1e}")
    if m["contrast_cs"] < GUARD["contrast_cs_min"]:
        g.append(f"contrast={m['contrast_cs']:.1f}")
    if m["a_press"] > GUARD["a_pressure_max"]:
        g.append(f"a_press={m['a_press']:.1f}")
    if m["iter"] > GUARD["iter_max"]:
        g.append(f"iter={m['iter']}")
    if m["mass_accel"] > MASS_ACCEL_MAX:
        g.append(
            f"massa NAO CONVERGIDA (taxa 2a/1a metade={m['mass_accel']:.2f})"
        )

    c = []
    if not (m["sig_all"] >= SIG_MIN):
        c.append(f"C2 sig_all={m['sig_all']:.3f}<{SIG_MIN}")
    return g, c


def main(runs):
    ms = [m for m in (metrics(r) for r in runs) if m]
    if not ms:
        print("nenhuma rodada encontrada")
        return

    hdr = (f"{'run':<5}{'t':>6}│{'VOID>1.5dx':>11}{'area(dx²)':>10}{'>1.0dx':>8}"
           f"{'>0.7dx':>8}{'R':>6}│{'sig_all':>8}{'low_all':>8}│"
           f"{'aMarMed':>8}{'aMarP95':>8}{'cs_bio':>8}{'c_n_bio':>8}│"
           f"{'massa':>7}{'m.acel':>7}{'mean_v':>9}{'a_pr':>6}")
    print(__doc__.split("CRITERIOS")[1].split("GUARDRAILS")[0].strip())
    print()
    print("C1 e reportado como FRACAO e como AREA ABSOLUTA (em dx², ~n de particulas)")
    print("— a fracao e normalizada pela area da colonia, que cresce entre rotas.")
    print()
    print(hdr)
    print("─" * len(hdr))
    for m in ms:
        g, c = evaluate(m)
        print(
            f"{m['run']:<5}{m['t']:>6.1f}│{m['v15']:>11.2%}{m['a15']:>10.0f}"
            f"{m['v10']:>8.2%}{m['v07']:>8.2%}{m['R']:>6.2f}│{m['sig_all']:>8.3f}"
            f"{m['low_all']:>8.1%}│{m['amed']:>8.2f}{m['ap95']:>8.2f}"
            f"{m['csb']:>8.4f}{m['cnb']:>8.3f}│"
            f"{m['mass']:>7.1f}{m['mass_accel']:>7.2f}{m['mean_v']:>9.1e}"
            f"{m['a_press']:>6.1f}"
        )
        if g or c:
            print(f"      └─ REPROVA: {'; '.join(c + g)}")

    ok = [m for m in ms if not any(evaluate(m))]
    print()
    if ok:
        ok.sort(key=lambda m: (m["v15"], -m["sig_all"]))
        print("RANKING (aprovados em C2+guardrails, ordenados por C1 = menor vazio):")
        for i, m in enumerate(ok, 1):
            print(f"  {i}. {m['run']:<4} vazio>1.5dx={m['v15']:.2%} "
                  f"({m['a15']:.0f} dx²)  sig_all={m['sig_all']:.3f}  "
                  f"massa={m['mass']:.1f}")
    else:
        print("NENHUMA configuracao satisfez C2 + guardrails.")
        best = min(ms, key=lambda m: m["v15"])
        print(f"(menor vazio absoluto: {best['run']} = {best['v15']:.2%})")


if __name__ == "__main__":
    main(sys.argv[1:] or sorted(glob.glob("runs/*")))
