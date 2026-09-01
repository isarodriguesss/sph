"""Veredito das rotas contra o topico 2 (licao #66): agar engolido nao vira colonia.

    python tools/verdict_topico2.py runs/E1_eosfix runs/C4 [outro...]

Imprime, no ultimo frame de cada run, as metricas pre-registradas: fracao de area
da colonia ocupada por agar morto, R99, coesao efetiva na junção, e os guardrails
(contrast_cs, a_pressure, massa, nutriente, pin).
"""
import csv
import glob
import sys

import numpy as np
from pysph.solver.utils import load

DX = 0.0538


def fade(v):
    """fade_rep = fade_att da BiomassEOS: smoothstep em rho_b [0.1, 0.5]."""
    s = np.zeros_like(v)
    m = (v >= 0.1) & (v < 0.5)
    t = (v[m] - 0.1) / 0.4
    s[m] = t * t * (3 - 2 * t)
    s[v >= 0.5] = 1.0
    return s


def medir(run):
    fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
    if not fs:
        return None
    d = load(fs[-1])
    pa = d["arrays"]["fluid"]
    r = np.hypot(pa.x, pa.y)
    rb = pa.rho_b_grown
    fil = pa.is_filler > 0.5
    R99 = np.percentile(r[rb > 0.1], 99)

    agar = (r < R99) & (~fil) & (rb < 1e-12)
    junc = (r >= 0.35) & (r < 0.65)
    occ = junc.sum() * DX * DX / (np.pi * (0.65**2 - 0.35**2))

    lg = list(csv.DictReader(open(f"{run}/log.csv")))[-1]
    return {
        "t": d["solver_data"]["t"],
        "R99": R99,
        "%agar": 100 * agar.sum() * DX * DX / (np.pi * R99 * R99),
        "n_agar": int(agar.sum()),
        "coes": occ * np.mean(fade(rb[junc])),
        "viva_j": int((junc & (~fil) & (rb >= 0.1)).sum()),
        "contr_cs": float(lg["constrast_cs"]),
        "a_press": float(lg["a_pressure"]),
        "massa": float(lg["mass_total"]),
        "min_c_n": float(lg["min_c_n"]),
        "n_pin": int(lg["n_pinned"]),
    }


COLS = ["t", "R99", "%agar", "n_agar", "coes", "viva_j",
        "contr_cs", "a_press", "massa", "min_c_n", "n_pin"]
ALVO = {"%agar": "< 60", "R99": ">= 3.9", "contr_cs": ">= 11",
        "a_press": "<= 4", "massa": "<= 210"}


def main(runs):
    print(f"{'run':>16} " + " ".join(f"{c:>9}" for c in COLS))
    print(f"{'ALVO':>16} " + " ".join(f"{ALVO.get(c,''):>9}" for c in COLS))
    for run in runs:
        v = medir(run)
        if v is None:
            print(f"{run.split('/')[-1]:>16}   (sem frames)")
            continue
        cells = []
        for c in COLS:
            x = v[c]
            cells.append(f"{x:>9.2f}" if isinstance(x, float) else f"{x:>9d}")
        print(f"{run.split('/')[-1]:>16} " + " ".join(cells))


if __name__ == "__main__":
    main(sys.argv[1:] or ["runs/C4", "runs/E1_eosfix"])
