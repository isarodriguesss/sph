"""Diagnostico do ALVO DE RECRUTAMENTO de biomassa (Passo 0.1 do PLANO_K3_JUNCAO).

Uso:  python tools/diag_recrut.py runs/C4 [runs/K3 ...] [--t 48]
      python tools/diag_recrut.py .              # rodada corrente (main_output/ na raiz)

Complementa `diag_juncao.py`, que mede o DEFEITO (vale `V` na crista) e o MOTOR
(`a_mar_front`). Este mede a ALAVANCA: para onde um mecanismo de recrutamento
levaria as particulas com `rho_b = 0`.

A `BiomassColonization` relaxa para a media SHEPARD da vizinhanca
(`rho_b_smooth / sigma_a`). Num campo 92% vazio essa media e a media entre a mae e
o vacuo — mede o vacuo, nao a mae. O alvo biologicamente correto e a densidade de
QUEM DOA (a filha nasce com a densidade da mae, corolario da licao #48), que e a
media de `rho_b` PONDERADA POR `rho_b`:

    rho_b_doador = sum_j V_j rho_b_j^2 W_ij / sum_j V_j rho_b_j W_ij

Medido no C4 (t=21.6), nos buracos: Shepard da 0.006-0.025 e o doador da
0.185-0.381 — 7x a 60x. O primeiro recruta ABAIXO do quorum 0.1, onde a particula
e mecanicamente invisivel (EOS `fade=0`, gate flagelar [0.1,0.6], `ParticleShift`
>= 0.1) e so quimicamente ativa: e a licao #53.

FILLER COMO DOADOR (§5 do plano): `BiomassColonization.loop` soma todas as fontes
sem checar `s_is_filler`, e o filler entra com `rho_b >= 0.5` herdado e CONGELADO.
Sob o alvo-doador ele vira doador forte e imovel — valor historico que nao decai,
a armadilha da licao #39. Por isso as duas colunas `com/sem filler`.

Le direto do HDF5; nao depende de colunas de log.
"""

import argparse
import glob
import os

import h5py
import numpy as np
from scipy.spatial import cKDTree

R_MAX = 3.6
HOLE_RHO_B = 0.01
QUORUM = 0.1
RINGS = np.arange(0.0, 3.2, 0.4)
DONOR_GATES = (0.1, 0.2, 0.3, 0.5)


def frame_time(f):
    sd = f["solver_data"]
    if "t" in sd.attrs:
        return float(sd.attrs["t"])
    return float(np.array(sd["t"]))


def pick_frame(run, t_target):
    out = os.path.join(run, "main_output")
    if not os.path.isdir(out):
        out = run
    files = sorted(glob.glob(os.path.join(out, "*.hdf5")))
    if not files:
        raise SystemExit(f"sem HDF5 em {out}")
    if t_target is None:
        return files[-1]
    best, best_d = files[-1], 1e30
    for fn in files:
        with h5py.File(fn, "r") as f:
            t = frame_time(f)
        if abs(t - t_target) < best_d:
            best, best_d = fn, abs(t - t_target)
    return best


def load(fn):
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        d = {
            k: a[k][:]
            for k in (
                "x",
                "y",
                "m",
                "rho",
                "h",
                "rho_b_grown",
                "cs",
                "c_n",
                "is_filler",
            )
        }
        d["t"] = frame_time(f)
    d["r"] = np.hypot(d["x"], d["y"])
    return d


def cubic_spline(rij, h):
    q = rij / h
    sig = 10.0 / (7.0 * np.pi * h * h)
    w = np.zeros_like(q)
    a = q <= 1.0
    b = (q > 1.0) & (q <= 2.0)
    w[a] = sig * (1.0 - 1.5 * q[a] ** 2 + 0.75 * q[a] ** 3)
    w[b] = sig * 0.25 * (2.0 - q[b]) ** 3
    return w


def targets(d):
    """Alvo Shepard e alvo-doador (com e sem filler) para cada particula em r<R_MAX."""
    x, y, rb = d["x"], d["y"], d["rho_b_grown"]
    h = float(np.median(d["h"]))
    vol = d["m"] / np.maximum(d["rho"], 1e-9)
    real = d["is_filler"] < 0.5

    sel = np.where(d["r"] < R_MAX)[0]
    tree = cKDTree(np.column_stack([x, y]))
    nbrs = tree.query_ball_point(np.column_stack([x[sel], y[sel]]), r=2.0 * h)

    n = len(x)
    shep = np.zeros(n)
    donor = np.zeros(n)
    donor_real = np.zeros(n)

    for k, i in enumerate(sel):
        j = np.asarray(nbrs[k], dtype=int)
        w = cubic_spline(np.hypot(x[j] - x[i], y[j] - y[i]), h)
        vw = vol[j] * w
        rbj = rb[j]

        den_s = vw.sum()
        if den_s > 1e-12:
            shep[i] = (vw * rbj).sum() / den_s

        den_d = (vw * rbj).sum()
        if den_d > 1e-12:
            donor[i] = (vw * rbj * rbj).sum() / den_d

        rj = real[j]
        den_r = (vw[rj] * rbj[rj]).sum()
        if den_r > 1e-12:
            donor_real[i] = (vw[rj] * rbj[rj] * rbj[rj]).sum() / den_r

    return shep, donor, donor_real


def report(run, t_target):
    fn = pick_frame(run, t_target)
    d = load(fn)
    rb, r, cs = d["rho_b_grown"], d["r"], d["cs"]
    real = d["is_filler"] < 0.5
    bio = (rb > QUORUM) & real
    R99 = float(np.percentile(r[bio], 99)) if bio.any() else 0.0
    disc = (r < max(R99, 1e-6)) & real

    print(f"\n{'=' * 78}")
    print(f"{run}   {os.path.basename(fn)}   t = {d['t']:.2f} s")
    print(
        f"N={len(rb)} ({int(real.sum())} reais + {int((~real).sum())} filler)  "
        f"n_bio(rho_b>{QUORUM}, real)={int(bio.sum())}  "
        f"filler com rho_b>{QUORUM}={int(((rb > QUORUM) & ~real).sum())}  "
        f"R99={R99:.2f}  zeros no disco (reais)={np.mean(rb[disc] <= 1e-12):.1%}"
    )

    shep, donor, donor_real = targets(d)
    hole = (rb < HOLE_RHO_B) & (r < R_MAX) & real

    print("\n[A] ALVO DE RECRUTAMENTO nos buracos (rho_b < 0.01), por anel")
    print(
        f"{'r':>5} {'n_hole':>7} {'shepard':>9} {'doador':>9} {'razao':>7} "
        f"{'doador s/filler':>16} {'>quorum?':>9}"
    )
    for r0 in RINGS:
        m = hole & (r >= r0) & (r < r0 + 0.4)
        if m.sum() < 5:
            continue
        s, dn, dr = np.median(shep[m]), np.median(donor[m]), np.median(donor_real[m])
        ratio = dn / s if s > 1e-9 else np.inf
        print(
            f"{r0 + 0.2:5.1f} {int(m.sum()):7d} {s:9.4f} {dn:9.4f} "
            f"{ratio:7.1f} {dr:16.4f} "
            f"{'shep NAO' if s < QUORUM else 'ambos':>9}"
        )

    print("\n[B] RECRUTAVEIS por gate de doador (buracos em r < R99)")
    hd = hole & (r < R99)
    for g in DONOR_GATES:
        for name, field in (
            ("doador", donor),
            ("dr s/filler", donor_real),
            ("shepard", shep),
        ):
            s = hd & (field > g)
            tgt = np.median(field[s]) if s.any() else 0.0
            print(
                f"  {name:>7} > {g:.1f}: {int(s.sum()):6d} buracos   "
                f"alvo mediano = {tgt:.3f}   "
                f"{'ACIMA do quorum' if tgt > QUORUM else 'sub-quorum'}"
            )

    print("\n[C] OCUPACAO por anel — populacao REAL (filler contado a parte)")
    print(
        f"{'r':>5} {'n_real':>7} {'frac>quorum':>12} {'frac==0':>9} "
        f"{'rho_b p90':>10} {'n_filler':>9} {'fil>quorum':>11} "
        f"{'cs med':>8} {'cs/cs_max':>10} {'c_n med':>8}"
    )
    cs_max = float(np.max(cs)) if cs.size else 1.0
    for r0 in RINGS:
        m = (r >= r0) & (r < r0 + 0.4) & real
        mf = (r >= r0) & (r < r0 + 0.4) & ~real
        if not m.any():
            continue
        print(
            f"{r0 + 0.2:5.1f} {int(m.sum()):7d} {np.mean(rb[m] > QUORUM):12.3f} "
            f"{np.mean(rb[m] <= 1e-12):9.3f} {np.percentile(rb[m], 90):10.3f} "
            f"{int(mf.sum()):9d} {int(((rb > QUORUM) & mf).sum()):11d} "
            f"{np.median(cs[m]):8.4f} {np.median(cs[m]) / cs_max:10.3f} "
            f"{np.median(d['c_n'][m]):8.3f}"
        )

    return {
        "run": run,
        "t": d["t"],
        "zeros": float(np.mean(rb[disc] <= 1e-12)),
        "n_bio": int(bio.sum()),
        "junc_donor": float(np.median(donor[hole & (r >= 0.4) & (r < 1.2)]))
        if (hole & (r >= 0.4) & (r < 1.2)).any()
        else 0.0,
        "junc_shep": float(np.median(shep[hole & (r >= 0.4) & (r < 1.2)]))
        if (hole & (r >= 0.4) & (r < 1.2)).any()
        else 0.0,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", default=["."])
    ap.add_argument("--t", type=float, default=None, help="instante alvo (s)")
    args = ap.parse_args()

    rows = [report(run, args.t) for run in (args.runs or ["."])]

    print(f"\n{'=' * 78}\nRESUMO (junção r em [0.4, 1.2), buracos)")
    print(
        f"{'run':<16} {'t':>7} {'zeros':>8} {'n_bio':>7} {'shepard':>9} {'doador':>9}"
    )
    for m in rows:
        print(
            f"{m['run']:<16} {m['t']:7.1f} {m['zeros']:8.1%} {m['n_bio']:7d} "
            f"{m['junc_shep']:9.4f} {m['junc_donor']:9.4f}"
        )
    print(
        "\nalvo Shepard abaixo de 0.1 = recruta sub-quorum (licao #53). "
        "O doador e o alvo correto."
    )


if __name__ == "__main__":
    main()
