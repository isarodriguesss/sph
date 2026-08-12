"""Diagnostico do VALE DE BIOMASSA na juncao nucleo-braco.

Uso:  python tools/diag_juncao.py runs/C4 runs/K2 [--t 48] [--out FIG.png]

O defeito nao e vacuo de PARTICULA (C1/C2 de §2.5 ja estao satisfeitos: void15
0.16%, sigma_a estavel). E um vale no campo escalar `rho_b`: a biomassa vale ~1.0
no nucleo, cai para ~0.3 num anel intermediario e volta a ~0.55 nos bracos. A
colonia e materialmente continua e biologicamente descontinua.

Por que a metrica e o perfil AO LONGO DA CRISTA e nao a media azimutal: a media
mistura braco com baia e reporta um vale que existe so porque a baia esta vazia
por construcao (e a morfologia dendritica que se quer). O perfil da crista segue
o maximo de `rho_b` dentro do cone angular de cada braco dominante — ali um vale
e um vale de verdade.

METRICA PRIMARIA — profundidade normalizada do vale:

    V = 1 - min_r( rho_b_crista(r) ) / mediana( rho_b_crista(r > r_arm) )

V = 0 significa crista sem depressao; V = 0.5 significa que o vale tem metade da
biomassa do braco. Normalizar pelo braco (e nao por 1.0) e o que torna a metrica
insensivel ao nivel absoluto de biomassa — sem isso, qualquer alavanca que apenas
INFLE a biomassa global "melhora" o vale sem reconectar nada (foi assim que a K2
pontuou bem em enchimento enquanto matava o motor).

GUARDRAIL OBRIGATORIO — `a_mar_bio_med >= 2.5` (§2.5.1 C3, baseline C4 = 2.89).
A forca e `-beta*grad(cs)`: encher biomassa satura `cs` (produção `qs(rho_b)` sobe
contra o teto `(1-cs/cs_max)`), o gradiente some e o motor morre. Medido na K2:
biomassa 3.6x -> `cs/cs_max` 0.885->0.975 na juncao -> `a_mar_bio_med` 2.89->0.13.
Qualquer rodada que baixe o vale as custas do motor esta REPROVADA.

Le direto do HDF5 — nao depende das colunas de log, que so existem nas rodadas
posteriores a 2026-08-11.
"""

import sys
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import h5py

ARM_R0, ARM_R1 = 2.0, 3.0
ARM_HALF_ANGLE = np.deg2rad(9)
N_ARMS = 4
R_BINS = np.arange(0.3, 4.6, 0.2)
JUNC_LO, JUNC_HI = 0.4, 1.2
CS_MAX = 0.5
PIN_RHO_B = 0.8
GATE_LO, GATE_HI = 0.4, 0.8
A_MAR_MIN = 2.5
A_MAR_FRONT_MIN = 4.0


def run_cs_max(run):
    """cs_max da rodada, lido do log — a saturacao de T1 crava `max_cs` no teto.

    Nao pode ser constante: a serie J varia `cs_max` entre rodadas, e normalizar
    tudo por 0.5 inverteria a comparacao.
    """
    path = os.path.join(run, "log.csv")
    if not os.path.exists(path):
        return CS_MAX
    try:
        with open(path) as fh:
            head = fh.readline().strip().split(",")
            i = head.index("max_cs")
            vals = [float(ln.split(",")[i]) for ln in fh if ln.strip()]
        return max(vals) if vals else CS_MAX
    except (ValueError, IndexError):
        return CS_MAX


def cn_gate(c):
    t = np.clip((c - GATE_LO) / (GATE_HI - GATE_LO), 0.0, 1.0)
    return np.where(c < GATE_LO, 0.0, t * t * (3.0 - 2.0 * t))


def frame_time(f):
    sd = f["solver_data"]
    if "t" in sd.attrs:  # PySPH grava t como atributo; versoes antigas, como dataset
        return float(sd.attrs["t"])
    return float(np.array(sd["t"]))


def pick_frame(run, t_target):
    files = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
    if not files:
        raise SystemExit(f"sem HDF5 em {run}/main_output")
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
            for k in ("x", "y", "rho_b_grown", "cs", "c_n", "is_filler", "au_mar")
        }
        d["t"] = frame_time(f)
    d["r"] = np.hypot(d["x"], d["y"])
    d["th"] = np.arctan2(d["y"], d["x"])
    return d


def arm_azimuths(d):
    sel = (d["r"] > ARM_R0) & (d["r"] < ARM_R1) & (d["rho_b_grown"] > 0.3)
    if int(sel.sum()) < 8:
        sel = (d["r"] > ARM_R0) & (d["r"] < ARM_R1) & (d["rho_b_grown"] > 0.1)
    if int(sel.sum()) < 4:
        return []
    h, e = np.histogram(d["th"][sel], bins=72, range=(-np.pi, np.pi))
    return [0.5 * (e[b] + e[b + 1]) for b in np.argsort(h)[-N_ARMS:] if h[b] > 0]


def ridge_profile(d, azimuths):
    """Maximo de rho_b dentro do cone de cada braco, media entre bracos."""
    prof = np.full(len(R_BINS), np.nan)
    if not azimuths:
        return prof
    per_arm = []
    for c in azimuths:
        dth = np.abs(((d["th"] - c + np.pi) % (2 * np.pi)) - np.pi)
        cone = dth < ARM_HALF_ANGLE
        vals = []
        for r0 in R_BINS:
            m = cone & (np.abs(d["r"] - r0) < 0.1)
            vals.append(np.max(d["rho_b_grown"][m]) if m.any() else np.nan)
        per_arm.append(vals)
    with np.errstate(all="ignore"):
        prof = np.nanmean(np.array(per_arm, dtype=float), axis=0)
    return prof


def radial_stat(d, mask, field, agg=np.median):
    out = np.full(len(R_BINS), np.nan)
    for i, r0 in enumerate(R_BINS):
        m = mask & (np.abs(d["r"] - r0) < 0.1)
        if m.sum() >= 2:
            out[i] = agg(field[m])
    return out


def analyse(run, t_target):
    fn = pick_frame(run, t_target)
    d = load(fn)
    cs_max = run_cs_max(run)
    real = d["is_filler"] < 0.5
    bio = (d["rho_b_grown"] > 0.1) & real
    az = arm_azimuths(d)
    prof = ridge_profile(d, az)

    arm_level = np.nanmedian(prof[R_BINS >= ARM_R0])
    valley_zone = (R_BINS >= JUNC_LO) & (R_BINS < ARM_R0)
    with np.errstate(all="ignore"):
        valley = np.nanmin(prof[valley_zone])
    # Sem crista mensuravel (colonia homogenea de baixa densidade) `V` nao tem
    # significado — a razao explode. Reportar NaN em vez de um numero inventado.
    if not np.isfinite(arm_level) or arm_level < 0.05 or not np.isfinite(valley):
        V = np.nan
    else:
        V = 1.0 - valley / arm_level

    ju = (d["r"] >= JUNC_LO) & (d["r"] < JUNC_HI)
    pinned = (d["rho_b_grown"] >= PIN_RHO_B) | (d["c_n"] < 0.6) | (d["is_filler"] > 0.5)

    # Continuidade do corpo: fracao do disco da colonia que carrega biomassa.
    # `frac_zero` isola o estado absorvente — rho_b EXATAMENTE 0 nao cresce nunca,
    # porque BiomassGrowth e multiplicativo (rate*rho_b).
    R99 = float(np.percentile(d["r"][d["rho_b_grown"] > 0.1], 99))
    disc = d["r"] < R99
    n_disc = max(int(disc.sum()), 1)
    F = float(np.sum(disc & (d["rho_b_grown"] > 0.01)) / n_disc)
    frac_zero = float(np.sum(disc & (d["rho_b_grown"] == 0.0)) / n_disc)

    return {
        "F": F,
        "frac_zero": frac_zero,
        "run": os.path.basename(run.rstrip("/")),
        "file": os.path.basename(fn),
        "t": d["t"],
        "d": d,
        "bio": bio,
        "prof": prof,
        "V": V,
        "valley": valley,
        "arm_level": arm_level,
        "n_bio": int(bio.sum()),
        "R99": float(np.percentile(d["r"][d["rho_b_grown"] > 0.1], 99)),
        "rho_b_junc": float(np.percentile(d["rho_b_grown"][ju], 90)),
        "c_n_junc": float(np.mean(d["c_n"][ju])),
        "a_mar_med": float(np.median(d["au_mar"][bio])),
        "a_mar_p95": float(np.percentile(d["au_mar"][bio], 95)),
        # Licao #48: a mediana sobre TODA a biomassa nao e comparavel entre rodadas
        # que mudam a populacao (167 -> 712 particulas). A frente e o motor real.
        "a_mar_front": float(np.median(d["au_mar"][bio & (d["r"] > 0.75 * R99)]))
        if int((bio & (d["r"] > 0.75 * R99)).sum()) > 2
        else float("nan"),
        "pin_frac": float(np.mean(pinned[bio])),
        "cs_max": cs_max,
        "cs_junc": float(np.median(d["cs"][ju]) / cs_max),
        "prof_cn": radial_stat(d, real, d["c_n"]),
        "prof_cs": radial_stat(d, real, d["cs"]) / cs_max,
        "prof_amar": radial_stat(d, bio, d["au_mar"]),
    }


def main():
    args = [a for a in sys.argv[1:]]
    t_target, out = None, "runs/juncao.png"
    if "--t" in args:
        i = args.index("--t")
        t_target = float(args[i + 1])
        del args[i : i + 2]
    if "--out" in args:
        i = args.index("--out")
        out = args[i + 1]
        del args[i : i + 2]
    if not args:
        raise SystemExit(__doc__)

    res = [analyse(r, t_target) for r in args]

    print(
        f"\n{'rodada':>12} {'t':>6} {'F(rb>.01)':>10} {'rb==0':>7} {'R99':>5} "
        f"{'cs/max':>7} {'a_mar_med':>10} {'a_mar_FRENTE':>13} {'n_bio':>6} {'veredito':>10}"
    )
    for r in res:
        motor_ko = r["a_mar_front"] < A_MAR_FRONT_MIN
        drown = r["cs_junc"] > 0.85
        ok = (r["F"] >= 0.60) and not motor_ko and not drown
        verd = (
            "APROVA"
            if ok
            else ("MOTOR KO" if motor_ko else ("AFOGOU" if drown else "halo"))
        )
        print(
            f"{r['run']:>12} {r['t']:6.1f} {100 * r['F']:9.1f}% {100 * r['frac_zero']:6.1f}% "
            f"{r['R99']:5.2f} {r['cs_junc']:7.3f} {r['a_mar_med']:10.3f} "
            f"{r['a_mar_front']:13.3f} {r['n_bio']:6d} {verd:>10}"
        )
    print(
        f"\nalvo: F >= 60%  E  a_mar_FRENTE >= {A_MAR_FRONT_MIN}  E  cs/cs_max(junção) <= 0.85"
        f"   (J0 @t=48: F=11.8%, frente=7.23, cs/max=0.73)\n"
    )

    n = len(res)
    fig = plt.figure(figsize=(4.6 * max(n, 2), 10.5))
    gs = fig.add_gridspec(
        3, max(n, 2), height_ratios=[1.15, 1.0, 1.25], hspace=0.34, wspace=0.24
    )

    ax = fig.add_subplot(gs[0, :])
    for r in res:
        ax.plot(
            R_BINS,
            r["prof"],
            "o-",
            ms=3.5,
            lw=1.6,
            label=f"{r['run']}  V={r['V']:.2f}  a_mar={r['a_mar_med']:.2f}",
        )
        if np.isfinite(r["valley"]):
            i = int(
                np.nanargmin(
                    np.where((R_BINS >= JUNC_LO) & (R_BINS < ARM_R0), r["prof"], np.inf)
                )
            )
            ax.plot(
                R_BINS[i],
                r["prof"][i],
                "v",
                ms=11,
                mfc="none",
                mew=2,
                color=ax.lines[-1].get_color(),
            )
    ax.axvspan(JUNC_LO, ARM_R0, color="0.85", zorder=0)
    ax.axhline(0.8, ls=":", c="k", lw=1)
    ax.text(0.42, 0.82, "limiar do pin (rho_b=0.8)", fontsize=8)
    ax.text(0.45, 0.03, "zona da juncao", fontsize=9, color="0.35")
    ax.set_xlabel("r")
    ax.set_ylabel(r"$\rho_b$ na crista do braço")
    ax.set_title("Perfil de biomassa ao longo da crista — o vale é o marcador (▽)")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[1, 0])
    for r in res:
        ax.plot(R_BINS, r["prof_cn"], "-", lw=1.6, label=r["run"])
    ax.axhline(0.4, c="crimson", ls="--", lw=1.2)
    ax.axhline(0.6, c="darkorange", ls="--", lw=1.2)
    ax.text(3.0, 0.42, "gate de crescimento", fontsize=8, color="crimson")
    ax.text(3.0, 0.62, "pin cinemático", fontsize=8, color="darkorange")
    ax.axvspan(JUNC_LO, ARM_R0, color="0.85", zorder=0)
    ax.set_xlabel("r")
    ax.set_ylabel(r"$c_n$")
    ax.grid(alpha=0.3)
    ax.set_title("Nutriente: <0.4 não cresce, <0.6 não se move", fontsize=10)
    ax.legend(fontsize=8)

    ax = fig.add_subplot(gs[1, 1])
    for r in res:
        ax.plot(R_BINS, r["prof_cs"], "-", lw=1.6, label=r["run"])
    ax.axhline(0.95, c="crimson", ls="--", lw=1.2)
    ax.text(
        2.6, 0.965, "saturação: sem gradiente, sem motor", fontsize=8, color="crimson"
    )
    ax.axvspan(JUNC_LO, ARM_R0, color="0.85", zorder=0)
    ax.set_xlabel("r")
    ax.set_ylabel(r"$c_s/c_{s,max}$")
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05)
    ax.set_title("Surfactante: encher biomassa empurra a curva para 1", fontsize=10)
    ax.legend(fontsize=8)

    for j, r in enumerate(res):
        ax = fig.add_subplot(gs[2, j])
        d, bio = r["d"], r["bio"]
        ax.scatter(
            d["x"][d["is_filler"] > 0.5],
            d["y"][d["is_filler"] > 0.5],
            s=1,
            c="0.88",
            lw=0,
        )
        sc = ax.scatter(
            d["x"][bio],
            d["y"][bio],
            s=14,
            c=d["rho_b_grown"][bio],
            cmap="viridis",
            vmin=0,
            vmax=1,
            lw=0,
        )
        for rr in (JUNC_LO, ARM_R0):
            ax.add_patch(
                plt.Circle((0, 0), rr, fill=False, ls="--", ec="crimson", lw=1)
            )
        ax.set_aspect("equal")
        ax.set_xlim(-5, 5)
        ax.set_ylim(-5, 5)
        ax.set_title(
            f"{r['run']}  t={r['t']:.0f}  n_bio={r['n_bio']}  V={r['V']:.2f}",
            fontsize=10,
        )
        plt.colorbar(sc, ax=ax, fraction=0.046, label=r"$\rho_b$")

    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=130, bbox_inches="tight")
    print(f"figura: {out}")


if __name__ == "__main__":
    main()
