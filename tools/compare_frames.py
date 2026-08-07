"""Comparativo VISUAL entre rodadas isoladas (Protocolo CLAUDE.md §11).

Uso:  python tools/compare_frames.py runs/E runs/R0 runs/R1 ...
      python tools/compare_frames.py --t 50 runs/*

Renderiza o mesmo instante de cada rodada com escalas de cor FIXAS (o plot.py usa
percentil por frame, o que impediria comparacao entre rodadas) e ancora nas duas
referencias obrigatorias (assets/reference.jpg e assets/reference_result.png).

Colunas: swarm(rho) | biomassa vs inseridas | sigma_a | cs
"""

import sys
import os
import glob
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import h5py

DOMAIN = None  # inferido por rodada (rodadas antigas [-5,5], novas [-7,7])
MARKER = 15

# Escalas FIXAS — identicas em todas as rodadas (requisito da comparacao)
RHO_VMIN, RHO_VMAX = 0.97, 1.13  # cristas dos braços (estilo reference.jpg)
SIG_VMIN, SIG_VMAX = 0.60, 1.05
CS_VMAX = 0.40
SIG_TRIG = 0.85

OUT = "runs/comparativo.png"


def load_arrays(fn):
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]

        def g(n):
            return np.array(a[n]) if n in a and np.array(a[n]).size else None

        d = {
            k: g(k)
            for k in (
                "x",
                "y",
                "rho",
                "rho_b_grown",
                "cs",
                "sigma_a",
                "is_filler",
                "is_wake",
            )
        }
        t = None
        try:
            t = float(np.array(f["solver_data"]["t"]))
        except Exception:
            pass
        d["t"] = t
    return d


def pick_frame(run, t_target):
    fs = sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5")))
    if not fs:
        return None
    if t_target is None:
        return fs[-1]
    best, best_dt = fs[-1], 1e9
    for fn in fs:
        d = load_arrays(fn)
        if d["t"] is None:
            continue
        if abs(d["t"] - t_target) < best_dt:
            best, best_dt = fn, abs(d["t"] - t_target)
    return best


def inserted_mask(d):
    n = d["x"].size
    for k in ("is_filler", "is_wake"):
        v = d.get(k)
        if v is not None and v.size == n:
            return v > 0.5
    return np.zeros(n, bool)


def main(argv):
    t_target = None
    runs = []
    i = 0
    while i < len(argv):
        if argv[i] == "--t":
            t_target = float(argv[i + 1])
            i += 2
        else:
            runs.append(argv[i])
            i += 1
    runs = [r for r in runs if os.path.isdir(r)]
    if not runs:
        print("nenhuma rodada; uso: python tools/compare_frames.py runs/R0 runs/R1")
        return

    nrow = len(runs) + 1  # +1 para a linha das referencias
    fig, axes = plt.subplots(nrow, 4, figsize=(22, 5.4 * nrow))
    if nrow == 1:
        axes = axes[None, :]

    # --- linha 0: referencias obrigatorias ---
    for j, (path, title) in enumerate(
        [
            ("assets/reference.jpg", "REFERENCIA experimental (PA14, Michiels)"),
            ("assets/reference_result.png", "REFERENCIA numerica (Trinschek 2018)"),
        ]
    ):
        ax = axes[0, j]
        if os.path.exists(path):
            ax.imshow(mpimg.imread(path))
            ax.set_title(title, fontsize=11, fontweight="bold")
        ax.axis("off")
    for j in (2, 3):
        axes[0, j].axis("off")
    axes[0, 2].text(
        0.5,
        0.5,
        "escalas de cor FIXAS\nentre rodadas",
        ha="center",
        va="center",
        fontsize=12,
        style="italic",
    )

    for r, run in enumerate(runs, start=1):
        name = os.path.basename(run.rstrip("/"))
        fn = pick_frame(run, t_target)
        if fn is None:
            for j in range(4):
                axes[r, j].axis("off")
            axes[r, 0].text(0.5, 0.5, f"{name}: sem frames", ha="center")
            continue
        d = load_arrays(fn)
        x, y = d["x"], d["y"]
        ins = inserted_mask(d)
        bio = ~ins & (d["rho_b_grown"] > 0.1)
        tstr = f"t={d['t']:.1f}s" if d["t"] is not None else ""

        ax = axes[r, 0]
        ax.scatter(
            x,
            y,
            c=d["rho"],
            cmap="inferno",
            s=MARKER,
            vmin=RHO_VMIN,
            vmax=RHO_VMAX,
            linewidths=0,
        )
        ax.set_facecolor("black")
        ax.set_title(f"{name} — swarm (rho) {tstr}", fontweight="bold")

        ax = axes[r, 1]
        ax.scatter(x, y, c="#101020", s=2, linewidths=0)
        if ins.any():
            ax.scatter(
                x[ins],
                y[ins],
                c="#ff4d4d",
                s=MARKER,
                linewidths=0,
                label=f"inseridas ({int(ins.sum())})",
            )
        if bio.any():
            ax.scatter(
                x[bio],
                y[bio],
                c="#4dff88",
                s=MARKER,
                linewidths=0,
                label=f"biomassa ({int(bio.sum())})",
            )
        ax.set_facecolor("black")
        ax.legend(loc="upper right", fontsize=8, framealpha=0.6)
        ax.set_title(f"{name} — biomassa vs inseridas")

        ax = axes[r, 2]
        sa = d["sigma_a"]
        if sa is not None:
            colony = d["rho_b_grown"] > 0.1
            ax.scatter(x, y, c="#dddddd", s=1, linewidths=0)
            sc = ax.scatter(
                x[colony],
                y[colony],
                c=sa[colony],
                cmap="RdYlGn",
                s=MARKER,
                vmin=SIG_VMIN,
                vmax=SIG_VMAX,
                linewidths=0,
            )
            plt.colorbar(sc, ax=ax, fraction=0.046)
            low = colony & (sa < SIG_TRIG)
            ax.scatter(
                x[low],
                y[low],
                facecolors="none",
                edgecolors="blue",
                s=MARKER * 3,
                linewidths=0.5,
            )
            frac = low.sum() / max(colony.sum(), 1)
            ax.set_title(f"{name} — sigma_a (azul: <{SIG_TRIG}, {frac:.0%})")
        else:
            ax.axis("off")

        # cs: fillers quimicamente transparentes tem cs CONGELADO no valor herdado
        # — nao sao campo, sao historico. Renderiza-los pinta braços falsos.
        ax = axes[r, 3]
        chem = ~ins
        sc = ax.scatter(
            x[chem],
            y[chem],
            c=d["cs"][chem],
            cmap="hot",
            s=MARKER,
            vmin=0,
            vmax=CS_VMAX,
            linewidths=0,
        )
        plt.colorbar(sc, ax=ax, fraction=0.046)
        ax.set_facecolor("black")
        ax.set_title(f"{name} — cs (só partículas químicas)")

        lim = max(float(np.max(np.abs(x))), float(np.max(np.abs(y))))
        for j in range(4):
            axes[r, j].set_xlim(-lim, lim)
            axes[r, j].set_ylim(-lim, lim)
            axes[r, j].set_aspect("equal")

    plt.tight_layout()
    os.makedirs("runs", exist_ok=True)
    plt.savefig(OUT, dpi=100, bbox_inches="tight")
    print(f"comparativo visual salvo em {OUT}  ({len(runs)} rodadas)")


if __name__ == "__main__":
    main(sys.argv[1:])
