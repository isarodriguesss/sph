import numpy as np
import h5py
import glob
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def ft(f):
    sd = f["solver_data"]
    return float(sd.attrs["t"]) if "t" in sd.attrs else float(np.array(sd["t"]))


fig, axes = plt.subplots(2, 2, figsize=(12, 11.5))
for j, run in enumerate(("runs/J0", "runs/J7")):
    fs = sorted(glob.glob(f"{run}/main_output/*.hdf5"))
    fn = min(fs, key=lambda p: abs(ft(h5py.File(p, "r")) - 48))
    g = h5py.File(fn, "r")["particles"]["fluid"]["arrays"]
    rb = g["rho_b_grown"][:]
    x = g["x"][:]
    y = g["y"][:]
    m = (np.abs(x) < 5) & (np.abs(y) < 5)
    ax = axes[0, j]
    ax.scatter(x[m], y[m], s=3, c=rb[m], cmap="viridis", vmin=0, vmax=1, lw=0)
    ax.set_title(f"{run[-2:]} — escala LINEAR 0-1 (o que voce ve hoje)", fontsize=11)
    ax = axes[1, j]
    sc = ax.scatter(
        x[m],
        y[m],
        s=3,
        c=np.maximum(rb[m], 1e-9),
        cmap="viridis",
        norm=LogNorm(vmin=1e-6, vmax=1),
        lw=0,
    )
    frac0 = 100 * np.mean(rb[m] == 0)
    ax.set_title(
        f"{run[-2:]} — escala LOG 1e-6..1  ({frac0:.0f}% em zero exato)", fontsize=11
    )
    plt.colorbar(sc, ax=ax, fraction=0.046)
    for a in (axes[0, j], axes[1, j]):
        a.set_aspect("equal")
        a.set_xlim(-5, 5)
        a.set_ylim(-5, 5)
fig.suptitle(
    "O halo e visivel ou e a escala de cor? J0 (zeros exatos) vs J7 (k_col=0.03)",
    fontsize=13,
)
fig.tight_layout()
fig.savefig("runs/juncao_J7_escala.png", dpi=125, bbox_inches="tight")
print("runs/juncao_J7_escala.png")
