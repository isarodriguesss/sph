"""Pares de particulas a menos de 0.05 dx, dt de saida e rho/rho0 p99 por frame (t >= 60).

    python tools/pares_dt.py runs/A runs/B ...

Instrumento da licao #89: so compare runs contados por ESTE script.
"""

import glob
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

DX = 14.0 / 260
for run in sys.argv[1:]:
    print(f"== {run}")
    print("   t      dt_saida   pares<0.05dx  rho/rho0 p99")
    for f in sorted(glob.glob(f"{run}/main_output/main_*.hdf5")):
        with h5py.File(f, "r") as h:
            sd = h["solver_data"].attrs
            t, dt = float(sd["t"]), float(sd["dt"])
            if t < 60:
                continue
            a = h["particles"]["fluid"]["arrays"]
            x, y, rho = (np.asarray(a[k]) for k in ("x", "y", "rho"))
        n = len(cKDTree(np.c_[x, y]).query_pairs(0.05 * DX))
        print(f"  {t:6.2f}  {dt:.5f}   {n:6d}        {np.percentile(rho, 99):.2f}")
