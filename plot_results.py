import h5py
import matplotlib.pyplot as plt

import os
import glob

# Find the latest HDF5 file
files = glob.glob("main_output/main_*.hdf5")
if not files:
    print("No HDF5 files found in main_output/")
    exit(1)
filename = max(files, key=os.path.getctime)
print(f"Plotting latest file: {filename}")

with h5py.File(filename, "r") as f:
    # PySPH stores particles under 'particles'/'fluid'/'arrays'
    fluid = f["particles"]["fluid"]["arrays"]
    x = fluid["x"][:]
    y = fluid["y"][:]
    rho_b = fluid["rho_b_grown"][:]

    plt.figure(figsize=(8, 8), dpi=150)
    plt.scatter(x, y, c=rho_b, cmap="viridis", s=2, alpha=0.8)
    plt.colorbar(label="Biomass (rho_b)")
    plt.xlim(-3, 3)
    plt.ylim(-3, 3)
    plt.title("Pseudomonas Swarm Simulation (Iter: 4600)")
    plt.tight_layout()
    plt.savefig("render_latest.png")
    print("Saved render_latest.png successfully!")
