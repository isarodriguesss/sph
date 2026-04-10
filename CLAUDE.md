# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A 2D Smoothed Particle Hydrodynamics (SPH) simulation of bacterial biomass dynamics using [PySPH](https://pysph.readthedocs.io/). It models biomass growth, Marangoni surface forces driven by surfactant gradients, viscous flow, and linear drag in a fluid domain.

## Commands

```bash
make run          # Clear output dir and run simulation (python main.py)
make view         # Visualize output with PySPH viewer
make run_view     # Run then view
make format       # Lint and format with Ruff
make paraview     # Export HDF5 output to VTK for ParaView
```

The simulation outputs HDF5 files to `main_output/` at intervals defined by `print_freq`.

## Environment setup

```bash
conda create -n pysph_env python=3.10 numpy scipy matplotlib -c conda-forge
conda activate pysph_env
conda install pysph cython mako -c conda-forge
pip install PySPH
conda install mayavi -c conda-forge
conda install mpi4py -c conda-forge
```

## Architecture

**Entry point:** `main.py` — defines `SwarmApp(Application)`, sets all physical parameters as module-level globals, creates particles, scheme, and solver. The `post_step` hook prints per-step diagnostics and optionally handles cell division (`use_splitting = False` by default).

**`src/particles.py`** — `create_initial_state()` builds the initial 2D grid: fluid particles with a Gaussian `rho_b_grown` (biomass) field plus perturbations, and 2-layer solid ghost particles forming the boundary walls.

**`src/scheme.py`** — `MyBiomassScheme(Scheme)` wires all equations into two PySPH `Group`s:
1. Pre-step (non-real): `SummationDensity` → `BiomassEOS`
2. Main step: `MomentumEquation` (with strong Monaghan `alpha`) → `BiomassGrowth` → `BiomassGradient` → `MarangoniForce` → `ViscousForce` → `LinearDrag` → `SurfactantEquation`

`BiomassGradient` **must** run before `MarangoniForce` — Marangoni uses `grad_rho_b_mag` as an interface gate.

`CustomEulerStep` extends `EulerStep` to also integrate `rho_b_grown` (biomass density) and `cs` (surfactant concentration), with clamping: `rho_b_grown` ∈ [0, 1], `cs` ≥ 1e-9.

**`src/equations.py`** — All custom SPH equations as `Equation` subclasses:
- `BiomassGrowth`: logistic growth `drho_b/dt = r * rho_b * (1 - rho_b/rho_max)`
- `BiomassGradient`: SPH gradient of biomass field; magnitude used by Marangoni gate
- `MarangoniForce`: `f = -β ∇cs · gate(|∇rho_b|)` — only active at the colony interface (smoothstep gate on `grad_rho_b_mag ∈ [0.3, 1.5]`), preventing spurious bulk-interior pushes from cs noise
- `SurfactantEquation`: `Dcs/Dt = σ·rho_b·(0.3 + grad_rho_b)·noise + D∇²cs - λ·cs`
- `LinearDrag`: `f = -γ·v` (background damping — controls terminal velocity)
- `ViscousForce`: standard viscous SPH term, with density-weighted `mu_eff = mu·min(rho_avg, 1)`
- `BiomassEOS`: "Soft Interior, Cohesive Edge" — quadratic repulsion when `ρ > ρ₀`, **slight negative pressure when `ρ < ρ₀`** (surface tension / cohesion via `MomentumEquation`). `edge_fade` smoothstep on `rho_b ∈ [0.1, 0.5]` ensures P=0 at the free colony border. `tension_ratio=0.3` → attractive branch is 30% of repulsive strength.
- `OsmoticForce`: **disabled** (kept in file for reference). Replaced by EOS cohesive branch, which is more stable than `∇(rho_SPH)`-based pulling.

## Key parameters (tuned in `main.py`)

| Parameter | Variable | Current value | Notes |
|-----------|----------|:-------------:|-------|
| Viscosity | `mu` | 0.05 | Medium — maintains dendrite arm continuity |
| Linear drag | `gamma` | 120.0 | a_drag = γ·v_term = 12 at v=0.1 |
| Marangoni coeff | `beta` | 1.5 | Gated by interface — only fires where \|∇rho_b\| > 0.1 |
| Surfactant production | `sigma` | 1.5 | cs_eq = σ/λ = 30 |
| Diffusion | `D` | 1e-3 | Elevated → allows thicker and longer fingers |
| Surfactant decay | `lambda_` | 0.05 | Slower decay → more range for instability |
| Growth rate | `r_growth` | 1.0 | Slow → finer dendrites |
| Monaghan artificial visc | `alpha_mon` | 0.5 | Lowered to allow MS instability |
| Speed of sound (EOS) | `c0` | 0.8 | B ≈ 0.14; tension branch 10% |
| Kernel smoothing | `h_factor` (in `particles.py`) | 1.8·dx | ~35 neighbors per particle |
| Timestep | `dt` in `create_solver` | 5e-5 (adaptive, CFL=0.4) | |
| Grid | `x_dim, y_dim` | 100×100 | |
| Domain | `x/y_min/max_domain` | [-3, 3]² | 6×6 centered at origin |

## Target acceleration budget (v_term = 0.1)

| Term | Target | Formula |
|------|:------:|---------|
| a_marangoni (líquida) | ~10 | β · \|∇cs\| · gate (only at tips/interface) |
| a_drag | ~12 | γ · v_term = 120·0.1 |
| a_pressão EOS (repulsão) | <2 | B · excess² · edge_fade |
| a_pressão EOS (coesão) | <1 | −B·0.3 · deficit · edge_fade (attractive) |
| a_viscosa | ~3 | μ · v / h² |
| **Total \|a\|** | **10–30** | Equilibrium: Marangoni ≈ Drag at tips |
