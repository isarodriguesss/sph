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
- `MarangoniForce`: `f = -β ∇cs · gate(|∇rho_b|)` — only active at the colony interface (smoothstep gate on `grad_rho_b_mag`, default `[0.05, 0.6]`). Uses the **symmetric SPH gradient** `∇cs_i ≈ (1/ρ_i)·Σ_j m_j·(cs_j − cs_i)·∇W_ij` (Pass E, ram-8827-v1, [equations.py:134](src/equations.py#L134)). The earlier non-symmetric form `(1/ρ_i)·Σ_j m_j·cs_j·∇W_ij` was found to capture only **transient kernel-deficit** at sharp boundaries and decay as the colony "wets" — see resolved failure mode below. The symmetric form requires a cs field with a transition layer **resolved across at least ~1.5·h** to produce a non-zero pairwise difference `cs_ij`.
- `SurfactantEquation`: `Dcs/Dt = σ·qs(rho_b)·(1 − rho_b)·noise + D·∇²cs − λ·cs`, where `qs = rho_b²/(rho_b² + K²)` is a **Hill function (quorum-sensing model)** with `K² = 0.01` ([equations.py:75-78](src/equations.py#L75-L78)). The `(1 − rho_b)` factor (Pass F applied) gates production OFF in the mature core (rho_b → 1) so the core becomes a cs SINK and the active growth ring (rho_b ≈ 0.5) is the production peak. This creates a moving production front rather than a static saturated source. **Confirmed working topologically** (Pass F log: `mean_cs` no longer saturates, decays to zero — but motor amplitude is now subcritical, see "production-front extinction" failure mode). The diffusive Brookshaw Laplacian uses the symmetric form `cs_ij = d_cs − s_cs` weighted by `(x_ij·∇W_ij)/|r_ij|²`.
- `LinearDrag`: `f = -γ_eff·v` with `γ_eff = γ_base + γ_mature·rho_b²` — biological mobility gradient: motile edge swarmers vs. immobile EPS-embedded core. **Currently `γ_mature = 0.3·γ_base`** in [scheme.py:125](src/scheme.py#L125) (was `2.0·γ_base` — caused core freezing).
- `ViscousForce`: standard viscous SPH term, with density-weighted `mu_eff = mu·min(rho_avg, 1)`
- `BiomassEOS`: "Soft Interior, Cohesive Edge" — quadratic repulsion when `ρ > ρ₀`, **slight negative pressure when `ρ < ρ₀`** (surface tension / cohesion via `MomentumEquation`). `edge_fade` smoothstep on `rho_b ∈ [0.1, 0.5]` ensures P=0 at the free colony border. **Default `tension_ratio=0.1`** (note: not overridden in [scheme.py:93](src/scheme.py#L93)) → attractive branch is 10% of repulsive strength.
- `OsmoticForce`: **disabled** (kept in file for reference). Replaced by EOS cohesive branch, which is more stable than `∇(rho_SPH)`-based pulling.

## Key parameters (tuned in `main.py`)

Current calibration (post Passes A+B+C+D+E.1+E.2+F+G, ram-8827-v1). Goal: dendritic fingering pattern (Michiels et al.). **MAJOR PROGRESS**: Pass G achieved a sustained motor for the first time — `mean_v` peak 0.019, `n_fast` peak 657, `a_mar` peak 49, active phase 27s, colony radius doubled from 0.65 → 1.7, **dendritic fingering morphology confirmed in frame 023**. Remaining issue (Pass H pending): late "tail extinction" at t≈37s — once the entire colony saturates uniformly with `rho_b → 1`, the `(1 − rho_b)` gate kills production globally and cs drains exponentially. Need a residual production floor to enable indefinite expansion. See "Known failure modes" below.

| Parameter | Variable | Current value | Notes |
|-----------|----------|:-------------:|-------|
| Viscosity | `mu` | 0.025 | Halved from 0.05 — allow small-scale perturbations to grow |
| Linear drag (base) | `gamma` | 60.0 | a_drag = γ·v_term = 6 at v=0.1 (was 120) |
| Drag (mature core) | `gamma_mature` (in scheme) | `0.3·γ` | Was `2.0·γ` — core was frozen at γ_eff=360 |
| Marangoni coeff | `beta` | 4.0 | Hold while testing Pass E; reduce to ~2.5 if a_mar > 30 sustained after D-cut |
| Surfactant production | `sigma` | 1.2 | Pass G applied. ×3 boost compensated the (1−rho_b) suppression and pushed motor above the bulk drag floor. Validated: a_mar peak 49, mean_v peak 0.019 |
| Diffusion | `D` | 4.0e-3 | L_D = 0.163 ≈ 1.5·h. Pass E.2 applied. Confirmed: produces correct spatial profile (a_mar peak = 16.28 vs predicted 17), but does NOT fix temporal saturation lock |
| Surfactant decay | `lambda_` | 0.15 | Hold — reducing it shortens L_D but proportionally lowers cs_eq, no net |∇cs| gain |
| Growth rate | `r_growth` | 0.8 | Pass G applied. Slow growth extends front residence τ_sat 0.8s → 2.5s, raises per-particle cs accumulation 11% → 31% of equilibrium. Validated: max_cs reached 3.18 |
| Production model | `qs · (1−rho_b)` | `[rho_b²/(rho_b²+0.01)] · (1−rho_b)` | Pass F applied. Pass H target: change to `(1.2 − rho_b)` to provide a 17% production floor in the mature core, preventing late extinction at t≈37s |
| Monaghan artificial visc | `alpha_mon` | 0.15 | Lowered to allow MS instability seeds to grow |
| Speed of sound (EOS) | `c0` | 0.8 | B ≈ 0.09; tension branch 10% (default `tension_ratio=0.1`) |
| Kernel smoothing | `h_factor` (in `particles.py`) | 1.8·dx | ~35 neighbors per particle |
| Timestep | `dt` in `create_solver` | 5e-5 (adaptive, CFL=0.4) | |
| Grid | `x_dim, y_dim` | 100×100 | |
| Domain | `x/y_min/max_domain` | [-3, 3]² | 6×6 centered at origin |
| Initial perturbation | `noise` in `main.py:77` | `1.0 + 0.25·sin(12θ) + 0.1·rand()` | Coherent N=12 azimuthal mode seed |

## Target acceleration budget (v_term = 0.1)

| Term | Target | Formula |
|------|:------:|---------|
| a_marangoni (líquida) | ~6 | β · \|∇cs\| · gate (only at tips/interface) |
| a_drag | ~6 | γ · v_term = 60·0.1 |
| a_pressão EOS | <3 | B · excess² · edge_fade (computed in main.py post_step from `au − ax_mar − ax_drag`) |
| a_viscosa | ~1.5 | μ · v / h² |
| **Total \|a\|** | **5–15** | Equilibrium: Marangoni ≈ Drag at tips |

## Known failure modes

- **Marangoni "thermal death" (t ≈ 30s) — RESOLVED**: After initial expansion, the colony spread into a low-density disk and `cs` production collapsed because it was linearly proportional to `rho_b`. **Fixed** by Hill quorum-sensing model in `SurfactantEquation.post_loop` (saturates production for `rho_b > 0.3`). Confirmed working: log shows `mean_cs` stable at ~0.057 indefinitely, no longer decaying.

- **Marangoni "diffusive death" (t ≈ 12s) — RESOLVED**: Was caused by `D = 3e-4` (too small). **Fixed** by raising `D` to 1.5e-2 and removing the `(0.3 + grad_rho_b)` production bias, so cs accumulates in the colony interior with smooth profile.

- **Marangoni "kernel asymmetry collapse" (t ≈ 11s) — RESOLVED (Pass E.1)**: The non-symmetric SPH gradient operator captured only transient kernel-deficit at sharp free-surface boundaries and decayed to zero once the boundary "wet" via EOS cohesion. **Fixed** by reverting to the symmetric form `cs_ij = s_cs - d_cs` in [equations.py:134](src/equations.py#L134). Confirmed: with the symmetric form, a_mar reached **14.6 sustained** for ~3s (iter 400, t=7.7s) — a real, well-formed signal — vs. the non-symmetric form's transient spike of 28→0.024.

- **Marangoni "diffusive flatness" (t ≈ 15s) — RESOLVED-PARTIAL (Pass E.2)**: With L_D = 2.9·h the spatial spreading was confirmed excessive. **Fixed** by reducing `D` from 1.5e-2 → 4.0e-3 (L_D = 1.5·h). Spatial signature confirmed: a_mar peak rose from 14.6 → **16.28** at iter 200 — quantitatively matching the prediction (~17). However, the underlying TEMPORAL saturation was not addressed: collapse still occurs at t≈15.8s (vs t=15.3s before) — only marginally better. The bottleneck is now reframed below.

- **Marangoni "production-saturation lock" (t ≈ 15s) — RESOLVED (Pass F)**: The Hill production saturated uniformly across the entire mature ring, creating a flat cs blanket with ∇cs → 0. **Fixed** by multiplying production by `(1 − rho_b/rho_max)` in [equations.py:77-78](src/equations.py#L77-L78). Confirmed: `mean_cs` no longer locks at σ/λ — it now decays to zero in the post-Pass-F log, proving the saturation attractor is gone. The mature core successfully became a sink.

- **CSV column drift — RESOLVED (Pass F instrument repair)**: writer at [main.py:185-200](main.py#L185-L200) now correctly emits all 13 columns including `mean_v` and `n_fast`. Confirmed: post-Pass-F log shows `mean_v` and `n_fast` populated, allowing bulk-vs-outlier disambiguation. **Critical finding enabled by repair**: in steady state, `n_fast = 1` confirms the persistent "single outlier particle" pollution that was masking diagnostics in all prior runs — `max_v = 0.222` and `a_drag = 17.3` in the final state are entirely due to ONE rogue particle, not the bulk.

- **Marangoni "production-front extinction" (t ≈ 17s) — RESOLVED (Pass G)**: Pass F's moving-front motor was too weak (τ_sat/τ_decay = 0.12) and too short to bootstrap expansion. **Fixed** by combined `sigma: 0.4 → 1.2` (×3 amplitude) + `r_growth: 2.5 → 0.8` (extends τ_sat from 0.8s → 2.5s, raises cs accumulation per particle from 11% → 31%). Validated: a_mar peak 49.0 (predicted 18, undersold by gate factor), max_cs peak 3.18 (predicted 2.5), mean_v peak 0.019 (predicted 0.03), n_fast peak **657** (predicted ≥ 50). The bootstrap paradox is closed.

- **Pass G dendritic morphology — ACHIEVED**: frame 023 of [main_output/movie/](main_output/movie/) (t≈30s) shows the colony at r≈1.7 with an irregular ameboid perimeter and ~30 small protrusions — clearly a Mullins-Sekerka fingering pattern in the noise-dominated regime. Not the clean N=12 organized dendrites of Michiels et al. (the azimuthal seed was overwhelmed by stochastic noise in the production term), but qualitatively a true dendritic colony for the first time. Active phase 3–37s (~34s of sustained motor), peak activity at iter 4600 (t=30s) with n_fast=657.

- **Marangoni "tail extinction" (t ≈ 37s) — CURRENT (Pass H)**: Pass G achieved sustained expansion but eventually hits the same wall as Pass F, just delayed: once the colony stops expanding (because it has filled the available domain or exhausted its biomass headroom), `rho_b → 1` everywhere, the `(1 − rho_b)` factor zeros production globally, and cs drains exponentially to zero with τ = 1/λ ≈ 6.7s. Confirmed in [log.csv](log.csv):
  - iter 4600 (t=30s): `mean_v = 0.019` ★ peak, `n_fast = 657`, `mean_cs = 0.226` (still climbing slightly)
  - iter 5400 (t=32.5s): `mean_v = 0.012` (declining), `mean_cs = 0.252` (stalled — saturation onset)
  - iter 6800 (t=37s): `mean_v = 0.0019`, `n_fast = 1` (back to outlier-only), `a_mar = 6.7` (declining)
  - iter 7400 (t=51s): `mean_cs = 0.030`, `max_cs = 0.40` (cs draining)
  - iter 9600 (t=100s): `mean_cs ≈ 0`, `a_mar ≈ 0`, locked in outlier-only steady state
  - `contrast_cs` collapses from peak 67 (iter 200) → 13 (iter 5400) → 10 (final), the same homogenization signature as Pass E.2.

  **Mechanistic difference from Pass F**: Pass F died at t=10s because the motor never even started (insufficient amplitude). Pass G dies at t=37s because the motor RAN OUT of fresh territory. The colony grew, expanded, but eventually filled the domain at the limit allowed by its finite biomass + the (1−rho_b) silencing.

  **Fix direction (Pass H)**: replace `growth_headroom = 1.0 - rho_b` with `growth_headroom = 1.2 - rho_b` in [equations.py:77](src/equations.py#L77). Single-character change. At rho_b=1.0 the factor is now 0.2 (was 0), maintaining 17% of peak production indefinitely in the mature core. Predicted: active phase preserved (factor at rho_b=0.5 only changes 0.5 → 0.7, a 40% boost in front production), late phase becomes a sustained "weak motor" steady state rather than extinction (cs_steady ≈ 1.6 in core, |∇cs| ≈ 10 at outer boundary, a_mar ≈ 15 indefinitely). Biophysical justification: rhamnolipid secretion is constitutively expressed in stationary biofilm at reduced rate (Lequette & Greenberg 2005), not zero. The binary `(1 − rho_b)` model is unphysical at rho_b → 1.
