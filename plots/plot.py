import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysph.solver.utils import load
import glob
import os

TARGET_FRAMES = 40  # nao gera frame para TODO snapshot — subamostra ~40 + o ultimo
rho0 = 1.0  # densidade de referência SPH
DOMAIN = (-7.0, 7.0)  # janela = dominio cheio (expandido 2026-08-06)
MARKER = 15  # tamanho do ponto: 1-particula-de-largura le como braço conectado

# --- Footprint temporal (Parte 2): acumula a biomassa ocupada ao longo do run ---
# A reference.jpg e ela mesma uma foto do RASTRO integrado do swarm. O esqueleto
# instantaneo (~300 particulas) e esparso; acumular max(rho_b) ao longo dos frames
# preenche os braços de forma continua = pegada dendritica cheia.
FOOT_N = 300  # celulas da grade do footprint sobre DOMAIN
FOOT_SIGMA = 3.0  # suavizacao (celulas): preenche o tubo fino do braço
FOOT_RHO_TRIG = (
    1.05  # so as cristas densas dos braços (agar/comprimido<1.05, braço~1.1)
)
_fedges = np.linspace(DOMAIN[0], DOMAIN[1], FOOT_N + 1)


def deposit_footprint(x, y, rho, rho_b):
    """Ocupacao suavizada de braços+nucleo neste frame, para acumular via max no tempo.

    Braços: cristas de rho (rho>trig) — CONTINUAS a cada frame (a biomassa e esqueleto
    esparso demais, ~2-3 particulas/braço, p/ preencher por acumulacao). Nucleo: biomassa
    densa (rho_b>0.5), que tem rho abaixo da crista mas deve aparecer coeso como na
    reference.jpg. Suavizado, preenche o tubo do braço; o max no tempo traça o dendrito
    varrido completo — pegada dendritica cheia estilo reference.jpg.
    """
    from scipy.ndimage import gaussian_filter

    mask = (rho > FOOT_RHO_TRIG) | (rho_b > 0.5)
    hist, _, _ = np.histogram2d(y[mask], x[mask], bins=[_fedges, _fedges])
    return gaussian_filter((hist > 0).astype(float), sigma=FOOT_SIGMA)


# --- Reconstrucao de campo SPH (Shepard / σ-normalizado) para os paineis fisicos ---
# Price 2007 (splash, PASA 24, 159): campo continuo A(r) = Σ_j V_j A_j W / Σ_j V_j W.
# A normalizacao pelo denominador (Shepard) fecha os buracos do scatter — onde ha
# poucas particulas, reescala em vez de mostrar vazio. Como e viz (nao solver), os
# vacuos dos dendritos somem na figura sem tocar na fisica (licao #38).
GRID_N = 300  # resolucao da grade de reconstrucao sobre [-3, 3]²
KNN = 30  # vizinhos por celula (suporte 2·h_recon ≈ 2·dx → ~12 vizinhos; folga p/
# regioes densas com filler inserido, onde o espacamento cai abaixo de dx)
# Comprimento de suavizacao da reconstrucao, em unidades de dx (espacamento das
# particulas). h_recon ≈ dx preserva os braços finos (1-2 dx de largura) e ainda
# preenche gaps sub-dx; usar o h fisico (1.8·dx) borraria os dendritos em blobs.
H_RECON_FACTOR = 1.0

_gx, _gy = np.meshgrid(np.linspace(-3, 3, GRID_N), np.linspace(-3, 3, GRID_N))
_grid_pts = np.column_stack([_gx.ravel(), _gy.ravel()])


def _cubic_spline_shape(q):
    # Forma radial do cubic spline de Monaghan (sem a constante de normalizacao —
    # ela cancela na razao Shepard). Suporte compacto q = r/h < 2.
    w = np.zeros_like(q)
    m1 = q < 1.0
    m2 = (q >= 1.0) & (q < 2.0)
    w[m1] = 1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3
    w[m2] = 0.25 * (2.0 - q[m2]) ** 3
    return w


def reconstruct_field(px, py, values, vol, h):
    """Campo Shepard-normalizado dos `values` (das particulas) na grade regular."""
    from scipy.spatial import cKDTree

    tree = cKDTree(np.column_stack([px, py]))
    dist, idx = tree.query(_grid_pts, k=KNN)
    w = _cubic_spline_shape(dist / h) * vol[idx]
    num = np.sum(w * values[idx], axis=1)
    den = np.sum(w, axis=1)
    field = num / (den + 1e-12)
    return field.reshape(_gx.shape)


def render_all():
    os.makedirs("main_output/movie", exist_ok=True)
    all_files = sorted(glob.glob("main_output/main_*.hdf5"))
    stride = max(1, len(all_files) // TARGET_FRAMES)
    files = all_files[::stride]
    if all_files and all_files[-1] not in files:
        files.append(all_files[-1])
    footprint = np.zeros((FOOT_N, FOOT_N))  # acumulador do rastro de biomassa (Parte 2)
    for i, fpath in enumerate(files):
        _render_frame(i, fpath, footprint)
    print(f"Done. {len(files)} total files.")


def _render_frame(i, fpath, footprint):
    data = load(fpath)
    fluid = data["arrays"]["fluid"]
    x, y = fluid.x, fluid.y
    cs = fluid.cs
    rho = (
        fluid.rho
    )  # empacotamento SPH: revela a estrutura dendritica (cristas de braço)
    rho_b = fluid.rho_b_grown  # biomassa (nucleo denso p/ o footprint)

    # tempo do snapshot para o titulo (se disponivel)
    t = data.get("solver_data", {}).get("t", None)
    tstr = f"t={t:.1f}s" if t is not None else os.path.basename(fpath)

    # Footprint temporal: acumula o max da estrutura de braço (rho) ate este frame
    footprint[:] = np.maximum(footprint, deposit_footprint(x, y, rho, rho_b))

    # O swarm dendritico e revelado por rho (empacotamento SPH): agar ~0.85, braços
    # ~1.1 (mais densos), canais entre braços ~0.5. rho_b (biomassa) e um esqueleto
    # esparso (~300 particulas), ruim para campo — por isso renderizamos por PONTOS
    # grandes no dominio CHEIO, como o PySPH viewer faz.
    fig, axes = plt.subplots(1, 4, figsize=(32, 8))

    # --- Painel 1: rho — cientifico (fiel ao viewer) ---
    sc1 = axes[0].scatter(
        x, y, c=rho, cmap="viridis", s=MARKER, vmin=0.4, vmax=1.15, linewidths=0
    )
    axes[0].set_title(f"rho (empacotamento SPH) — {tstr}")
    plt.colorbar(sc1, ax=axes[0])

    # --- Painel 2: rho — estilo reference.jpg (biologico) ---
    # Fundo escuro (agar), dendritos claros. inferno com vmin no nivel do agar
    # mapeia o agar para ~preto e os braços densos para laranja/amarelo.
    axes[1].set_facecolor("black")
    # vmin ACIMA do agar (~0.85): agar e voids -> preto, so as cristas densas dos
    # braços (rho~1.0-1.15) acendem -> dendritos claros sobre fundo escuro.
    axes[1].scatter(
        x, y, c=rho, cmap="inferno", s=MARKER, vmin=0.97, vmax=1.13, linewidths=0
    )
    axes[1].set_title(f"swarm (estilo reference.jpg) — {tstr}")
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    # --- Painel 3: footprint temporal ∫rho_b (Parte 2) — braços preenchidos ---
    # Rastro acumulado da biomassa: preenche o esqueleto esparso num padrao
    # dendritico continuo, casando com a pegada integrada da reference.jpg.
    axes[2].set_facecolor("black")
    fp_pos = footprint[footprint > 0]
    fp_vmax = max(float(np.percentile(fp_pos, 98)), 1e-6) if fp_pos.size else 1.0
    axes[2].imshow(
        footprint,
        origin="lower",
        extent=[DOMAIN[0], DOMAIN[1], DOMAIN[0], DOMAIN[1]],
        cmap="inferno",
        vmin=0,
        vmax=fp_vmax,
        aspect="equal",
    )
    axes[2].set_title(f"footprint ∫rho_b (estilo reference.jpg) — {tstr}")
    axes[2].set_xticks([])
    axes[2].set_yticks([])

    # --- Painel 4: cs (surfactante) — halo de Marangoni alem da biomassa ---
    cs_vmax = max(float(np.percentile(cs, 99.5)), 0.1)
    sc4 = axes[3].scatter(
        x, y, c=cs, cmap="hot", s=MARKER, vmin=0, vmax=cs_vmax, linewidths=0
    )
    axes[3].set_title(f"cs (surfactante) — {tstr}")
    plt.colorbar(sc4, ax=axes[3])

    for a in (axes[0], axes[1], axes[3]):
        a.set_xlim(DOMAIN)
        a.set_ylim(DOMAIN)
        a.set_aspect("equal")

    plt.tight_layout()
    out = f"main_output/movie/frame_i7_{i:03d}.png"
    plt.savefig(out, dpi=120, bbox_inches="tight")
    plt.close()
    print(f"Frame {i:03d} — {tstr}")


if __name__ == "__main__":
    render_all()
