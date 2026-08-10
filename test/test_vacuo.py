import sys
import glob
import os
import numpy as np
from scipy.spatial import cKDTree
from pysph.solver.utils import load

RHO_B_COLONY = 0.05  # onde ha colonia
HOLE_THRESH = 0.7  # densidade SPH baixa no scatter
VOID_THRESH = 1.5  # em dx: distancia acima da qual a area conta como vazia
SIG_MIN = 0.85  # Violeau §3.6: abaixo disso ~15% de erro no operador
VOID_MAX = 0.01  # C1: no maximo 1% da area da colonia vazia
SIG_FRAC_MAX = 0.15  # C2: no maximo 15% da colonia abaixo de SIG_MIN
DEFAULT_RUN = "runs/C4"  # baseline validado (§2.5); sobrescreva pelo argumento


def find_snapshot(argv):
    """Ultimo HDF5 de uma rodada. Sem argumento, usa DEFAULT_RUN."""
    run = argv[1] if len(argv) > 1 else DEFAULT_RUN
    cands = sorted(glob.glob(os.path.join(run, "**", "main_*.hdf5"), recursive=True))
    if cands:
        return cands[-1]
    print(f"ERRO: nenhum main_*.hdf5 em {run}/")
    if len(argv) == 1:
        print(f"       ({DEFAULT_RUN} e o default; passe outra rodada como argumento)")
    return None


def void_fraction(x, y, rho_b, dx, thresholds=(0.7, 1.0, 1.5), n_grid=600):
    """Fracao da AREA da colonia sem nenhuma particula dentro de thr*dx.

    A grade e recortada pelo DOMINIO: quando a colonia passa da parede, o disco
    de raio R cobre regiao sem particula POR CONSTRUCAO, e isso apareceria como
    vacuo fisico. Esse bug ja produziu uma "explosao do vacuo" de 0.45% -> 21.6%
    que era 100% artefato (a fracao do disco fora do dominio em R=6.40 e 23.77%).
    """
    colony = rho_b > 0.1
    out = {t: float("nan") for t in thresholds}
    out["R"] = float("nan")
    if int(np.sum(colony)) < 10:
        return out
    r = np.hypot(x, y)
    R = float(np.percentile(r[colony], 99))
    xlim = float(np.max(np.abs(x)))
    ylim = float(np.max(np.abs(y)))
    g = np.linspace(-R, R, n_grid)
    GX, GY = np.meshgrid(g, g)
    ins = (
        (GX * GX + GY * GY <= R * R)
        & (np.abs(GX) <= xlim)
        & (np.abs(GY) <= ylim)
    )
    if not np.any(ins):
        return out
    d, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX[ins], GY[ins]])
    )
    res = {t: float(np.mean(d > t * dx)) for t in thresholds}
    res["R"] = R
    return res


def main(argv):
    hdf5 = find_snapshot(argv)
    if hdf5 is None:
        return 1
    f = load(hdf5)["arrays"]["fluid"]
    x, y, rho, rho_b = f.x, f.y, f.rho, f.rho_b_grown
    cs, m, sigma_a = f.cs, f.m, f.sigma_a
    v = np.sqrt(f.u**2 + f.v**2)
    h = float(np.median(f.h))
    dx = h / 1.8  # particles.py: h = 1.8*dx
    print(f"Snapshot: {hdf5}")
    print(f"  N={len(x)}   colonia={int((rho_b > RHO_B_COLONY).sum())}   dx={dx:.5f}\n")

    tree = cKDTree(np.column_stack([x, y]))
    colony = np.where(rho_b > RHO_B_COLONY)[0]
    holes = colony[rho[colony] < HOLE_THRESH]

    # ── (1) os "buracos" do scatter tem vizinhos? ────────────────────────
    nbr = (
        np.array([len(tree.query_ball_point([x[i], y[i]], 2.0 * h)) - 1 for i in holes])
        if len(holes)
        else np.array([0])
    )
    n_isolado = int((nbr < 5).sum())
    print("(1) Os 'buracos' do scatter sao queda de densidade, nao isolamento:")
    print(f"    holes (rho<{HOLE_THRESH}): {len(holes)} ({len(holes)/len(colony):.0%} da colonia)")
    print(f"    vizinhos em 2h: mediana={int(np.median(nbr))}  min={int(nbr.min())}")
    print(f"    quase-isoladas (<5 vizinhos): {n_isolado}")

    # ── (2) C2 — suporte do kernel ───────────────────────────────────────
    col = rho_b > 0.1
    sig_med = float(np.median(sigma_a[col]))
    sig_mean = float(np.mean(sigma_a[col]))
    frac_baixo = float((sigma_a[col] < SIG_MIN).mean())
    print("\n(2) C2 — suporte do kernel (Violeau §3.4/3.6):")
    print(f"    sigma_a na colonia: mediana={sig_med:.3f}  media={sig_mean:.3f}")
    print(f"    fracao com sigma_a < {SIG_MIN}: {frac_baixo:.1%}  (alvo <= {SIG_FRAC_MAX:.0%})")

    # ── (3) C1 — cobertura areal (o que sigma_a NAO ve) ──────────────────
    vf = void_fraction(x, y, rho_b, dx)
    area_dx2 = vf[VOID_THRESH] * np.pi * vf["R"] ** 2 / (dx * dx)
    print(f"\n(3) C1 — Fracao de Vazio AREAL (R_colonia={vf['R']:.2f}):")
    print("    fracao da AREA da colonia sem nenhuma particula dentro de:")
    for t in (0.7, 1.0, 1.5):
        print(f"      > {t:.1f} dx : {vf[t]:7.3%}")
    print(f"    vazio > {VOID_THRESH}dx em area absoluta: {area_dx2:.0f} dx^2"
          f"   (alvo <= {VOID_MAX:.0%})")
    print("    -> sigma_a nao enxerga isto: onde nao ha particula, nao ha amostra")

    # ── (4) fisica saudavel ──────────────────────────────────────────────
    # Filler e quimicamente transparente (S4): seu `cs` fica CONGELADO no valor
    # herdado e nao representa o campo. Incluí-lo infla mean_cs e deprime
    # contrast_cs artificialmente — licao #39.
    isf = getattr(f, "is_filler", None)
    chem = (isf < 0.5) if isf is not None and len(isf) == len(cs) else np.ones_like(cs, bool)
    cs_real = cs[chem]
    contrast_cs = (cs_real.max() - cs_real.min()) / (cs_real.mean() + 1e-9)
    contrast_all = (cs.max() - cs.min()) / (cs.mean() + 1e-9)
    print("\n(4) Fisica saudavel:")
    print(f"    contrast_cs (com filler, ENGANOSO) = {contrast_all:.1f}")
    print(f"    contrast_cs (so particulas quimicas) = {contrast_cs:.1f}")
    print(f"    mean_v = {float(np.mean(v)):.5f}")
    print(f"    massa = {float(np.sum(m)):.1f}")

    # ── Asseroes ─────────────────────────────────────────────────────────
    a1 = n_isolado == 0
    a2 = frac_baixo <= SIG_FRAC_MAX
    a3 = vf[VOID_THRESH] <= VOID_MAX
    a4 = contrast_cs > 10 and float(np.mean(v)) > 0
    print("\n" + "=" * 66)
    print(f"(1) nenhuma particula quase-isolada:              {'OK' if a1 else 'FALHA'}")
    print(f"C2  suporte do kernel (frac<{SIG_MIN} <= {SIG_FRAC_MAX:.0%}):        {'OK' if a2 else 'FALHA'}")
    print(f"C1  cobertura areal (vazio>{VOID_THRESH}dx <= {VOID_MAX:.0%}):        {'OK' if a3 else 'FALHA'}")
    print(f"(4) fisica saudavel (motor vivo):                 {'OK' if a4 else 'FALHA'}")
    ok = a1 and a2 and a3 and a4
    print("-" * 66)
    print(
        "✅ ESTADO DE PREENCHIMENTO DENSO — C1 e C2 satisfeitos simultaneamente"
        if ok
        else "❌ revisar: C1 e C2 sao INDISSOCIAVEIS (§2.5)"
    )
    print("=" * 66)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
