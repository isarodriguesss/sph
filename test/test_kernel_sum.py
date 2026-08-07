"""
Teste-ouro: KernelSum (σ_a = Σ_j (m_j/ρ_j)·W_ij) — a particao da unidade discreta.

Reproduz EXATAMENTE a matematica de KernelSum (src/equations.py) com o mesmo
kernel CubicSpline(dim=2) do solver (main.py), e mede σ_a em dois cenarios de
referencia — os dois que o comentario no proprio KernelSum pede:

  σ_a ≈ 1  → SUPORTE COMPLETO: particula no interior de uma rede uniforme dx.
             Vale a particao da unidade (Violeau §3.6, Liu §3.3.3): a soma dos
             volumes ponderados pelo kernel reconstroi a unidade. Operadores SPH
             consistentes (~0% erro de ordem zero).

  σ_a ≈ 0  → VACUO: um ponto de agar sem particula nenhuma dentro do raio 2h.
             Nenhum vizinho → a soma e zero. E a assinatura numerica da "regiao
             que nao parece ter particulas".

Cenarios intermediarios (contexto, nao criterio de pass/fail):
  - particula ISOLADA real (so ela mesma, conta o auto-termo W(0)): σ_a = piso
    de auto-suporte ≈ 0.14. E o MINIMO que uma particula real (que se conta)
    atinge — NAO chega a 0. Documentado para nao confundir "isolada" com "vacuo".
  - SUPERFICIE LIVRE (particula SOBRE uma borda plana, meio-plano de vizinhos,
    kernel truncado, Liu §6.5): σ_a ≈ 0.7 — mantem o auto-termo + a metade
    inferior do suporte (NAO cai a 0.5 porque a propria particula se conta).

Parametros identicos ao solver: dx = 10/186, h = 1.8·dx, m = dx², ρ0 = 1.0.

── PARTE 2 (dados reais) ────────────────────────────────────────────────────
Alem dos cenarios sinteticos, mede σ_a em um HDF5 REAL da simulacao (o ultimo
frame de main_output/). Confirma o caso do cliente: o MIOLO de um braco
dendritico e vazio de verdade — σ_a ≈ 0.03 (0 particulas dentro de h; so as
paredes, a ~2h, onde o peso do kernel e minusculo) — enquanto uma baia de agar
preenchida tem σ_a alto. Cruza tambem a reimplementacao contra o `sigma_a` que
o proprio solver (KernelSum) gravou no arquivo.
"""

import sys
import glob
import numpy as np
from pysph.base.kernels import CubicSpline

DX = 10.0 / 186.0  # main.py: (x_max-x_min)/(x_dim-1) = 10/186
H = 1.8 * DX  # particles.py: h = 1.8·dx
RHO0 = 1.0  # densidade de referencia (rede uniforme)
M = DX * DX  # particles.py: m = dx·dx
VOL = M / RHO0  # volume por particula = m/ρ
KERNEL = CubicSpline(dim=2)
SUPPORT = KERNEL.radius_scale * H  # raio de suporte do kernel = 2h


def sigma_a(px, py, xs, ys):
    """Replica KernelSum.loop: σ_a = Σ_j (m_j/ρ_j)·W(r_aj, h).

    xs, ys: posicoes dos vizinhos (INCLUINDO a propria particula, como no
    solver — o loop PySPH percorre self quando dest e source sao o mesmo array).
    """
    s = 0.0
    for xj, yj in zip(xs, ys):
        dxv, dyv = px - xj, py - yj
        r = (dxv * dxv + dyv * dyv) ** 0.5
        w = KERNEL.kernel([dxv, dyv, 0.0], r, H)
        s += VOL * w
    return s


def uniform_lattice(cx, cy, half_extent):
    """Rede uniforme dx cobrindo [c-half, c+half]² (folga > 2h para suporte cheio)."""
    n = int(np.ceil(half_extent / DX))
    offs = np.arange(-n, n + 1) * DX
    gx, gy = np.meshgrid(cx + offs, cy + offs)
    return gx.ravel(), gy.ravel()


def sigma_a_real(px, py, xs, ys, ms, rhos, hq):
    """σ_a = Σ_j (m_j/ρ_j)·W num ponto de consulta sobre PARTICULAS reais.

    Usa m_j/ρ_j de cada particula (nao o VOL uniforme) — fiel a KernelSum.loop.
    hq = h de consulta (mediana das particulas; e o h do campo reconstruido).
    """
    s = 0.0
    for xj, yj, mj, rj in zip(xs, ys, ms, rhos):
        dxv, dyv = px - xj, py - yj
        r = (dxv * dxv + dyv * dyv) ** 0.5
        s += (mj / rj) * KERNEL.kernel([dxv, dyv, 0.0], r, hq)
    return s


def real_data_section():
    """Mede σ_a num HDF5 real: vazio do braco (~0.03) vs agar/nucleo preenchido.

    Retorna True (passou), False (falhou) ou None (sem HDF5 → skip, nao falha).
    """
    from pysph.solver.utils import load
    from scipy.spatial import cKDTree

    files = sorted(glob.glob("main_output/main_*.hdf5"))
    if not files:
        print("\n[PARTE 2] sem HDF5 em main_output/ — pulando (rode `make run`).")
        return None

    f = load(files[-1])["arrays"]["fluid"]
    x, y, rho, rho_b, m = f.x, f.y, f.rho, f.rho_b_grown, f.m
    sig_solver = f.sigma_a  # σ_a que o proprio KernelSum gravou
    hq = float(np.median(f.h))
    tree = cKDTree(np.column_stack([x, y]))
    print(
        f"\n[PARTE 2] dados reais: {files[-1].split('/')[-1]}  N={len(x)}  h={hq:.4f}"
    )

    def sig_at(px, py):
        idx = tree.query_ball_point([px, py], SUPPORT)  # vizinhos em 2h
        if not idx:
            return 0.0, 0
        s = sigma_a_real(px, py, x[idx], y[idx], m[idx], rho[idx], hq)
        return s, len(idx)

    # ── (A) cross-check: reimplementacao ≈ σ_a do solver (particula do bulk) ──
    # particula bem no interior do agar (muitos vizinhos) → σ_a ≈ 1 nos dois
    counts = tree.query_ball_point(np.column_stack([x, y]), SUPPORT, return_length=True)
    ib = int(np.argmax(counts))  # a particula com MAIS vizinhos = bulk denso
    sig_mine, _ = sig_at(x[ib], y[ib])
    sig_solv = float(sig_solver[ib])
    print(
        f"  (A) cross-check bulk: minha σ_a={sig_mine:.3f}  solver σ_a={sig_solv:.3f}"
    )

    # ── (B) MIOLO DO BRACO: maior vazio dentro da colonia → σ_a ≈ 0.03 ──────
    g = np.linspace(-4.5, 4.5, 300)
    gx, gy = np.meshgrid(g, g)
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    d1, _ = tree.query(pts, k=1)  # dist ao vizinho mais proximo = tamanho do vazio
    r = np.sqrt(pts[:, 0] ** 2 + pts[:, 1] ** 2)
    cand = np.where((d1 > 2 * DX) & (r < 4))[0]  # vazios largos, dentro da colonia
    if len(cand):
        k = cand[np.argmax(d1[cand])]
        pv = pts[k]
        n_in_h = len(tree.query_ball_point(pv, hq))  # particulas dentro de h
        sig_void, n_in_2h = sig_at(pv[0], pv[1])
        print(
            f"  (B) miolo do braco ({pv[0]:.2f},{pv[1]:.2f}): "
            f"{n_in_h} part em h, {n_in_2h} em 2h → σ_a={sig_void:.3f}"
        )
    else:
        sig_void, n_in_h = 1.0, -1
        print("  (B) nenhum vazio largo encontrado (run com inserção na frontier?)")

    # ── (C) AGAR/NUCLEO preenchido: σ_a alto (tem dado) ─────────────────────
    icore = int(np.argmax(rho_b))  # centro do nucleo, denso
    sig_fill, _ = sig_at(x[icore], y[icore])
    print(
        f"  (C) nucleo preenchido ({x[icore]:.2f},{y[icore]:.2f}): σ_a={sig_fill:.3f}"
    )

    # ── Asseroes ─────────────────────────────────────────────────────────
    a_xcheck = abs(sig_mine - sig_solv) < 0.05  # reimplementacao == solver
    a_void = (n_in_h == 0) and (sig_void < 0.15)  # braco = vazio real
    a_fill = sig_fill > 0.5  # nucleo = tem suporte
    print(f"  → cross-check reimplementacao==solver:  {'OK' if a_xcheck else 'FALHA'}")
    print(f"  → miolo do braco e vazio (σ_a<0.15):    {'OK' if a_void else 'FALHA'}")
    print(f"  → nucleo preenchido (σ_a>0.5):          {'OK' if a_fill else 'FALHA'}")
    return a_xcheck and a_void and a_fill


def main():
    pad = 1.5 * SUPPORT  # folga generosa para nao truncar o kernel no interior

    # ── Cenario σ_a ≈ 1: suporte completo ────────────────────────────────
    lx, ly = uniform_lattice(0.0, 0.0, pad)
    sig_full = sigma_a(0.0, 0.0, lx, ly)

    # ── Cenario σ_a ≈ 0: vacuo (probe sem nenhum vizinho em 2h) ──────────
    # rede deslocada para longe: o probe em (0,0) nao tem particula em 2h
    fx, fy = uniform_lattice(10.0, 10.0, pad)
    sig_void = sigma_a(0.0, 0.0, fx, fy)

    # ── Contexto: particula isolada real (so o auto-termo) ───────────────
    sig_isolated = sigma_a(0.0, 0.0, [0.0], [0.0])

    # ── Contexto: superficie livre (meio-plano y<=0) ─────────────────────
    hx, hy = uniform_lattice(0.0, 0.0, pad)
    keep = hy <= 1e-12
    sig_edge = sigma_a(0.0, 0.0, hx[keep], hy[keep])

    print(
        f"dx={DX:.5f}  h={H:.5f}  suporte(2h)={SUPPORT:.5f}  W(0)={KERNEL.kernel([0, 0, 0], 0, H):.3f}\n"
    )
    print(f"σ_a SUPORTE COMPLETO (rede uniforme)  = {sig_full:.4f}   (alvo ≈ 1)")
    print(f"σ_a VACUO (nenhum vizinho em 2h)      = {sig_void:.4f}   (alvo = 0)")
    print(
        f"σ_a isolada (auto-termo W(0) apenas)  = {sig_isolated:.4f}   (piso, NAO chega a 0)"
    )
    print(
        f"σ_a superficie livre (borda plana)    = {sig_edge:.4f}   (~0.7, kernel truncado)"
    )

    # ── Asseroes ─────────────────────────────────────────────────────────
    a_one = abs(sig_full - 1.0) < 0.02  # particao da unidade satisfeita
    a_zero = sig_void == 0.0  # vacuo exato
    a_floor = 0.10 < sig_isolated < 0.20  # piso de auto-suporte, longe de 1
    a_edge = 0.60 < sig_edge < 0.80  # auto-termo + metade inferior

    print("\n" + "=" * 60)
    print(f"σ_a ≈ 1 no suporte completo (|σ-1|<0.02):  {'OK' if a_one else 'FALHA'}")
    print(f"σ_a = 0 no vacuo:                          {'OK' if a_zero else 'FALHA'}")
    print(f"isolada no piso ~0.14 (contexto):          {'OK' if a_floor else 'FALHA'}")
    print(f"superficie livre ~0.7 (contexto):          {'OK' if a_edge else 'FALHA'}")
    ok_synth = a_one and a_zero and a_floor and a_edge
    print("-" * 60)
    print(
        "✅ [PARTE 1] KernelSum: 1 no bulk, 0 no vacuo — particao da unidade correta"
        if ok_synth
        else "❌ [PARTE 1] revisar"
    )

    # ── PARTE 2: dados reais (skip se nao houver HDF5) ────────────────────
    real = real_data_section()

    print("\n" + "=" * 60)
    if real is None:
        print("✅ PARTE 1 OK (PARTE 2 pulada — sem HDF5)" if ok_synth else "❌ revisar")
        return 0 if ok_synth else 1
    ok = ok_synth and real
    print(
        "✅ σ_a distingue VAZIO REAL do braco (~0.03) de suporte pleno (~1)"
        if ok
        else "❌ revisar"
    )
    print("=" * 60)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
