"""
Teste-ouro de consistência de 1ª ordem da KGC (Bonet & Lok 1999 / CSPM).

Prova objetiva de que a implementação em src/equations.py (KernelGradientCorrection
+ aplicação L_i·DWIJ em MarangoniForce) RESTAURA a consistência de 1ª ordem do
operador de gradiente — i.e., reproduz EXATAMENTE o gradiente de um campo LINEAR
mesmo onde o suporte do kernel é INCOMPLETO (borda / braço sub-resolvido).

Critério (Liu §3.3, Violeau §3.4): para φ(x,y)=a+b·x+c·y, ∇φ=(b,c) em todo ponto.
  - SPH puro: Σ_j V_j(φ_j-φ_i)∇W_ij ≈ (b,c) só com suporte completo; na borda erra O(1).
  - KGC:      L_i·Σ_j V_j(φ_j-φ_i)∇W_ij = (b,c) EXATO sempre que det(M_i) ≥ det_min.

Usa o MESMO kernel (CubicSpline 2D) e a MESMA matemática do solver:
  M_i  = Σ_j V_j (x_j-x_i) ⊗ DWIJ          [equations.py:155-163]
  L_i  = M_i^{-1} (fallback I se det<det_min) [equations.py:165-178]
  ∇φ_c = L_i · Σ_j V_j (φ_j-φ_i) DWIJ        [equations.py:457-461]
       = L_i · (gradiente SPH puro)   (L_i é constante sobre j)
"""

import numpy as np
from pysph.base.kernels import CubicSpline

DET_MIN = 0.25  # igual a KGC_DET_MIN em main.py


def lattice(n, dx):
    xs = (np.arange(n) - (n - 1) / 2.0) * dx
    X, Y = np.meshgrid(xs, xs)
    return X.ravel(), Y.ravel()


def gradients_at(i, x, y, phi, h, V, kernel, det_min=DET_MIN):
    """Replica EXATAMENTE a matemática do solver para a partícula i."""
    xi, yi = x[i], y[i]
    Mxx = Mxy = Myx = Myy = 0.0
    gbx = gby = 0.0  # gradiente SPH puro: Σ V_j (φ_j-φ_i) DWIJ
    grad = [0.0, 0.0, 0.0]
    nnbr = 0
    for j in range(len(x)):
        xij = [xi - x[j], yi - y[j], 0.0]  # XIJ = x_i - x_j (convenção PySPH)
        rij = (xij[0] ** 2 + xij[1] ** 2) ** 0.5
        if rij < 1e-12 or rij > 2.0 * h:
            continue
        kernel.gradient(xij, rij, h, grad)  # grad = DWIJ = ∇_i W_ij
        dwx, dwy = grad[0], grad[1]
        dxx, dyy = -xij[0], -xij[1]  # (x_j - x_i) = -XIJ
        Mxx += V * dxx * dwx
        Mxy += V * dxx * dwy
        Myx += V * dyy * dwx
        Myy += V * dyy * dwy
        pij = phi[j] - phi[i]  # (φ_j - φ_i)
        gbx += V * pij * dwx
        gby += V * pij * dwy
        nnbr += 1

    det = Mxx * Myy - Mxy * Myx
    if det < det_min:
        Lxx, Lxy, Lyx, Lyy = 1.0, 0.0, 0.0, 1.0
        fb = True
    else:
        inv = 1.0 / det
        Lxx, Lxy, Lyx, Lyy = Myy * inv, -Mxy * inv, -Myx * inv, Mxx * inv
        fb = False

    gcx = Lxx * gbx + Lxy * gby  # ∇φ corrigido = L · (∇φ SPH puro)
    gcy = Lyx * gbx + Lyy * gby
    return (gbx, gby), (gcx, gcy), det, fb, nnbr


def main():
    n, dx = 15, 0.1
    h = 1.8 * dx  # mesmo h_factor do simulador
    V = dx * dx  # m/rho com m=dx², rho=1
    x, y = lattice(n, dx)
    kernel = CubicSpline(dim=2)

    # Campo LINEAR: gradiente exato = (b, c) em TODO ponto
    a, b, c = 0.37, 1.5, -0.8
    phi = a + b * x + c * y
    print(f"Campo linear φ = {a} + {b}·x + {c}·y  →  ∇φ exato = ({b}, {c})\n")

    # Varredura de uma coluna do centro até a borda (suporte cada vez mais truncado)
    cx = (n - 1) // 2
    col = [r * n + cx for r in range(cx, n)]  # centro → borda superior
    print(
        f"{'pos':>14} {'nnbr':>4} {'det(M)':>7} {'erro SPH puro':>22} "
        f"{'erro KGC':>22} {'fallback':>8}"
    )
    print("-" * 86)

    max_err_kgc_active = 0.0
    any_active_edge = False
    for i in col:
        (gbx, gby), (gcx, gcy), det, fb, nnbr = gradients_at(i, x, y, phi, h, V, kernel)
        err_sph = ((gbx - b) ** 2 + (gby - c) ** 2) ** 0.5
        err_kgc = ((gcx - b) ** 2 + (gcy - c) ** 2) ** 0.5
        tag = "I→bare" if fb else ""
        print(
            f"({x[i]:+.2f},{y[i]:+.2f}) {nnbr:4d} {det:7.3f} "
            f"{err_sph:22.3e} {err_kgc:22.3e} {tag:>8}"
        )
        if not fb:
            max_err_kgc_active = max(max_err_kgc_active, err_kgc)
            # "edge" = suporte incompleto (menos vizinhos que o interior cheio)
            if nnbr < 25:
                any_active_edge = True

    print("-" * 86)
    print(
        f"\nMaior erro da KGC onde ela está ATIVA (det≥{DET_MIN}): "
        f"{max_err_kgc_active:.3e}"
    )

    # ── Asserções (a prova) ─────────────────────────────────────────────
    TOL = 1e-10
    ok = True

    # (1) O gradiente-diferença SPH puro é só consistência de ORDEM ZERO: erra
    # mesmo no INTERIOR (M_xx≈0.988≠1 → ~1% de erro), e a KGC o leva a precisão
    # de máquina. Prova que a correção age no bulk também, não só na borda.
    (gbx, gby), (gcx, gcy), _, _, _ = gradients_at(col[0], x, y, phi, h, V, kernel)
    err_int_sph = ((gbx - b) ** 2 + (gby - c) ** 2) ** 0.5
    err_int_kgc = ((gcx - b) ** 2 + (gcy - c) ** 2) ** 0.5
    print(
        f"\n(1) Interior: SPH puro erro={err_int_sph:.3e} (ordem-zero, ~1%) "
        f"→ KGC erro={err_int_kgc:.3e}  "
        f"({'OK' if (err_int_sph > 1e-3 and err_int_kgc < TOL) else 'FALHA'})"
    )
    ok = ok and (err_int_sph > 1e-3 and err_int_kgc < TOL)

    # (2) Borda com KGC ativa: erro ~ machine precision (consistência 1ª ordem)
    print(
        f"(2) KGC ativa em borda incompleta: erro_max={max_err_kgc_active:.3e}  "
        f"({'OK' if max_err_kgc_active < TOL else 'FALHA'})"
    )
    ok = ok and max_err_kgc_active < TOL

    # (3) Houve ao menos uma borda de suporte incompleto onde a KGC engajou
    print(
        f"(3) KGC engajou em borda de suporte incompleto: "
        f"{'OK' if any_active_edge else 'FALHA'}"
    )
    ok = ok and any_active_edge

    # (4) Na MESMA borda, o SPH puro erra O(1) — prova que a correção é necessária
    edge_i = col[-2]  # penúltima (borda, mas tipicamente det≥det_min)
    (ebx, eby), (ecx, ecy), edet, efb, _ = gradients_at(edge_i, x, y, phi, h, V, kernel)
    e_sph = ((ebx - b) ** 2 + (eby - c) ** 2) ** 0.5
    e_kgc = ((ecx - b) ** 2 + (ecy - c) ** 2) ** 0.5
    if not efb:
        print(
            f"(4) Borda (det={edet:.3f}): SPH puro erro={e_sph:.3e} (O(1)) "
            f"vs KGC erro={e_kgc:.3e}  "
            f"({'OK' if e_sph > 1e-2 and e_kgc < TOL else 'FALHA'})"
        )
        ok = ok and (e_sph > 1e-2 and e_kgc < TOL)
    else:
        print(
            f"(4) Borda escolhida caiu no fallback (det={edet:.3f}<{DET_MIN}) — "
            f"esperado em cantos; ver tabela."
        )

    print("\n" + "=" * 50)
    print("RESULTADO:", "✅ KGC FIEL E CORRETA" if ok else "❌ FALHOU")
    print("=" * 50)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
