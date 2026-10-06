"""Anisotropia azimutal de R(theta): a rede de pilares imprime a propria simetria?

Rede triangular de rotacao 0 com o inoculo num sitio VAGO: a primeira coroa tem
pilar de frente em 0/60/120 graus e GARGANTA em 30/90/150. Se o lider trava na
geometria (transporte superamortecido em rede de obstaculos, enquadramento de
caminhada enviesada), R(theta) ganha componente m=6 com MAXIMOS em 30 graus.

Medido no R1 (secao 12 de docs/CRITERIOS_RUGOSIDADE.md): a_6 pico 0.223 contra
0.014 do liso no mesmo instante, com a FASE travada em 27-30 graus por nove
quadros, enquanto no liso ela passeia de 2 a 42. A fase e o que distingue
assinatura da rede de flutuacao do modelo — reportar sempre as duas.

A metrica primaria da serie (`dR/dt`) e media azimutal e e CEGA a isto: no R1 a
modulacao vale 22% e a media difere 3% do liso (licao #101 por outro caminho —
a regua escondia o maior efeito da rodada).

    python tools/anisotropia.py runs/swarm/P2R23_conduz runs/rugosidade/R1_lambda28
    python tools/anisotropia.py --autoteste

Colonia = `rho_b >= 0.1` ou filler (CLAUDE.md §2.2), R(theta) = maior raio com
colonia no bin, direto das PARTICULAS e nao do raster (licao #92).
"""

import argparse
import glob
import os

import h5py
import numpy as np

NTH = 720
M_MAX = 40


def frame(fn):
    with h5py.File(fn, "r") as f:
        a = f["particles"]["fluid"]["arrays"]
        d = {k: a[k][:] for k in ("x", "y", "rho_b_grown", "is_filler")}
        try:
            t = float(f["solver_data"].attrs["t"])
        except Exception:
            t = float(f.attrs["solver_data"]["t"])
    return d, t


def r_theta(d, nth=NTH):
    """R(theta) com indice 0 em theta=0 (eixo +x), passo 2*pi/nth."""
    col = (d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)
    r = np.hypot(d["x"][col], d["y"][col])
    th = np.arctan2(d["y"][col], d["x"][col])
    ib = ((th + np.pi) / (2.0 * np.pi) * nth).astype(int) % nth
    R = np.zeros(nth)
    np.maximum.at(R, ib, r)
    vaz = R == 0.0
    if vaz.any() and not vaz.all():
        i = np.arange(nth)
        R[vaz] = np.interp(i[vaz], i[~vaz], R[~vaz], period=nth)
    return np.roll(R, -nth // 2)


def modo(R, m):
    """Amplitude relativa a R medio e azimute do MAXIMO, em graus mod 360/m."""
    F = np.fft.rfft(R)
    amp = 2.0 * np.abs(F[m]) / np.abs(F[0])
    th = np.degrees(-np.angle(F[m]) / m) % (360.0 / m)
    return amp, th


def autoteste():
    """Calibra a regua antes de ranquear por ela (licoes #67-F, #101)."""
    th = np.arange(NTH) * 2.0 * np.pi / NTH
    ok = True
    for m, a, t0 in ((6, 0.20, 30.0), (6, 0.05, 12.0), (4, 0.10, 15.0), (12, 0.08, 7.0)):
        R = 2.0 * (1.0 + a * np.cos(m * (th - np.radians(t0))))
        am, tm = modo(R, m)
        bom = abs(am - a) < 1e-9 and abs(tm - t0) < 1e-6
        ok &= bom
        print(f"  m={m:2d} a={a:.2f} fase={t0:5.1f} -> {am:.6f} {tm:8.4f}  "
              f"{'ok' if bom else 'FALHOU'}")
    print("autoteste:", "ok" if ok else "FALHOU")
    return ok


def main():
    p = argparse.ArgumentParser()
    p.add_argument("runs", nargs="*")
    p.add_argument("--t-min", type=float, default=8.0)
    p.add_argument("--janela", type=float, default=35.0,
                   help="inicio da janela de media (CLAUDE.md §2.2)")
    p.add_argument("--autoteste", action="store_true")
    a = p.parse_args()

    if a.autoteste:
        raise SystemExit(0 if autoteste() else 1)

    for run in a.runs:
        fns = sorted(glob.glob(os.path.join(run, "main_output", "*.hdf5")))
        print(f"\n=== {os.path.basename(run)} ===")
        print(f"{'t':>6} {'R_med':>6} {'a6':>7} {'fase6':>6} | pico m>=4")
        acc = []
        for fn in fns:
            d, t = frame(fn)
            if t < a.t_min:
                continue
            R = r_theta(d)
            a6, th6 = modo(R, 6)
            F = np.fft.rfft(R)
            A = 2.0 * np.abs(F) / np.abs(F[0])
            top = np.argsort(A[4:M_MAX])[::-1][:3] + 4
            s = "  ".join(f"m{m}={A[m]:.3f}" for m in top)
            print(f"{t:6.1f} {R.mean():6.3f} {a6:7.4f} {th6:6.1f} | {s}")
            if t >= a.janela:
                acc.append((a6, th6))
        if acc:
            M = np.array(acc)
            print(f"  media t>={a.janela:.0f}: a6={M[:, 0].mean():.4f}+-{M[:, 0].std():.4f}"
                  f"  fase6={M[:, 1].mean():.1f}+-{M[:, 1].std():.1f}")


if __name__ == "__main__":
    main()
