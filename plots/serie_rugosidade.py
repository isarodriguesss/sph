"""Figuras de comparacao da serie de rugosidade (Pass L): liso contra todas as rugosas.

    python plots/serie_rugosidade.py                       # serie padrao
    python plots/serie_rugosidade.py liso=runs/swarm/P2R23_conduz \
        R1pp=runs/rugosidade/R1pp_lambda28_anel_shift R2=runs/rugosidade/R2_lambda20
    python plots/serie_rugosidade.py --refaz               # ignora o cache

Pos-processamento puro dos quadros HDF5 e do log.csv de cada rodada (nao roda simulacao).
A primeira rodada da lista e o LISO; as outras precisam de main_output/pilares.json. Uma
rodada nova entra so acrescentando `nome=caminho`; as escalas sao fixas.

Cada figura responde a uma pergunta da tese; tres seguem figuras do artigo do professor
(Cesgrav 2026): crescimento e area x raio (Fig. 6), painel por porosidade com as particulas
(Fig. 3) e grandeza x porosidade em rede e em meio sorteado (Fig. 7 dir.).

  fig_serie_crescimento  area contra t; area contra raio de giro (log-log)
  fig_serie_dose         raio, expansao e motor relativos ao liso, e D_b, contra a porosidade
  fig_serie_lideres      numero de lideres e motor contra t
  fig_serie_quimica      sombra de c_s atras do pilar, halo L/R
  fig_serie_painel       colonias em t=50 lado a lado, particulas sobre o nutriente (Fig. 3)

DEFINICOES:
  colonia      rho_b >= 0.1 ou filler (§2.2); R99 = percentil 99 do raio dela
  area         area da mascara da fig_tese (regiao sem agar) com os pilares descontados
  raio de giro raiz da media de r^2 das particulas da colonia
  lider        viva (nao filler, rho_b >= 0.1) com r >= 0.7*R99 (tools/desfecho_contato.py) e
               avancando (u_r >= 0.02); sem o filtro de velocidade, recrutas recem-convertidas
               junto aos pilares inflam a contagem em picos transitorios
  D_b          contagem de caixas da mascara da colonia, caixas de 1 a 16 dx
  controle     o liso e medido contra a geometria de CADA rodada rugosa (pilares virtuais):
               e a hipotese nula da sombra quimica
  ruido        par de realizacoes do P2 (P2_fillerdonor x P2_t100, §20.7)
"""

# pass-l-aprovado: pos-processamento da serie de rugosidade (desbloqueada em 2026-09-17); le
# a geometria gravada em pilares.json e nao cria contorno nenhum.

import argparse
import csv
import glob
import json
import os
import pickle
import sys

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
sys.path.insert(0, os.path.join(AQUI, "..", "tools"))
from fig_tese import CMAP_CS, campo, spline  # noqa: E402
from lit_rank import campo_stats, dedos, rtheta  # noqa: E402

DX = 14.0 / 260.0
L_JANELA = 6.0
N_GRADE = 444  # celula = 0.5 dx
JANELA = (35.0, 50.0)
JANELA_R = (35.0, 48.0)
JANELA_DRDT = (25.0, 48.0)
CAIXAS_DX = (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0)
VERSAO = 10
N_PM, W_PM = 49, 4.0  # grade do pilar medio: 49 x 49 pontos em +-4 raios
UR_ATIVO = 0.02  # lider que avanca; recruta recem-convertida anda 0.003-0.01
CACHE = os.path.join(AQUI, "_cache_serie")

SERIE = [
    ("liso", "runs/swarm/P2R23_conduz"),
    ("R1″", "runs/rugosidade/R1pp_lambda28_anel_shift"),
    ("R2", "runs/rugosidade/R2_lambda20"),
]
RUIDO = ("runs/swarm/P2_fillerdonor", "runs/swarm/P2_t100")

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "font.size": 8.5,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 7,
    }
)


# --------------------------------------------------------------------------- leitura


def quadros(run):
    out = []
    for f in sorted(glob.glob(os.path.join(run, "main_output", "main_*.hdf5"))):
        with h5py.File(f, "r") as h:
            out.append((float(h["solver_data"].attrs["t"]), f))
    return out


def geometria(run):
    p = os.path.join(run, "main_output", "pilares.json")
    if not os.path.exists(p):
        return None
    j = json.load(open(p))
    c = np.asarray(j["centros"], float)
    if not len(c):
        return None
    sorteado = j.get("tipo") == "booleano" or (j.get("jitter") or 0) > 0
    return dict(
        c=c,
        R=float(j["raio"]),
        Rf=float(j.get("fluido_raio") or j["raio"]),  # onde o fluido comeca
        rot=float(j.get("rotacao_graus") or 0.0),
        tipo="sorteado" if sorteado else "rede",
        coroa=float(np.min(np.hypot(*c.T))),
        phi=porosidade(c, float(j["raio"])),
    )


def porosidade(c, R, r0=1.2, r1=5.5):
    g = np.arange(-r1, r1, 0.25 * DX)
    X, Y = np.meshgrid(g, g)
    rr = np.hypot(X, Y)
    anel = (rr > r0) & (rr < r1)
    d, _ = cKDTree(c).query(np.c_[X[anel], Y[anel]])
    return float(np.mean(d > R))


def carrega(f):
    with h5py.File(f, "r") as h:
        a = h["particles"]["fluid"]["arrays"]
        d = {
            k: np.asarray(a[k], float)
            for k in (
                "x",
                "y",
                "u",
                "v",
                "m",
                "rho",
                "rho_b_grown",
                "is_filler",
                "cs",
                "c_n",
            )
        }
    return d


def le_log(run):
    p = os.path.join(run, "log.csv")
    if not os.path.exists(p):
        return {}
    with open(p) as fp:
        rows = list(csv.DictReader(fp))
    out = {}
    for k in rows[0]:
        try:
            out[k] = np.array([float(r[k]) for r in rows])
        except (TypeError, ValueError):
            pass
    return out


# --------------------------------------------------------------------------- medidas


def box_counting(mask, cel):
    out = []
    for s_dx in CAIXAS_DX:
        k = max(1, int(round(s_dx * DX / cel)))
        n = mask.shape[0] // k
        b = mask[: n * k, : n * k].reshape(n, k, n, k).any(axis=(1, 3))
        out.append((k * cel, int(b.sum())))
    s, N = np.array(out, float).T
    D = -np.polyfit(np.log(s), np.log(np.maximum(N, 1)), 1)[0]
    return D, s, N


def shepard_pts(d, pts, h=1.8 * DX, k=32, den=False):
    dist, idx = cKDTree(np.c_[d["x"], d["y"]]).query(pts, k=k)
    w = spline(dist / h) * (d["m"] / d["rho"])[idx]
    v = (w * d["cs"][idx]).sum(1) / np.maximum(w.sum(1), 1e-12)
    return (v, w.sum(1)) if den else v


def pilar_medio(d, geo, r99, ids=None, furar=True):
    """Campo de c_s em volta de cada pilar do corpo, no referencial do pilar.

    Eixo u = direcao radial (u > 0 aponta para FORA da colonia), v = tangencial, os dois em
    unidades do raio do pilar. Pilares do corpo: centro alem de r = 1.2 e borda antes de
    0.8 R99 da PROPRIA rodada; `ids` forca outra lista (o liso usa a do rugoso). Ponto
    dentro do pilar ou sem particula no suporte (sigma < 0.5) vira NaN.
    """
    c, R = geo["c"], geo["R"]
    rc = np.hypot(*c.T)
    if ids is None:
        ids = np.where((rc > 1.2) & (rc + R < 0.8 * r99))[0]
    g = np.linspace(-W_PM, W_PM, N_PM)
    U, V = np.meshgrid(g, g)
    # entre a superficie do pilar (R) e o inicio do fluido (Rf) nao ha particula de fluido:
    # ali a interpolacao so ve um lado e sai enviesada, entao o furo vai ate Rf
    fora = np.hypot(U, V) > geo["Rf"] / R if furar else np.ones_like(U, bool)
    mapas = np.full((len(ids), N_PM, N_PM), np.nan, np.float32)
    for k, i in enumerate(ids):
        er = c[i] / rc[i]
        et = np.array([-er[1], er[0]])
        pts = c[i] + R * (U[..., None] * er + V[..., None] * et)
        v, den = shepard_pts(d, pts.reshape(-1, 2), den=True)
        sigma = den * 10.0 / (7.0 * np.pi * (1.8 * DX) ** 2)  # spline 2D normalizado
        ok = fora & (sigma.reshape(U.shape) > 0.5)
        mapas[k] = np.where(ok, v.reshape(U.shape), np.nan)
    return np.asarray(ids), mapas


def sombra(d, geo, r99):
    """c_s logo a FRENTE do pilar (lado do centro) sobre c_s logo ATRAS, nos pilares do corpo."""
    c, R = geo["c"], geo["R"]
    rc = np.hypot(*c.T)
    sel = (rc > 1.2) & (rc + R < 0.8 * r99)
    if sel.sum() < 2:
        return np.nan
    u = c[sel] / rc[sel, None]
    off = R + 2.0 * DX
    frente = shepard_pts(d, c[sel] - off * u)
    tras = shepard_pts(d, c[sel] + off * u)
    return float(np.median(frente) / np.median(tras))


def mede_quadro(f, t, geo_proprio, geos_ctrl, liso):
    d = carrega(f)
    col = (d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)
    viva = (d["rho_b_grown"] >= 0.1) & (d["is_filler"] < 0.5)
    r = np.hypot(d["x"], d["y"])
    r99 = float(np.percentile(r[col], 99))
    out = dict(
        r99=r99, rg=float(np.sqrt(np.mean(r[col] ** 2))), n_col=int(col.sum())
    )

    # mascara da colonia na grade fixa (fig_tese.campo, pilares perfurados)
    d["pilar"] = None if geo_proprio is None else (geo_proprio["c"], geo_proprio["R"])
    g, img, _, _ = campo(d, L_JANELA, N_GRADE, False, agar_iso=0.0)
    cel = g[1] - g[0]
    corpo = img > 0.5
    out["area"] = float(corpo.sum() * cel * cel)
    out["Db"], out["box_s"], out["box_N"] = box_counting(corpo, cel)
    Rt, centro, _, _ = rtheta(corpo, so_maior=False)
    out["dedos"] = dedos(Rt)

    # halo de c_s (lit_rank.campo_stats, todas as particulas: o filler conduz no P2R23)
    GX, GY = np.meshgrid(g, g)
    cs_g = shepard_pts(d, np.c_[GX.ravel(), GY.ravel()]).reshape(corpo.shape)
    out["halo_L"] = campo_stats(cs_g, corpo, centro, Rt.max())["L_R"]

    # lideres
    lid = viva & (r >= 0.7 * r99)
    ur = (d["x"] * d["u"] + d["y"] * d["v"]) / np.maximum(r, 1e-12)
    out["n_lid"] = int(lid.sum())
    out["n_lid_ativo"] = int(np.sum(lid & (ur >= UR_ATIVO)))

    # sombra quimica: geometria propria e, no liso, as geometrias de controle
    geos = geos_ctrl if liso else {"proprio": geo_proprio}
    out["sombra"] = {nome: sombra(d, geo, r99) for nome, geo in geos.items()}
    if t >= JANELA[0] - 1e-6:
        # no liso, todos os pilares (virtuais) com centro alem de r=1.2; a selecao do
        # rugoso e aplicada depois, quadro a quadro
        out["pilar_medio"] = {
            nome: pilar_medio(
                d,
                geo,
                r99,
                None if not liso else np.where(np.hypot(*geo["c"].T) > 1.2)[0],
                furar=not liso,  # no liso nao ha pilar: o campo existe dentro do circulo
            )
            for nome, geo in geos.items()
        }
    return out


def mede_run(nome, run, geos_ctrl, refaz):
    os.makedirs(CACHE, exist_ok=True)
    qs = [(t, f) for t, f in quadros(run) if t <= JANELA[1] + 0.5]
    chave = (VERSAO, run, [os.path.getmtime(f) for _, f in qs], sorted(geos_ctrl))
    arq = os.path.join(CACHE, os.path.basename(run.rstrip("/")) + ".pkl")
    if not refaz and os.path.exists(arq):
        dados = pickle.load(open(arq, "rb"))
        if dados["chave"] == chave:
            dados["nome"] = nome
            return dados
    geo = geometria(run)
    print(f"  medindo {nome} ({len(qs)} quadros)...", flush=True)
    fr = []
    for t, f in qs:
        q = mede_quadro(f, t, geo, geos_ctrl, geo is None)
        q["t"] = t
        fr.append(q)
        print(
            f"    t={t:5.1f}  R99={q['r99']:.3f}  Db={q['Db']:.2f}  lideres={q['n_lid']}",
            flush=True,
        )
    dados = dict(chave=chave, nome=nome, run=run, geo=geo, fr=fr, log=le_log(run))
    pickle.dump(dados, open(arq, "wb"))
    return dados


def r99_serie(run):
    t, r = [], []
    for tt, f in quadros(run):
        if tt > JANELA[1] + 0.5:
            break
        with h5py.File(f, "r") as h:
            a = h["particles"]["fluid"]["arrays"]
            col = (np.asarray(a["rho_b_grown"]) >= 0.1) | (
                np.asarray(a["is_filler"]) > 0.5
            )
            rr = np.hypot(np.asarray(a["x"])[col], np.asarray(a["y"])[col])
        t.append(tt)
        r.append(np.percentile(rr, 99))
    return np.array(t), np.array(r)


def janela_r(t, r):
    """R99 medio em JANELA_R e dR/dt ajustado em JANELA_DRDT, na grade comum de 0.5."""
    tg = np.arange(JANELA_R[0], JANELA_R[1] + 1e-9, 0.5)
    tg2 = np.arange(JANELA_DRDT[0], JANELA_DRDT[1] + 1e-9, 0.5)
    return float(np.mean(np.interp(tg, t, r))), float(
        np.polyfit(tg2, np.interp(tg2, t, r), 1)[0]
    )


def ruido():
    arq = os.path.join(CACHE, "ruido.pkl")
    if os.path.exists(arq):
        return pickle.load(open(arq, "rb"))
    vals = [janela_r(*r99_serie(run)) for run in RUIDO if os.path.isdir(run)]
    if len(vals) < 2:
        out = (0.011, 0.058)  # §20.7
    else:
        (a1, b1), (a2, b2) = vals
        out = (abs(a1 - a2) / (0.5 * (a1 + a2)), abs(b1 - b2) / (0.5 * (b1 + b2)))
    os.makedirs(CACHE, exist_ok=True)
    pickle.dump(out, open(arq, "wb"))
    return out


# --------------------------------------------------------------------------- utilidades


def serie(D, k):
    return np.array([q["t"] for q in D["fr"]]), np.array(
        [q[k] for q in D["fr"]], float
    )


def na_janela(D, k, jan=JANELA):
    t, v = serie(D, k)
    m = (t >= jan[0] - 1e-6) & (t <= jan[1] + 0.5)
    return float(np.nanmean(v[m])), float(np.nanstd(v[m]))


def estilo(D, cores):
    if D["geo"] is None:
        return dict(color="k", ls="-", marker="o")
    rede = D["geo"]["tipo"] == "rede"
    return dict(
        color=cores[D["nome"]], ls="-" if rede else "--", marker="o" if rede else "s"
    )


def phi(D):
    return 1.0 if D["geo"] is None else D["geo"]["phi"]


def salva(fig, out, nome):
    for e in ("pdf", "png"):
        fig.savefig(os.path.join(out, f"{nome}.{e}"), dpi=250)
    plt.close(fig)
    print(f"  {out}/{nome}.png")






# --------------------------------------------------------------------------- figuras


def fig_crescimento(DS, cores, out):
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1), constrained_layout=True)
    a = ax[0]
    for D in DS:
        t, v = serie(D, "area")
        a.plot(t, v, **estilo(D, cores), ms=2.5, label=D["nome"])
    a.set_xlabel(r"$t$")
    a.set_ylabel("área da colônia")
    a.set_title("(a) área contra o tempo", loc="left")
    a.legend()

    a = ax[1]
    for D in DS:
        _, rg = serie(D, "rg")
        _, ar = serie(D, "area")
        a.loglog(rg, ar, **estilo(D, cores), ms=2.5, label=D["nome"])
    x = np.array([0.4, 3.0])
    a.loglog(x, 2.0 * x**2, ":", color="0.5", lw=0.8)
    a.text(1.9, 2.0 * 1.9**2 * 0.55, r"$\propto R_g^2$" + "\n(disco cheio)", color="0.4", fontsize=7, va="top")
    a.set_xlabel(r"raio de giro $R_g$")
    a.set_ylabel("área da colônia")
    a.set_title("(b) área contra o tamanho", loc="left")
    salva(fig, out, "fig_serie_crescimento")


def fig_dose(DS, cores, rui, out):
    Rl, Vl = janela_r(*serie(DS[0], "r99"))
    p95 = {}
    for D in DS:
        lg = D["log"]
        if "a_mar_bio_p95" in lg:
            p95[D["nome"]] = float(np.median(lg["a_mar_bio_p95"][lg["t"] > 30]))
    ml = p95.get(DS[0]["nome"], np.nan)
    # razao rodada/liso (o liso vale 1 por construcao): forma usual de dose-resposta em
    # meio poroso, como v/v0 ou D_eff/D0 contra a porosidade
    paineis = [
        (
            r"(a) raio $R_{99}/R_{99}^{\rm liso}$",
            lambda D: janela_r(*serie(D, "r99"))[0] / Rl,
            rui[0],
        ),
        (
            r"(b) expansão $\dot R/\dot R^{\rm liso}$",
            lambda D: janela_r(*serie(D, "r99"))[1] / Vl,
            rui[1],
        ),
        (
            r"(c) motor $a_{\rm mar}^{p95}/a_{\rm mar}^{p95,\rm liso}$",
            lambda D: p95.get(D["nome"], np.nan) / ml,
            None,
        ),
        (r"(d) $D_b$ médio, t ∈ [35, 50]", lambda D: na_janela(D, "Db")[0], None),
    ]
    erros = [None, None, None, lambda D: na_janela(D, "Db")[1]]
    fig, ax = plt.subplots(1, 4, figsize=(9.8, 2.8), constrained_layout=True)
    for a, (tit, f, banda), ferr in zip(ax, paineis, erros):
        if banda is not None:
            a.axhspan(1 - banda, 1 + banda, color="0.85", lw=0, label="ruído (P2 × P2)")
        for tipo, ls, mk in (("rede", "-", "o"), ("sorteado", "--", "s")):
            sel = [D for D in DS if D["geo"] is None or D["geo"]["tipo"] == tipo]
            if tipo == "sorteado" and len(sel) < 2:
                continue
            xs = np.array([phi(D) for D in sel])
            ys = np.array([f(D) for D in sel])
            o = np.argsort(xs)
            a.plot(xs[o], ys[o], ls, color="0.4", lw=0.8, zorder=1)
            for D, x, y in zip(sel, xs, ys):
                a.errorbar(
                    x,
                    y,
                    ferr(D) if ferr else None,
                    fmt=mk,
                    color=estilo(D, cores)["color"],
                    ms=5,
                    capsize=2,
                    zorder=2,
                    mfc="none" if tipo == "sorteado" else None,
                )
                a.annotate(
                    D["nome"],
                    (x, y),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=6.5,
                )
        a.set_xlabel(r"porosidade $\phi$")
        a.set_title(tit, loc="left", fontsize=8)
        a.invert_xaxis()
    for a in ax[:3]:
        lo, hi = a.get_ylim()
        a.set_ylim(min(lo, 0.78), max(hi, 1.03))
        a.axhline(1, color="0.6", lw=0.6, zorder=0)
    ax[0].legend(loc="lower left")
    salva(fig, out, "fig_serie_dose")


def fig_lideres(DS, cores, out):
    fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8), constrained_layout=True)
    for D in DS:
        st = estilo(D, cores)
        t, v = serie(D, "n_lid_ativo")
        ax[0].plot(t[t > 0], v[t > 0], **st, ms=2.5, label=D["nome"])  # em t=0 u=0
        t, v = serie(D, "n_lid")
        ax[0].plot(t, v, color=st["color"], ls=":", lw=0.9, alpha=0.6)
        lg = D["log"]
        if "a_mar_bio_p95" in lg:
            m = lg["t"] <= JANELA[1] + 0.5
            v = ndi.uniform_filter1d(lg["a_mar_bio_p95"][m], 5)
            ax[1].plot(lg["t"][m], v, color=st["color"], ls=st["ls"], label=D["nome"])
    ax[0].plot([], [], ":", color="0.4", lw=0.9, label="todas as vivas na frente")
    ax[0].set_title(r"(a) líderes que avançam ($r \geq 0.7R_{99}$, $u_r \geq 0.02$)", loc="left", fontsize=8)
    ax[1].set_title("(b) motor a_mar_bio_p95 (média móvel)", loc="left")
    for a in ax:
        a.set_xlabel(r"$t$")
    ax[0].legend(loc="upper right")
    salva(fig, out, "fig_serie_lideres")


def pares_pilar(D, liso):
    """Mapas normalizados dos MESMOS pilares no rugoso e no liso, quadro a quadro.

    Para cada quadro do rugoso, toma o quadro do liso mais proximo no tempo e os pilares que o
    rugoso selecionou. Cada mapa e dividido pela media de c_s na faixa do lado do centro
    (u em [-4, -3], |v| <= 1): o que sobra e o efeito do pilar sobre o que chega do centro.
    """
    M, Ml = [], []
    fl = [q for q in liso["fr"] if "pilar_medio" in q]
    for q in D["fr"]:
        if "pilar_medio" not in q:
            continue
        ids, mapas = q["pilar_medio"]["proprio"]
        ql = min(fl, key=lambda x: abs(x["t"] - q["t"]))
        ids_l, mapas_l = ql["pilar_medio"][D["nome"]]
        pos = {int(i): k for k, i in enumerate(ids_l)}
        for i, m in zip(ids, mapas):
            if int(i) in pos:
                M.append(m)
                Ml.append(mapas_l[pos[int(i)]])
    if not M:
        return None, None, 0
    g = np.linspace(-W_PM, W_PM, N_PM)
    U, V = np.meshgrid(g, g)
    ref = (U <= -3.0) & (np.abs(V) <= 1.0)
    M, Ml = np.array(M), np.array(Ml)
    norm = lambda A: A / np.nanmean(A[:, ref], 1)[:, None, None]
    with np.errstate(invalid="ignore"):
        return np.nanmean(norm(M), 0), np.nanmean(norm(Ml), 0), len(M)


def em_volta(F, raio=1.4, nth=37):
    """c_s num circulo de `raio` (em raios do pilar) em volta do pilar medio, contra o angulo.

    0 graus = lado do centro da colonia, 180 = lado de fora; as metades de cima e de baixo
    sao promediadas (o pilar medio e simetrico por construcao). E a forma usual de mostrar o
    efeito de um obstaculo sobre um campo (Sherwood/Nusselt local contra o angulo).
    """
    th = np.linspace(0.0, np.pi, nth)
    u = -raio * np.cos(th)
    v = raio * np.sin(th)
    esc = (N_PM - 1) / (2 * W_PM)
    amostra = lambda vv: ndi.map_coordinates(
        F, [(vv + W_PM) * esc, (u + W_PM) * esc], order=1, cval=np.nan
    )
    return np.degrees(th), 0.5 * (amostra(v) + amostra(-v))


def fig_quimica(DS, cores, out):
    """Surfactante em volta do pilar medio (mesmos pilares no liso) e o halo."""
    liso, rug = DS[0], DS[1:]
    linhas = [[f"a{i}", f"p{i}", "h"] for i in range(len(rug))]
    fig, ax = plt.subplot_mosaic(
        linhas, figsize=(10.0, 2.9 * len(rug)), constrained_layout=True
    )
    g = np.linspace(-W_PM, W_PM, N_PM)
    faixa = np.abs(g) <= 0.5
    for i, D in enumerate(rug):
        M, Ml, n = pares_pilar(D, liso)
        if M is None:
            continue
        c = cores[D["nome"]]
        a = ax[f"a{i}"]
        th, vl = em_volta(Ml)
        _, vr = em_volta(M)
        a.plot(th, vl, color="k", lw=1.2, label="liso, mesmo círculo (sem pilar)")
        a.plot(th, vr, color=c, lw=1.6, label=f"{D['nome']}, em volta do pilar")
        a.set_xlim(0, 180)
        a.set_xticks([0, 45, 90, 135, 180])
        a.set_xticklabels(["0°\nfrente\n(centro)", "45°", "90°\nlado", "135°", "180°\ntrás\n(fora)"],
                          fontsize=7)
        a.set_ylabel(r"$c_s$ / $c_s$ que chega do centro")
        a.set_title(f"{D['nome']}: em volta do pilar ({n} pilares × quadros)", loc="left", fontsize=8)
        a.legend(fontsize=6.5)
        ins = a.inset_axes([0.04, 0.06, 0.22, 0.32])
        ins.add_patch(plt.Circle((0, 0), 1.0, color="#989894", lw=0))
        ins.add_patch(plt.Circle((0, 0), 1.4, fill=False, ls="--", lw=0.8, color=c))
        ins.annotate("", xy=(2.1, -1.75), xytext=(0.6, -1.75),
                     arrowprops=dict(arrowstyle="->", lw=0.7))
        ins.text(1.35, -1.55, "fora", fontsize=5.5, ha="center")
        ins.text(-1.65, 0, "0°", fontsize=5.5, ha="right", va="center")
        ins.text(1.65, 0, "180°", fontsize=5.5, ha="left", va="center")
        ins.set_xlim(-2.4, 2.4)
        ins.set_ylim(-2.1, 2.1)
        ins.set_aspect("equal")
        ins.axis("off")

        a = ax[f"p{i}"]
        a.axvspan(-1, 1, color="#989894", alpha=0.45, lw=0, label=f"pilar (só no {D['nome']})")
        with np.errstate(invalid="ignore"):
            a.plot(g, np.nanmean(Ml[faixa], 0), color="k", lw=1.2, label="liso, mesmo lugar")
            a.plot(g, np.nanmean(M[faixa], 0), color=c, lw=1.6, label=D["nome"])
        a.set_xlim(-W_PM, W_PM)
        a.set_xlabel("← centro      raios do pilar      fora →", fontsize=7.5)
        a.set_ylabel(r"$c_s$ / $c_s$ que chega do centro")
        a.set_title("corte pelo centro do pilar", loc="left", fontsize=8)
        a.legend(fontsize=6.5)
    a = ax["h"]
    for D in DS:
        t, v = serie(D, "halo_L")
        a.plot(t, v, **estilo(D, cores), ms=2.5, label=D["nome"])
    a.set_xlabel(r"$t$")
    a.set_ylabel(r"$L/R$")
    a.set_title(r"halo de $c_s$: alcance fora da colônia", loc="left", fontsize=8)
    a.legend()
    salva(fig, out, "fig_serie_quimica")


def fig_painel(DS, out, t_alvo=50.0, L=5.0, n=500):
    """Fig. 3 do artigo: particulas da colonia sobre o campo de nutriente, pilares em cinza.

    Fundo = nutriente `c_n` (verde = intacto, claro = consumido), como o verde do artigo.
    Particulas = colonia (rho_b >= 0.1 ou filler), cor = rho_b; o agar nao e desenhado.
    """
    cmap_n = LinearSegmentedColormap.from_list("nutriente", ["#ffffff", "#6f9e86"])
    fig, ax = plt.subplots(
        1, len(DS), figsize=(2.6 * len(DS) + 0.6, 3.0), constrained_layout=True
    )
    g = np.linspace(-L, L, n)
    GX, GY = np.meshgrid(g, g)
    for a, D in zip(np.atleast_1d(ax), DS):
        t, f = min(quadros(D["run"]), key=lambda p: abs(p[0] - t_alvo))
        d = carrega(f)
        janela = (np.abs(d["x"]) < L + 0.3) & (np.abs(d["y"]) < L + 0.3)
        dist, idx = cKDTree(np.c_[d["x"][janela], d["y"][janela]]).query(
            np.c_[GX.ravel(), GY.ravel()], k=16
        )
        w = spline(dist / (1.8 * DX)) * (d["m"] / d["rho"])[janela][idx]
        cn = (w * d["c_n"][janela][idx]).sum(1) / np.maximum(w.sum(1), 1e-12)
        a.imshow(
            cn.reshape(n, n),
            origin="lower",
            extent=[-L, L, -L, L],
            cmap=cmap_n,
            vmin=0.35,
            vmax=1.0,
            interpolation="bilinear",
        )
        geo = D["geo"]
        if geo is not None:
            for cx, cy in geo["c"]:
                if abs(cx) < L + 0.5 and abs(cy) < L + 0.5:
                    a.add_patch(plt.Circle((cx, cy), geo["R"], color="#989894", lw=0))
        col = ((d["rho_b_grown"] >= 0.1) | (d["is_filler"] > 0.5)) & janela
        sc = a.scatter(
            d["x"][col],
            d["y"][col],
            c=d["rho_b_grown"][col],
            cmap="viridis",
            vmin=0.1,
            vmax=0.8,
            s=0.35,
            lw=0,
            rasterized=True,
        )
        q = min(D["fr"], key=lambda q: abs(q["t"] - t))
        a.set_title(
            f"{D['nome']}, φ ≈ {phi(D):.2f}\n"
            f"N = {int(col.sum())}, $R_{{99}}$ = {q['r99']:.2f}, $D_b$ = {q['Db']:.2f}",
            fontsize=7.5,
        )
        a.set_xlim(-L, L)
        a.set_ylim(-L, L)
        a.set_aspect("equal")
        a.set_xticks([])
        a.set_yticks([])
    cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.01)
    cb.set_label(r"$\rho_b$ da partícula")
    fig.suptitle(f"t = {t:.0f}; fundo: nutriente $c_n$ (verde = intacto)", fontsize=8)
    salva(fig, out, "fig_serie_painel")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="*", help="nome=caminho; o primeiro e o liso")
    ap.add_argument("--out", default="plots")
    ap.add_argument("--refaz", action="store_true")
    a = ap.parse_args()
    lista = [tuple(s.split("=", 1)) for s in a.runs] if a.runs else SERIE

    if geometria(lista[0][1]) is not None:
        raise SystemExit("o primeiro run tem de ser o liso")
    geos = {}
    for nome, run in lista[1:]:
        g = geometria(run)
        if g is None:
            raise SystemExit(f"{run}: sem pilares.json (so o primeiro run pode ser liso)")
        geos[nome] = g

    print("medindo (cache em plots/_cache_serie)")
    DS = [
        mede_run(n, r, geos if i == 0 else {}, a.refaz)
        for i, (n, r) in enumerate(lista)
    ]
    rui = ruido()
    print(f"ruído entre realizações: R99 {100 * rui[0]:.1f}%, dR/dt {100 * rui[1]:.1f}%")

    cm = plt.get_cmap("plasma")
    ph = sorted({D["geo"]["phi"] for D in DS[1:]}, reverse=True)
    cores = {
        D["nome"]: cm(0.15 + 0.6 * ph.index(D["geo"]["phi"]) / max(1, len(ph) - 1))
        for D in DS[1:]
    }

    os.makedirs(a.out, exist_ok=True)
    fig_crescimento(DS, cores, a.out)
    fig_dose(DS, cores, rui, a.out)
    fig_lideres(DS, cores, a.out)
    fig_quimica(DS, cores, a.out)
    fig_painel(DS, a.out)


if __name__ == "__main__":
    main()
