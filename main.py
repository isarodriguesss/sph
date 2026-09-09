import csv
import numpy as np
from scipy import ndimage as ndi
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree
from pysph.solver.application import Application
from pysph.base.kernels import CubicSpline
from pysph.solver.solver import Solver

from src.particles import SEED_MODE, create_initial_state
from src.scheme import MyBiomassScheme


def cubic_spline_w(r, h):
    # W do CubicSpline 2D (PySPH) — usado so na logica de insercao em numpy
    q = np.asarray(r) / h
    fac = 10.0 / (7.0 * np.pi * h * h)
    w = np.zeros_like(q, dtype=float)
    m1 = q <= 1.0
    m2 = (q > 1.0) & (q <= 2.0)
    w[m1] = fac * (1.0 - 1.5 * q[m1] ** 2 + 0.75 * q[m1] ** 3)
    w[m2] = fac * 0.25 * (2.0 - q[m2]) ** 3
    return w


def void_fraction(x, y, rho_b, dx, thresholds=(0.7, 1.0, 1.5), n_grid=200):
    """Fracao da AREA da colonia sem nenhuma particula dentro de thr*dx.

    Metrica primaria do Criterio Obrigatorio de Validacao: ao contrario de sigma_a
    (que so existe onde ha particula), esta mede o vacuo geometrico e portanto
    enxerga a quebra de cobertura espacial exigida pela Particao da Unidade
    (Violeau §3.4; consistencia de interpolacao Liu §3.3).
    """
    colony = rho_b > 0.1
    if int(np.sum(colony)) < 10:
        return {t: 0.0 for t in thresholds}
    r = np.hypot(x, y)
    R = float(np.percentile(r[colony], 99))
    if R <= 0:
        return {t: 0.0 for t in thresholds}
    g = np.linspace(-R, R, n_grid)
    GX, GY = np.meshgrid(g, g)
    # Recorta ao DOMINIO: quando a colonia passa da parede, o disco de raio R cobre
    # regiao sem particula por construcao e isso seria contado como vacuo fisico.
    inside = (
        (GX * GX + GY * GY <= R * R)
        & (np.abs(GX) <= x_max_domain)
        & (np.abs(GY) <= y_max_domain)
    )
    if not np.any(inside):
        return {t: 0.0 for t in thresholds}
    d, _ = cKDTree(np.column_stack([x, y])).query(
        np.column_stack([GX[inside], GY[inside]])
    )
    return {t: float(np.mean(d > t * dx)) for t in thresholds}


LOG_FILE = "log.csv"
LOG_HEADER = [
    "t",
    "iteration",
    "max_v",
    "mean_v",
    "n_fast",
    "a_marangoni",
    "a_drag",
    "a_pressure",
    "a_flag",
    "a_total",
    "min_cs",
    "max_cs",
    "mean_cs",
    "constrast_cs",
    "min_c_n",
    "max_c_n",
    "mean_c_n",
    "contrast_c_n",
    "mass_total",
    "pass_n_spawned",
    "min_sig_bio",
    "mean_sig_bio",
    "frac_lowsig_bio",
    "mean_sig_all",
    "frac_lowsig_all",
    "n_bio_arms",
    "n_ins_arms",
    "void_07",
    "void_10",
    "void_15",
    "a_mar_bio_med",
    "a_mar_bio_p95",
    "cs_bio_arms",
    "c_n_bio_arms",
    "biomass_total",
    "biomass_arms",
    "n_pinned",
    "frac_clump",
    "nn_median",
    "n_shift_gate",
    "c_n_junc",
    "rho_b_junc",
    "rho_b_dip",
    "rho_b_dip_r",
    "n_junc_bio",
]

# Dominio expandido 2026-08-06: em [-5,5] a colonia rompia a parede em t~51s
# (R_p99=6.40, 228 particulas alem de 4.8) — Bloqueio H / licao M-B.10.
# 261 preserva dx: 14/260 = 0.05385 vs 10/186 = 0.05376.
x_dim, y_dim = 261, 261

x_min_domain, x_max_domain = -7.0, 7.0
y_min_domain, y_max_domain = -7.0, 7.0

dx = (x_max_domain - x_min_domain) / (x_dim - 1)


mu = 0.020
gamma = 60.0
# E4b: 5.0 -> 15.0. O E4 cortou `sigma` 10.3 -> 1.93 (para nao saturar pos-maturacao) e
# NAO compensou a amplitude: `a_mar = beta*|grad cs|` caiu de 11.50 para 3.71, a frente
# parou e `R99` foi de 4.58 para 1.85. Como a FORCA e o que se preserva, a razao
# `|F_mar|/B_tension` da licao #24 nao muda — o risco de fratura e o mesmo do E3.
beta = 5.0
# HILL_K — limiar do quorum sensing em `qs = rho_b^2/(rho_b^2 + K^2)`. Ficou em 0.1 desde
# o inicio do projeto: o joelho da curva cai em 10% da densidade de saturacao, ou seja
# `qs(0.1) = 0.5` — MEIA producao com um decimo da biomassa. E o que torna a banda
# sub-quorum quimicamente barulhenta e o que afogou todas as rotas de recrutamento
# (licoes #48, #52, #53, #59). Medido no C4 em t=50, com `sigma` recalibrado para
# preservar a producao da biomassa MADURA (rho_b>0.5, 64% do total em 106 particulas):
#   K=0.3, sigma=11.6 -> banda sub-quorum produz 0.18x; recruta em rho_b=0.05 sai da
#                        saturacao (cs_inf 0.490 -> 0.441)
#   K=0.5, sigma=14.7 -> 0.09x; cs_inf 0.389
# LAMBDA_BIO_RATIO — decaimento de cs no BIOFILME (agar fica em 0.5*lambda, halo
# preservado, licao #42). `cs_inf = P/(lambda_eff + P/cs_max)` satura quando
# `P > 4.5*lambda_eff`. Pos-maturacao P(0.40)=1.437 contra o limiar 1.350 em 2x — a
# colonia satura por 6% de margem. Em 3x o limiar vai a 2.025 e a fracao saturada zera.
# Bonus: L_D = sqrt(D/lambda) encurta, entao |grad cs| ~ cs/L_D fica 22% mais afiado.
# Risco: L_D_int 0.73h -> 0.60h; abaixo de ~1h o SPH resolve mal o gradiente (licao I.3).
# E5 mediu: com a producao do E3, `cs` gruda no teto em QUALQUER razao (53% saturado em
# 2x, 2.5x, 3x e 4x) — `lambda` so e alavanca quando `sigma` ja esta baixo, que e a
# armadilha do E4. Mantido em 2.0.
LAMBDA_BIO_RATIO = 2.0

# E5: 0.15 -> 0.25. O E3 recrutou 3266 particulas e 90% delas ficaram no limbo — elas
# nao contribuem com nada mecanico e sao a fonte da saturacao que derrubou o
# `contrast_cs` para 8.75. Em K=0.25 a banda sub-quorum produz 0.44x e a saturacao
# prevista cai de 53% para 37%, com `sigma` SUBINDO (10.3 -> 11.1), sem o problema de
# amplitude que colapsou o E4.
# Historico: H1 (K=0.3 sobre o C4) deu contraste +15% mas AR 4.97; E2 (K=0.15 sobre o
# N1) deu AR 6.31 — o K MELHORA o AR quando ha nutriente sustentado.
# PROPOSTA 1 do docs/PLANO_C5_CONTINUIDADE.md — PERFIL MONOTONO DE `cs`.
#
# O problema: `cs` e PLANO dentro da colonia (nucleo 0.486, juncao 0.490, mid-arm 0.490,
# rim 0.481). Como a forca e `-beta*grad(cs)`, campo plano da forca ZERO no interior — so
# o rim, onde `cs` cai para o agar, e empurrado. Por isso a expansao e uma vanguarda de 35
# particulas cavalgando 50-77 dx enquanto o corpo anda 7 (razao p90/p50 = 5.0).
#
# Duas causas, ambas necessarias de remover (§3.3.6, quatro configuracoes avaliadas):
#   TETO   — `P/cs_max ~ 14` contra `lambda_eff` = 0.30: o teto E o sumidouro dominante e
#            satura tudo. Sozinho removido, o perfil ainda nao serve (ver c_n abaixo).
#   c_n    — vale 0.1 no nucleo e 0.8 nos bracos, entao suprime a producao justo onde ela
#            deveria ser maxima; cria um bump na juncao (nucleo 17.4 < juncao 24.6).
#            Remover tambem APROXIMA o modelo do Xavier/CCR — ver a tensao ja documentada
#            em §3.2 Frente 2.
#
# Perfil previsto SEM os dois: nucleo 8.3 -> juncao 7.8 -> mid-arm 7.0 -> rim 3.4 ->
# agar 0.23. Unico MONOTONO decrescente: push outward em TODA a colonia (§3.3.4).
#
# `sigma` cai 4.2x porque sem teto o sistema e LINEAR em sigma e `a_mar = beta*|grad cs|`
# so depende do produto `beta*sigma`. Alvo: `a_mar` do corpo no orcamento §8 (~6).
# SUMIDOURO QUADRATICO (P10) — alternativa ao teto. `CS_SINK_K = 0` mantem o teto puro,
# que e o comportamento de toda a serie ate aqui. Ao ligar, poe-se `CS_CEILING = 0`.
#
#   producao = sigma*qs(rho_b)*(1-cs/CS_CEILING)*c_n_f      sumidouro = lambda*cs + K*cs^2
#
# `cs_inf` por zona, medido (licao #77): sob o TETO a razao interior/limbo e **1.05x** —
# uma particula com 20x mais biomassa produz o mesmo `cs` de equilibrio, entao biomassa
# nova nao gera gradiente, so preenche area. Com o sumidouro quadratico a razao vai a
# **2.98x** e o pico cai no MID-ARM (0.83 contra 0.58 do interior e 0.19 do limbo) — a
# linha da §3.3.4 marcada como "dendritica seletiva", o perfil do M-B.10.
#
# `K` fixa a AMPLITUDE sem tocar a ESTRUTURA (medido com K=3/10/30: razao 2.97/2.98/2.99),
# porque `cs ~ sqrt(P/K)`. Com K=10 o maximo e 0.83, contra os 14.3 da remocao pura do teto
# no Y3 (licao #56) — a faixa dinamica cresce como raiz, nao linear.
# P10 (K=10) REPROVADO na morfologia (2026-09-03): a fisica funcionou — mesma biomassa
# (1.1 nos dois em t=45) e `R99` 5.85 contra 3.33 do P5, `contrast_cs` 30 contra 5,
# `max_cs` destravado em 1.13 e ESTAVEL (contra 14.3 da remocao pura do teto no Y3).
# Mas o motor ficou 4x acima do orcamento do §8: `a_mar` mediana **24.34** (P5: 7.68),
# maximo 36.08, razao |F_mar|/B_tension = **887** (P5: 280). E o brittle neck da licao #24
# reproduzido — bracos viram filamentos. Reprovado na leitura visual (§11) antes da medicao.
#
# P11: K = 90. Como `cs ~ sqrt(P/K)`, K escala a AMPLITUDE sem tocar a ESTRUTURA (razao
# interior/limbo fica 2.98x para K=3/10/30). `a_mar` previsto = 24.34/sqrt(9) = **8.1**,
# o nivel do P5 (7.68), com o campo estruturado. Efeito colateral desejado: com o motor no
# ritmo do P5, o `R99` termina perto de 5.5 em t=100 e nao bate na parede util de 6.8 —
# o P10 cruzaria em t~55.
# P14 (2026-09-04) — CALIBRACAO PELO FOCO MORFOLOGICO, nao pela expansao. O P12 (K=90)
# tinha `a_mar` mediana **9.06** contra 4.23 do P5: motor 2x mais forte, e foi isso que
# rarefez (6 dedos em t=100 contra 18, disco 96% agar). O alvo aqui e IGUALAR a amplitude
# do P5 e ficar com o ganho estrutural, que nao depende de K (`cs ~ sqrt(P/K)` escala a
# amplitude; a razao `cs_inf` interior/limbo fica 2.98x para qualquer K).
#
#   K=200 -> max_cs 0.253, a_mar previsto 5.4, razao |F_mar|/B_tension 159 (P5: 154)
#
# HIPOTESE: sob o teto `cs_inf` e PLANO dentro da colonia (razao 1.05x), entao o gradiente
# so existe na fronteira colonia-agar — a forca atua na ponta e o corpo nao e puxado, dai o
# buraco no rastro. Com o sumidouro `cs_inf` varia DENTRO do braco (0.575 interior, **0.828
# mid-arm**, 0.193 limbo): o gradiente ao longo do dedo puxa material do nucleo para o
# mid-arm e empurra o mid-arm para fora. Esteira ao longo do braco, em vez de so a ponta.
CS_SINK_K = 0.0    # ROTA FECHADA (licao #80): 3 pontos, 20x em K, morfologia perdida
# em toda a faixa. Preservado parametrizado (§10) — o mecanismo esta certo, o custo e a forma.

CS_CEILING = 0.5   # base P5 restaurada.
# 0.5 era o valor de toda a serie ate aqui (P1 reprovou remover o teto SEM substituto).
CS_PROD_CN = 1.0   # P1 REPROVADA: c_n de volta. 0 = producao nao depende de c_n
# `cs_max` tinha DOIS papeis: teto da producao e escala do gate da colonizacao
# (`0.6*cs_max`). Sem o teto isso quebra — no P1 pus o gate em `cs > 5.0` e os bracos tem
# `cs` = 0.28, entao DESLIGUEI o recrutamento nos bracos sem querer (VIVAS 178 vs 276).
# Agora e limiar ABSOLUTO: 0.1 fica acima do agar (0.031) e abaixo dos bracos (0.28).
# BASE DE TRABALHO = P5 (decisao da usuaria, 2026-09-03). Sob a definicao de continuidade
# revisada na §2.2 — em que o limbo (`0 < rho_b < 0.1`) conta como tecido conectivo — o P5
# entrega **93.0%** de colonia ligada ao centro contra 59.8% do P4, alem de ocupacao dentro
# do dedo 0.60 contra 0.44 e 67% mais biomassa viva. Custo aceito: `contrast_cs` 10.0 -> 5.7
# e, sob a definicao ESTRITA (so quorum ou filler), continuidade 41.3% -> 19.4%.
#
# P6 (0.40) esta REPROVADO e fecha a direcao oposta: estrangula a colonizacao, que e o UNICO
# termo aditivo do modelo e portanto a fonte de todos os portadores — e sem portadores o
# wake fica sem mae (maes 4904 -> 2348, colonia r>1 2743 -> 1330). Licao #74.
# P12 (2026-09-03): 0.10 -> 0.034. NAO e alavanca nova — e a MESMA alavanca do sumidouro
# aplicada de forma consistente. `COL_CS_MIN` e ABSOLUTO e foi calibrado contra um campo com
# `cs_max`=0.5. Ao trocar o teto pelo sumidouro `K*cs^2`, o campo inteiro reescala por
# `sqrt(K)`: max_cs deu 1.130 (K=10) e 0.381 (K=90). No P11 eu mudei K e esqueci o gate —
# a colonizacao praticamente parou (biomassa cravada em 0.30 de t=10 a t=35, contra 0.78 do
# P5; 179 vivas contra 547) e a rodada testou uma colonia sem recrutamento, nao o sumidouro.
# Mesmo laco de realimentacao para baixo do P6 (licao #74): menos colonia -> menos `cs` ->
# menos elegiveis. 0.10 * (0.381/1.130) = 0.034 preserva a FRACAO do campo que e elegivel.
COL_CS_MIN = 0.10   # base P5. Limiar ABSOLUTO na escala de `cs`: se `CS_SINK_K` mudar,
# ESTE valor tem de ser reescalado junto (licao #78-b, erro do P11).
# Alvo da relaxacao = COL_TARGET_FRAC * doador. Preservado desligado em 1.0 (§10): o
# mecanismo esta certo — separa "nasce colonia" de "nasce limbo" — mas so cabe no orcamento
# de `cs` sob gate alto, e gate alto mata o wake.
COL_TARGET_FRAC = 1.0

# P9 (2026-09-03) — CRESCIMENTO SELETIVO POR MOTILIDADE. O meio-termo entre crescer e
# formar bracos NAO e um escalar entre r_growth 0.02 e 0.06: o P8 terminou com
# `frac(cs>0.45)` em 34% — SEM saturacao — e virou disco assim mesmo, porque o crescimento
# uniforme engorda o corpo inteiro e homogeneiza a fonte de `cs` (licao #76, mesma classe
# da #48/J5). O reforco tem de ser LOCALIZADO.
#
#     r_eff = r_growth * (1 + GROWTH_MOTILE_BOOST * min(|v|/GROWTH_V_SAT, 1))
#
# Discriminador: |v|, o invariante #3 do §9 ("swarmers ativos estao so nas pontas
# avancando"). Medido no P5, `|v|` mediano da biomassa VIVA por zona: nucleo 1.04e-3,
# meio 3.57e-3, 0.7-0.85 1.72e-2, ponta **2.87e-2 — 27.5x o nucleo**; as 10% mais rapidas
# estao em r/R99=0.73, dentro dos bracos. No P8 (blob) a mesma razao e **0.8x**: o
# discriminador some justamente no regime que queremos evitar, o que o torna autoverificavel.
#
# Com boost=3 e v_sat=0.01: nucleo cresce a 0.026 (~base), pontas a 0.08.
#
# NOVO em relacao ao K.7/T1: la o `motile_boost` multiplicava a producao de `cs` (removido
# no T1, e em tensao com Xavier/rhlAB via CCR). Aqui multiplica o CRESCIMENTO, o que e
# biologicamente mais defensavel — swarmers da borda ativa sao os que se dividem.
#
# RISCO declarado: crescer na ponta tem dois efeitos opostos — cria portadores densos (abre
# o gate `|grad rho_b|`, hoje 0.22-0.25 la) mas tambem produz `cs` na ponta, e como o `cs`
# medido la e 0.214 contra 0.459 no interior, isso ACHATA o perfil radial, que foi o que
# matou o P8. A diferenca de mecanismo: o P8 achatou em toda parte; aqui o acrescimo e uma
# casca fina e o corpo continua em 0.02.
# P9 REPROVADO (2026-09-03) com boost=3.0: 1 dedo (era 8), `R99` 3.84 -> 3.03, `dR/dt`
# 0.0742 -> 0.0495, e produziu MAIS biomassa que o P8 uniforme (3.44 vs 3.05) rodando a
# `r_growth` 3x menor no corpo. Realimentacao positiva: rapida cresce -> densa -> abre o
# gate da Marangoni -> acelera -> cresce mais. O discriminador denunciou a si mesmo — a
# razao |v| ponta/nucleo caiu de 27.5x para 7.2x. Licao #77. Preservado desligado (§10).
# P16 (2026-09-09) — DESACOPLAR VOLUME DE BIOMASSA. `d_am` ganha um multiplicador que
# NAO entra em `d_a_rho_b_grown`. Como `rho_b` e DENSIDADE e a producao de `cs` depende de
# `qs(rho_b)` — INTENSIVO —, ganhar massa nao custa `cs`. E a decupagem que as quatro
# refutacoes do crescimento (B2, D5, P8, P9) nunca tiveram: elas subiam `r_growth`, que
# move os dois juntos, e a biomassa extra achatava o campo.
#
# MOTIVO: no P15 a mitose so disparou em t=34 e rendeu 63 divisoes — os bracos se formam
# entre t=15 e t=35, entao ela nao participou da formacao. A massa cresce como
# `dm/dt ~ m*0.0045` (medido: 1.69x em 100 s) e dividir exige 1.3x. Com ganho 3 a massa
# dobra em 51 s e a mitose comeca em **t~19**, dentro da janela de formacao.
#
# Fisica: influxo de van't Hoff ([T2] Srinivasan) — osmolitos secretados puxam agua do
# agar e a colonia ganha VOLUME sem ganhar celulas. A licao #54 rejeitou o influxo porque
# ele "traz solvente, nao celulas" e diluiria; ali faltava o mecanismo que converte massa
# em particula nova, que e exatamente a mitose (P15).
# P17: 3.0 -> 5.0. A massa dobra em 31 s em vez de 51 s. Marcos previstos pela escala
# medida (mitose e dirigida por massa acumulada; `dm/dt` ~ ganho, entao o tempo escala por
# 3/5 = 0.60): 1a divisao t~11.7 (main_00600), ocupacao 0.65 t~21 (main_01200), regime
# preenchido t~39 (main_02400), **tip-splitting t~53 (main_03400)**.
GROWTH_MASS_GAIN = 5.0
# P18: o ganho so vale acima deste limiar — mesmo valor de `MITOSE_RHO_B_MIN`, para que
# ganhe massa exatamente quem pode converte-la em particula nova. Abaixo dele o ganho e 1
# (comportamento do baseline), o que impede o limbo de engordar sem poder se expandir.
GROWTH_GAIN_RHO_B_MIN = 0.05

GROWTH_MOTILE_BOOST = 0.0
GROWTH_V_SAT = 0.01
# Filler como DOADOR (licao #59 revisitada). A recruta nasce na densidade do doador, e nos
# bracos o filler E o material (`rho_b` ~0.3) — sem ele sobram poucas vivas esparsas e a
# recruta nasce em ~0.1, abaixo do quorum. Medido: recruta leva 29 s para cruzar 0.1
# enquanto a frente passa por um ponto em 0.73 s (40x), e 96% nunca cruzam. O requisito e
# NASCER no quorum, nao crescer ate ele.
COL_FILLER_DONOR = 1.0

HILL_K = 0.25
# PROPOSTA 3 do docs/PLANO_C5_CONTINUIDADE.md — distribuir a forca em vez de concentra-la
# numa vanguarda. A `FlagellarForce` tem magnitude CONSTANTE (`f0*gate*n`), entao e imune a
# saturacao de `cs` — precisa so da direcao. O gate `[0.2, 0.6]` seleciona 72 particulas de
# 276 vivas: sao elas que cavalgam 50-77 dx enquanto o corpo anda 7. Alargar distribui a
# forca por ~230. O K.2 so testou ESTREITAR; alargar nunca foi testado.
FLAG_GATE_LO = 0.2  # P3 REPROVADA: alargar p/ [0.1,0.8] levou a razao p90/p50 de 5.0 a 6.8
FLAG_GATE_HI = 0.6  # (forca constante nao move corpo, cria cavaleiros). Preservado parametrizado.
sigma = 11.1  # E5: recalibrado p/ producao da biomassa MADURA constante sob K=0.25
D = 1.5e-3
D_ext = 0.08
lambda_ = 0.15
# REFUTADO sob teto de `cs` (licao #64). Testado em 0.15 duas vezes: D3b (sem a correcao
# de massa e sem S3.3-ii) e D5 (com as duas). O D5 eliminou a causa de massa — `rho/rho0`
# 1.10, massa 203.7 contra 516.8 — e AINDA assim deu `frac(cs>0.45)`=89.6% e `R99`=1.77,
# contra 90.7% e 1.80 do D3b. Enquanto `cs_max` existir, biomassa densa satura `cs` por
# area e mata `∇cs` no interior, qualquer que seja a lei de producao. Religar so depois
# de decidir o teto — e a serie Y mostrou que remove-lo explode a faixa dinamica (#56).
# D6 (2026-08-14): delimita a fronteira entre 0.02 (D4, funciona) e 0.15 (D5, reprovado).
# 0.05 ja cumpre o objetivo da Frente 3 — logistica 0.3→0.8 em 44.7 s, dentro da janela,
# contra 112 s em 0.02 — sem ir ao extremo que satura `cs` em 90% da colonia.
# `R99` e a metrica decisiva: >3.5 significa que ha janela utilizavel sem mexer no teto;
# ~2 significa que o teto e parede dura (licao #64) e a quimica precisa ser redesenhada.
# E4: 0.02 -> 0.05. Com r_growth=0.02 nada amadurece na janela — de rho_b=0.03 ate 0.50
# leva 193 s. Em 0.05 leva 77 s a partir de 0.03 e 49 s a partir de 0.10. Refutado antes
# (D3b/D5) sob lambda_bio=2 e K=0.1, quando a biomassa madura saturava o cs.
# P8 (2026-09-03) testou 0.06 sobre o P5 e REPROVOU — TERCEIRA refutacao da mesma
# alavanca, depois de B2 (licao #43) e D5 (licao #64). Ela entregou tudo que prometia:
# `rho_b` p90 dos portadores no braco 0.19 -> **0.44** (alvo era 0.45), razao max/min da
# largura 2.01 -> **1.46** (cintura/barriga resolvida), continuidade a 1.05 dx 19.4% ->
# **82.9%**, `frac(cs>0.45)` 21% -> 34%, `contrast_cs` 5.67 -> 6.33.
#
# E a colonia deixou de ser dendritica: `R99` 3.84 -> **2.79**, `dR/dt` 0.0742 -> **0.0332**,
# dedos 8 -> **1**. Disco compacto com lobulos — (a) Modulated de Trinschek, criterio de
# falha da §2.2. As tres condicoes que eu julguei terem mudado desde B2/D5 (`k_src`
# sustentando `c_n`, headroom de `cs` nos bracos, `contrast_cs` invalidado como guardrail)
# NAO alteram o mecanismo: mais biomassa achata `grad cs` e a colonia compacta. Licao #76.
r_growth = 0.02
rho_max = 1.0
alpha_mon = 0.12

D_o = 0.04
k_o = 0.5
lambda_o = 0.05
Q0 = 5.0

D_n = 0.05
D_n_int = 1e-4
k_n = 0.5

# N1 — regime nutrient-rich [T2]: reposicao do nutriente pelo agar.
# Equilibrio c_n = k_src/(k_src + k_n*rho_b): 0.3 da c_n=0.67 nos braços (gate de
# crescimento 0.74) contra 0.40 do k_src=0.1, que e o ZERO do gate. 0.0 = desligado.
# P4 (2026-08-14): religado para repor o nutriente que a PROMOCAO consome. O P1 mediu
# `c_n` da juncao caindo 0.997 -> 0.186 com 2400 promovidas, o que fechou o gate do
# ParticleShift (3030 -> 263) e levou `frac_clump` a 0.54. `k_src=0.3` e o valor de
# N1/D3a: equilibrio `c_n = k_src/(k_src + k_n*rho_b)`, `min_c_n` nunca abaixo de 0.375.
# Nao e alavanca nova nem calibracao — e reposicao do sumidouro que a promocao cria.
k_src = 0.3  # E1 APROVADO (=N1): min_c_n 0.000 -> 0.376, a_mar +29%, biomassa +15%,
# AR 5.56, conectividade inalterada. Remove a trava de nutriente. Antes:
# A1 isolava o arrasto: k_src existia para compensar o sumidouro da
# promocao (P4). Sem promocao nao ha sumidouro, e mante-lo seria 2a alavanca.

# X1 — fluxo quimiotatico de biomassa [T3] Giverso: u = chi*(-grad cs), conservativo.
# chi=0.15: X0 calibrou 0.5 com o gradiente do PERFIL RADIAL (0.19), mas a equacao
# usa o gradiente LOCAL SPH, cuja mediana e 0.684 (3.6x maior). 0.099/0.684 = 0.145.
# 0.0 = desligado.
chi = 0.0  # X1/X1b REPROVADOS: transporte com orcamento fixo dilui (licao #54)

dt_global = 0.001
total_sim_time = 60.0   # P17: cobre o tip-splitting previsto em t~53 com 7 s de folga
# mostrou o motor do E5 ainda CRESCENDO em t=100. Checagem antes do sumidouro quadratico.
print_freq = 200

# None = cada rodada tem condicao inicial propria (producao).
# Fixar um inteiro torna a rodada bit-reproduzivel — OBRIGATORIO ao comparar rotas
# (§2.5): sem semente, diferencas de ate ~9 pontos percentuais em metricas de
# amostra pequena (n_bio ~70-90) nao sao atribuiveis ao mecanismo.
SEED = 20260806  # serie K: comparacao com C4 exige mesma semente (§2.5)

trajectory_store_interval = 20

c0 = 0.35

# ============================ MITOSE (P15) ============================
# A FONTE DE VOLUME QUE FALTA. Em [T1] Trinschek `dh/dt = -div(J) + crescimento`: o
# crescimento e termo-fonte de VOLUME e a Marangoni e UM DOS FLUXOS — orienta, nao origina.
# No nosso modelo a Marangoni faz os dois trabalhos e o crescimento nenhum: ele engorda a
# particula (massa sobe -> densidade sobe -> pressao sobe), e o P13 mediu que isso gera
# RESISTENCIA, nao escoamento. Dividir cria volume no espacamento local — a filha ocupa
# espaco vizinho, como uma bacteria que se divide.
#
# CONSERVA BIOMASSA, ao contrario de tudo que reprovou antes: as filhas herdam `rho_b` da
# mae (que e DENSIDADE) e dividem a MASSA, entao `sum(m*rho_b)` nao muda. Nao ha a inflacao
# de biomassa que achatou o `cs` no P8/P9.
#
# O bloco anterior (`use_splitting`) tinha quatro defeitos e nunca rodou de fato:
#   1. `daughter_data = parent_props.copy()` chamado DUAS vezes — a segunda apagava o
#      `m/2`, entao cada uma das 2 filhas levava a massa INTEIRA: a massa DOBRAVA.
#   2. limiar `m > 1.99*m0` nunca dispara — a massa maxima medida em t=100 e **1.695*m0**.
#      (o `m0=0` do HDF5 e artefato: `m0` nao esta em `add_output_arrays`.)
#   3. copiava 14 de 70 propriedades — filhas herdariam LIXO em `is_filler`, `c_n`, `noise`,
#      `x_dep`... e uma filha com `is_filler=1` nasceria congelada.
#   4. offsets aleatorios INDEPENDENTES: as duas filhas podiam nascer coincidentes.
#
# Sizing (P5_t100, t=100, so biomassa viva): com limiar 1.3 ha 2820 maduras, **96% NAO
# pinadas** (o pin por `rho_b>=0.8` pega ZERO). Distribuicao radial: 276 em r/R99<0.2,
# **1297 em 0.2-0.4**, 60 em 0.8-1.0 — a divisao ocorre no CORPO, nao na frente, porque a
# massa cresce como `m*rho_b*(1-rho_b)` e e CUMULATIVA.
# P17: mantido em 1.3 (nao 1.1). Medido: com limiar 1.1 a filha nasce em 0.55*m0 e leva
# ~125 s para voltar a ser elegivel — cada particula dividiria UMA vez na janela. Com 1.3 a
# filha nasce em 0.65*m0 e volta em ~62 s: DUAS divisoes por particula em 60 s. Antecipa o
# inicio em so 3 s a menos e dobra a cadencia, que e o que de fato preenche.
MITOSE_M_RATIO = 1.3      # divide quando m > MITOSE_M_RATIO * m0
MITOSE_FREQ = 200         # iteracoes entre chamadas
MITOSE_MAX = 150          # cap por chamada
MITOSE_EPS = 0.35         # separacao das filhas: +-MITOSE_EPS*h, direcoes OPOSTAS
MITOSE_RHO_B_MIN = 0.05   # nao divide traco
# P17: gate de DENSIDADE — nao divide onde `rho/rho0` ja esta alto, porque ali nao falta
# volume, SOBRA. Sem ele o P16 (ganho 3) levou `rho/rho0` p99 a 5.93 e o `dt` a 1/4 do
# inicial em t=58.8; escalado para ganho 5 isso cairia em **t~35**, quatro segundos ANTES
# do primeiro marco que queremos observar (regime preenchido, t~39) e dezoito antes do
# tip-splitting (t~53). Sem o gate a rodada empaca antes de responder a pergunta.
#
# Limiar 2.0: o baseline P5 (sem mitose) opera a `rho/rho0` p99 = 4.98 no fim e funciona,
# entao o gate nao pode ser agressivo. Em t=50 o P16 esta em p99 = 3.12 — o corpo tipico
# fica bem abaixo de 2 e a mitose segue livre onde ha espaco; o gate so morde nos
# aglomerados que puxam o p99 e estrangulam o `dt`.
# GATE DE DENSIDADE — DESLIGADO (P18). Medido no P16: as 118 particulas maduras tem `rho`
# mediano **1.69** e **ZERO** acima de 4, enquanto o limbo tem p99 = 8.22 com 42% acima de
# 4. **Quem divide nao e quem esta comprimido**, entao o gate nao toca a fonte da
# compressao: em 4.0 e no-op (passa 100%), em 2.0 barra 25% arbitrariamente (foi o P17,
# que cortou a mitose em 83% e nao testou o desenho). Substituido pelo limiar no GANHO
# (`GROWTH_GAIN_RHO_B_MIN`), que ataca o mecanismo medido. Preservado parametrizado (§10).
MITOSE_RHO_MAX = 99.0     # em unidades de rho0; 99 = desligado
RHO0_NOMINAL = 1.0        # = o `rho0` da BiomassEOS em src/scheme.py

prob_of_splitting = 0.03  # usado so pelo bloco antigo, preservado desligado
use_splitting = False     # bloco antigo (bugado) — mantido desligado, ver §10
use_mitose = True

use_shift = True
SHIFT_COEFF = 0.5
SHIFT_CAP = 0.0006
SHIFT_RHO_B_MIN = 0.1

use_kgc = True
KGC_DET_MIN = 0.25

# S4 = transparencia quimica de cs (sempre ativa em equations.py)
# S5 = + filler nao consome nutriente (alavanca isolada)
FILLER_NUTRIENT_TRANSPARENT = 0

# Rota D1 — difusao de biomassa DENTRO da fase densa, p/ curar o degrau da junção.
# Gate rho_b_i*rho_b_j anula o fluxo na interface colonia-agar (nao borra os braços).
# 0.0 = desligado. Estimativa: L_D = sqrt(D_b*t); 0.008 da ~0.5 em 30s no corpo denso.
D_b = 0.0  # D1 REPROVADO: destruiu o nucleo (rho_b 1.0->0.48, n_pinned 43->0)

# K3 — colonizacao com alvo-DOADOR (docs/PLANO_K3_JUNCAO.md, licao #59). O alvo e a
# densidade da MAE (media ponderada por rho_b), nao a Shepard; filler fora dos doadores;
# gate cs > 0.6*cs_max. Converte de fato (agar invadido 72.5%->52.1%) mas ~92% do
# convertido para no limbo sub-quorum, porque crescer de 0.15 a 0.5 leva 124 s contra a
# janela de 50 s. Religar so depois de resolver a maturacao (defeito D, licao #58).
# E3 (2026-08-28): religado sobre E1+E2. O gate de nutriente da colonizacao (`c_n>0.4`)
# estava FECHADO no C4 — so 36 particulas elegiveis em t=50, porque `c_n` colapsa. Com
# `k_src=0.3` sao 712. E o alvo-doador mediano e 0.233, ACIMA do quorum, com 100% dos
# alvos >=0.1: sob k_col=0.03 (tau=33s) a recruta chega a ~0.18 em 50 s.
# P6 testou 0.20 (tau 33s -> 5s, para a recruta cruzar o quorum antes da frente passar) e
# REPROVOU junto com o gate 0.40 — mas `k_col` nao e alavanca independente do gate: ele so
# acelera a chegada ao mesmo equilibrio, fixado pelo gate. Licao #74.
k_col = 0.03

# P1 (2026-08-14) — promocao do limbo sub-quorum ao quorum. Ver bloco em post_step.
# P2 (virar filler com rho_b=0.45) REPROVADO: `OxigenConsumption` nao isenta filler,
# entao levar 2529 particulas de rho_b~1e-91 para 0.45 criou um sumidouro de nutriente
# que nao existia; `mean_c_n` 0.999->0.941, o pin quimico `c_n<0.6` disparou e a colonia
# congelou em t<6 (mean_v 3.1e-4 -> 2e-6). Ver a tabela da licao #61.
# A1 (2026-08-14) — arrasto do agar. `gamma = 60` se aplicava a TODAS as particulas,
# inclusive as de `rho_b=0`, o que dava ao agar velocidade terminal 4% da colonia
# (medido: `|v|` mediano 4.9e-29). Somado a `|p| = 0` exato da EOS, o agar nao podia
# ser empurrado NEM arrastado — so engolido. `gamma` e friccao flagelo-substrato,
# propriedade da bacteria; aplica-la ao meio e o que o torna fundo rigido.
# 0.1 -> gamma=6 no agar; velocidade terminal ~10x maior.
# ROTA B, segunda metade — o arrasto do agar. Fisicamente `gamma` e friccao
# flagelo-substrato, propriedade da BACTERIA; aplica-la ao meio e o que o transforma em
# fundo rigido (ver docstring da `LinearDrag`).
#
# As duas metades NAO sao independentes, e cada uma sozinha ja falhou:
#   A1  — baixou para 0.1 SEM pressao no agar: o arrasto era o unico resistente a
#         compressao, entao `rho/rho0` do agar engolido foi de 4.8 para 16.7 e `R99` de
#         4.62 para 3.06.
#   B1  — pressao (AGAR_FADE=1) COM arrasto cheio: o agar passou a resistir (`|p|` 0 ->
#         0.12 em 21268 particulas, `rho/rho0` max 5.86 -> 3.63) mas NAO passou a ser
#         empurrado — adveccao em `r>2` ficou em 31 particulas, identica ao E5, e
#         `C5b` = 9.4% contra 11.3% do baseline.
# B2 testa a combinacao, que e a celula que sobrou: pressao resiste a compressao, arrasto
# baixo deixa o meio acompanhar a frente.
AGAR_DRAG_RATIO = 1.0  # B2 REPROVADO (ver acima): preservado desligado

# ROTA B do docs/PLANO_C5_CONTINUIDADE.md — ramo SO repulsivo no meio sem biomassa
# (`fade_rep = AGAR_FADE`, `fade_att = 0` para `rho_b < 0.1`).
#
# Licao #73: o agar NAO e empurrado — deslocamento mediano 0.005 dx, so 26% se move mais
# de 1 dx — porque tem `|p| = 0` EXATO mesmo esmagado a `rho/rho0 = 4.8`. Um meio sem
# pressao nao transmite empurrao, entao a colonia o ATRAVESSA em vez de desloca-lo, e a
# adveccao transporta apenas 31 das 1343 particulas dos bracos (98% chegam la por
# deposicao do wake). E o mecanismo por tras da violacao do C5.
#
# O ramo ATRATIVO fica em zero de proposito: dar coesao a regiao sub-densa CONTRAI a
# colonia (licao #66-C, P3: `R99` 4.62 -> 1.93), porque abaixo de `rho/rho0 = 1` a EOS so
# tem esse ramo. O que falta e resistencia a COMPRESSAO, nao coesao.
#
# Historico dos testes anteriores, ambos invalidos para esta pergunta:
#   A2  — rodou com o bug de Group da `BiomassEOS` (o ramo repulsivo NUNCA executava,
#         `p>0` em 0 de 70982 particulas), entao AGAR_FADE multiplicava algo que nao
#         rodava e o run saiu bit-identico ao C4. Nao foi refutado, nao foi testado.
#   A2b — ja sobre a EOS corrigida, mediu "zero efeito na desjuncao" — mas por OCUPACAO
#         em t~29, antes do C5 existir, e sobre config anterior ao HILL_K/k_col.
# P13 (2026-09-04) — O CRESCIMENTO COMO FONTE DE VOLUME. Em [T1] Trinschek a governante e
# `dh/dt = -div(J_conv + J_Mar) + crescimento`: o crescimento e termo-fonte de VOLUME, gera
# pressao, a pressao gera fluxo, e a Marangoni e UM DOS FLUXOS — ela orienta, nao origina.
# [T2] Srinivasan idem. No nosso modelo a Marangoni faz os DOIS trabalhos e o crescimento
# nenhum: `m -> rho -> p -> forca` existe, mas a `BiomassEOS` zera `fade_rep` abaixo de
# `rho_b=0.1`. Medido em t=50: **83.4% do corpo tem pressao EXATAMENTE zero**, e a fracao
# cai com o raio (62% no nucleo, **4.1% na ponta**). Crescer no limbo — 90% do braco — nao
# empurra nada.
#
# Isso explica as duas observacoes da usuaria de uma vez: sem fonte de volume a Marangoni
# move PARTICULA A PARTICULA (so as que passam no gate), e num continuo com fonte a
# conservacao OBRIGA o meio a acompanhar — dai os buracos no rastro; e com crescimento alto
# a biomassa achata `cs`, o motor morre e sobra so o inchaco (P8/P9).
#
# `agar_fade` liga o ramo REPULSIVO abaixo do quorum mantendo `fade_att = 0`. Distincao do
# que ja reprovou: a licao #66-C (P2/P3) deu COESAO (ramo ATRATIVO) a regiao sub-densa e a
# colonia contraiu (`R99` 4.62 -> 1.93); a rota B testou repulsao-apenas mas no AGAR
# (`rho_b=0`, que nao cresce) e julgada por C5, nunca por expansao.
#
# PRE-REQUISITO VERIFICADO: `p = B*excess^2` so age onde `rho/rho0 > 1`. Medido, o limbo
# esta a **1.506 de mediana com 98.4% acima de 1** — mola carregada. Tem de estar ligado
# desde t=0: num run formado liberaria de uma vez a energia de 14 000 particulas.
# P13 REPROVADO (2026-09-04): o meio pressurizado RESISTE em vez de escoar. Expoente em
# `dR/dt ~ R^p` foi de -2.15 para **-2.64** (criterio era subir acima de -1.6), e em
# t[75,100] a expansao caiu a METADE (0.0100 contra 0.0195). Ajuda enquanto a colonia e
# pequena (1.25x em t[50,75]) e cobra caro quando a frente fica longa — cruzamento em t~75.
# Morfologia praticamente identica ao P5. Ganhos que ficam: 21 dedos contra 18, e agar
# engolido 437 contra 740 (-41%) — o deslocamento acontece, so custa mais do que rende.
# `a_pressure` 3.07 vs 3.05, picos>4 em 2% nos dois: o risco de enrijecer NAO se
# materializou. Licao #79. Preservado desligado (§10).
AGAR_FADE = 0.0

# Conversao do agar ENGOLIDO em MATRIZ PASSIVA (licao #66-A/E). O agar dentro do
# envelope da colonia nunca vira colonia: 97% da area em t=50. Converter para filler
# (nao para biomassa viva) e a unica rota quimicamente grátis — o filler e gateado
# fora da SurfactantEquation como fonte E destino. [T2]: fase passiva = matriz/EPS.
# `is_matrix` marca as convertidas para que NAO sirvam de semente da proxima rodada,
# senao a conversao vira flood-fill para o agar aberto.
# PONTE DE CONECTIVIDADE. Medido no C4: o corpo (`rho_b>0.1`) parte em 26 componentes a
# partir de t=17.8, e ligar TODOS os bracos ao nucleo custa 31 particulas em t=50 (pico de
# 49 em t=41), mediana de 2 por braco, todas agar. E 1.1% do corpo — contra as 7553 que o
# preenchimento do envelope pedia. Recruta por CONECTIVIDADE (caminho minimo no grafo),
# nao por raio nem por vizinhanca, que e o que espalhava o custo do `k_col`.
# PISO DE ENVELOPE — troca o `rho_b` do agar que a colonia JA ocupou por um valor
# visivel em escala log, sem tocar na fisica. As particulas ja estao la; o zero exato
# vem do underflow da gaussiana inicial (licao #48), nao de fisica, e faz o painel
# mostrar um esqueleto de fios em vez da massa continua da `reference.jpg`.
#
# Inercia verificada por calculo (8712 particulas no envelope, eps=1e-3):
#   quimica  : +0.53% na producao de cs  (qs(1e-3, K=0.25) = 1.6e-5)
#   mecanica : GRATIS — EOS, ParticleShift e FlagellarForce ja excluem rho_b<0.1
#   metabolica: ZERO — `OxigenConsumption` isenta a banda [0, 0.01)
# O teste F1 da licao #48 ja mostrou que um piso inerte sai bit-identico.
# DESLIGADO (2026-09-01). O piso resolvia um problema de IMAGEM mexendo no SOLVER, o que
# a licao #38 proibe. O preenchimento topologico foi movido para a RENDERIZACAO
# (`tools/plot_piso.py --fill`), onde entrega a mesma continuidade visual com perturbacao
# ZERO. Codigo preservado desligado (§10): com `use_floor=False` o `is_env` fica 0 em toda
# parte e as sete isencoes em `src/equations.py` sao no-op, entao a fisica e a do E5.
# Medido no solver antes de desligar (E10, topologico 3.5dx, contra E9/E5):
#   buraco 650 -> 369 e agar engolido 1001 -> 216, mas `a_pressure` mediana 2.64 -> 3.24,
#   `mean_v` a menor da serie e `frac(cs>0.45)` de volta a 53.4%. Nao e de graca.
use_floor = False
# 0.1 EXATO: o smoothstep da BiomassEOS comeca em 0.1, entao `fade(0.1) = 0.0000` — a
# particula conta como corpo em qualquer metrica e tem coesao e pressao NULAS por
# construcao. O flag `is_env` isenta das outras seis: crescimento, gradiente (gate da
# Marangoni), producao de cs, consumo de nutriente, shifting e doacao na colonizacao.
# Sobra `LinearDrag` (gamma +0.9 em 60, 1.5%) e a difusao de cs, que e desejavel: o
# marcador conduz surfactante como meio, so nao produz.
RHO_B_FLOOR = 0.1
FLOOR_FREQ = 100
# PREENCHIMENTO TOPOLOGICO, nao criterio de vizinhanca. Rasteriza o corpo, fecha vaos
# de ate FLOOR_FECHA*dx, inunda a partir de FORA: o que a inundacao nao alcanca e buraco.
# Nao tem raio de busca, entao nao existe o vazamento que qualquer criterio local sofre —
# baia e ligada ao exterior por construcao e NUNCA e marcada. Medido em `runs/E9` contra
# o enclausuramento a 6/8 que ele substitui, com MENOS particulas (3565 vs 4152):
#   buraco cercado  715 -> 405      dedos 17 -> 20      largura do braco 0.760 -> 0.618
#   agar engolido  1090 -> 241      agar limpo na baia 44.5% -> 48.7% (sem piso: 49.8%)
# O enclausuramento vazava porque o raio (6.9 dx) e comparavel a largura da baia no anel
# medio; iterar so piorava (17 dedos -> 5 em um passo).
FLOOR_FECHA = 3.5     # vao maximo fechado antes da inundacao, em multiplos de dx
FLOOR_CELL = 0.5      # lado da celula da grade, em multiplos de dx

use_bridge = False  # REPROVADO: metrica circular (ver licao #68)
BRIDGE_FREQ = 100
BRIDGE_VALUE = 0.3      # fade da EOS = 0.5 (coesao real); em 0.1 exato o fade e ZERO
BRIDGE_MIN_COMP = 8     # componente menor que isto e ruido, nao braco
BRIDGE_MAX_COST = 6     # nao construir ponte longa: acima disso o braco esta solto mesmo
BRIDGE_LINK = 1.5       # multiplos de h para a aresta do grafo

use_matrix = False  # serie M REVERTIDA: reprovou no §11 (regressao fingering -> modulated)
MATRIX_VALUE = 0.3   # fade da EOS = 0.5; equilibrio c_n = k_src/(k_src+k_n*rho_b)
MATRIX_FREQ = 50
MATRIX_ENV = 1.5  # M4 (ponto de operacao). M6 testou 3.0: converteu 96% do agar e a
# colonia virou disco compacto (1 dedo, AR 0.13, ocupacao 2.05) — extremo do trade-off.

use_promo = False  # A1 isola o arrasto do agar: promocao DESLIGADA
PROMO_VALUE = 0.12
PROMO_RHO_B_MAX = 0.1  # alvo: 0 < rho_b <= este valor (o limbo)
PROMO_R_MAX = 1.8  # coroa uniforme; alem disso o limbo vira spokes nos braços

use_insert = True
INSERT_FREQ = 200
INSERT_SIGMA_TRIG = 0.85
INSERT_RHO_B_MIN = (
    0.5  # C3.4 validado: filler frozen so no nucleo estrutural (licao #31/#35)
)
INSERT_PROX = 0.7
INSERT_MAX = 100

use_wake = True
WAKE_FREQ = 100  # baseline E5; a serie C1-C3 nao achou efeito atribuivel entre 25 e 100
WAKE_DISP = 1.0
WAKE_PROX = 0.7
WAKE_RHO_B_MIN = 0.05
# Piso de coesao da fase passiva (matriz/EPS, [T2] Srinivasan): o filler herda o rho_b
# da mae, e mae sub-quorum gera filler com fade_rep=fade_att=0 na BiomassEOS — inerte
# para sempre, porque filler nao cresce. Medido no C4: fade mediano 0.016 na juncao.
FILLER_RHO_B_FLOOR = 0.0  # teste do piso concluido (runs/F_floor); desligado p/ isolar a conversao
WAKE_MAX = 150
WAKE_MODE = 2
WAKE_CLUSTER_MAX = 7
WAKE_RING_RATIO = 0.75
WAKE_MASS_BUDGET = 0.12
# ROTA C do docs/PLANO_C5_CONTINUIDADE.md — fecha o RASTRO entre a posicao antiga e a
# atual da ponta, nao so a antiga. O wake deposita 1 ponto (mais anel) enquanto a ponta
# percorre 2-7 dx, entao o braco nasce como colar de grumos: medido no E5, 2410 particulas
# de wake em 612 componentes de mediana 2 na escala de contato. Isso e o mecanismo por tras
# do C5 (§2.2) — a colonia nao expande, ela e construida em pulsos (licao #73).
#
# Diferenca para o B1 (licao #65, REPROVADO por AR -25%): o B1 depositava round(d/dx)
# pontos ao longo do segmento SEMPRE. Aqui cada ponto passa pelo mesmo teste de
# proximidade do resto do wake, entao so entra onde ha vao de fato, e o anel de
# adensamento (n_extra) fica so no ponto de origem.
# ANCORAGEM (proposta 2026-09-02). O `WAKE_PROX` rejeita deposito PERTO demais (evita
# empilhar) e nao havia limite SUPERIOR: se a ponta correu e deixou um vazio de 3 dx, o
# wake depositava no meio dele e nascia uma ilha. Medido no E5: **69% das particulas do
# wake nascem a mais de 1.05 dx de qualquer material** (mediana 1.44, cauda ate 3.7).
#
# E de la que vem a descontinuidade: entre t=3.7 e t=8.3 a conectividade cai de 76% para
# 44% e 185 particulas de wake aparecem desconectadas; dai em diante o wake e 71-76% de
# todo o material fora do corpo conexo.
#
# A leitura fisica: a celula-filha nasce ADJACENTE a mae, nao no meio do espaco vazio.
# Como `added_pts` tambem conta como ancora, o deposito ao longo do rastro (WAKE_SEG)
# constroi CADEIA: o primeiro ponto ancora no corpo, o segundo no primeiro.
WAKE_ATTACH = 1.05  # 0 = desligado

# P7 — REDISTRIBUIR o rastro (nao adicionar). Medido no P5: a barriga do braco fica
# cravada em r/R99 = 0.80-0.90 do t=12 ao t=50 enquanto o raio absoluto vai de 1.11 a
# 3.07 — ela ACOMPANHA a frente, a ~15% de R99 atras dela. Nao e evento historico num
# raio fixo: e a zona de deposicao viajando colada a ponta. Largura em r/R99=0.8 e
# 5.0 dx contra 2.2 dx em 0.6 (razao max/min 2.1x, igual nos 8 bracos).
#
# Hoje o wake poe 1 + n_extra particulas TODAS na posicao vagada (`x_dep`), em anel —
# um grumo. Com WAKE_SPREAD as MESMAS maes depositam a MESMA contagem, distribuida ao
# longo do segmento `x_dep -> x`. Mesma massa, mesmo gate, posicoes diferentes.
#
# NAO confundir com `WAKE_SEG` (licao #65, REPROVADO, AR 10.48 -> 7.83): aquele
# ADICIONAVA particulas ao longo do segmento alem do grumo. Este REALOCA. A diferenca
# entre adicionar e realocar e o que separa engrossar de uniformizar.
#
# P7 REPROVADO (2026-09-03): o mecanismo funcionou — as deposicoes ficaram visivelmente
# espalhadas — mas a razao max/min da largura PIOROU, 2.14x -> 2.47x (barriga 4.8 -> 5.5,
# cintura igual em 2.2). Causa: o wake e ~3200 particulas contra ~14000 de limbo no braco,
# e e o limbo que forma o manto cuja espessura varia. Redistribuir 3200 nao muda a
# distribuicao dos 14000. A barriga e engarrafamento do limbo (|v| 4.5e-3 empilhando
# contra 8.2e-4), problema de FORCA, nao de granularidade da deposicao. Licao #75.
# Preservado desligado (§10): ele foi o unico lever a subir a continuidade a 1.05 dx
# (19.4% -> 24.2%), porque depositar espalhado adiciona por ADJACENCIA em vez de criar ilhas.
WAKE_SPREAD = False

WAKE_SEG = False  # rota C REPROVADA (C5b 12->19 contra alvo 70). Preservado desligado:
# a deposicao dirigida e barata (AR -8% contra -25% do B1), mas nao resolve C5
WAKE_SEG_MAX = 6
# Rota C (ABORTADA 2026-08-06) — realocar agar ocioso conservaria massa, MAS o agar
# tem pressao ZERO (BiomassEOS: fade_rep=fade_att=0 para rho_b<0.1), entao o buraco
# deixado pelo doador NAO cicatriza: cada doacao e uma puncao permanente no campo.
# Codigo preservado para referencia; so reativar se o agar ganhar resposta de pressao.
WAKE_RECYCLE = False
WAKE_DONOR_MARGIN = 0.6

use_pass_n = False
PASS_N_FREQ = 100
PASS_N_MAX_PARENTS = 100
PASS_N_SIGMA_TRIG = 0.95
PASS_N_ALPHA = 0.75
PASS_N_EPSILON = 0.5
PASS_N_MAX_GEN = 1
PASS_N_RHO_B_MIN = 0.15
PASS_N_RHO_B_MAX = 0.95
PASS_N_PROXIMITY_MIN = 0.4


class SwarmApp(Application):
    def initialize(self):
        with open(LOG_FILE, "w", newline="") as f:
            csv.writer(f).writerow(LOG_HEADER)
        self._m_initial = (
            None  # snapshot da massa total em t=0 (preenchido em post_step)
        )
        self._pass_n_spawned_since_log = 0  # acumula spawns entre linhas de log
        self._wake_mass_added = 0.0  # so o que o WAKE adicionou (nao BiomassGrowth)
        self._promo_total = 0
        self._matrix_total = 0
        self._bridge_total = 0

    def create_particles(self):
        if SEED is not None:
            np.random.seed(SEED)
        fluid_solid = create_initial_state(
            x_dim,
            y_dim,
            rho_max,
            dt_global,
            x_min=x_min_domain,
            x_max=x_max_domain,
            y_min=y_min_domain,
            y_max=y_max_domain,
            seed=SEED,
        )

        for pa in fluid_solid:
            if pa.name == "fluid":
                pa.add_property("noise")
                pa.noise[:] = (
                    1.0
                    + 0.6 * np.sin(SEED_MODE * np.arctan2(pa.y, pa.x))
                    + 0.01 * np.random.rand(len(pa.x))
                )
                pa.add_property("dt_force")
                pa.add_property("dt_cfl")
                pa.add_property("au_mar")
                pa.add_property("ax_mar")
                pa.add_property("ay_mar")
                pa.add_property("au_drag")
                pa.add_property("grad_rho_b_x")
                pa.add_property("grad_rho_b_y")
                pa.add_property("grad_rho_b_mag")
                pa.add_property("grad_cs_x")
                pa.add_property("grad_cs_y")
                pa.add_property("grad_co_x")
                pa.add_property("grad_co_y")
                pa.add_property("sigma_a")
                pa.sigma_a[:] = 1.0
                pa.add_property("rho_b_smooth")
                pa.add_property("rho_b_w2")
                pa.add_property("shift_dC_x")
                pa.add_property("shift_dC_y")
                pa.add_property("shift_x")
                pa.add_property("shift_y")
                pa.shift_x[:] = 0.0
                pa.shift_y[:] = 0.0
                pa.add_property("Mxx")
                pa.add_property("Mxy")
                pa.add_property("Myx")
                pa.add_property("Myy")
                pa.add_property("Lxx")
                pa.add_property("Lxy")
                pa.add_property("Lyx")
                pa.add_property("Lyy")
                pa.Lxx[:] = 1.0
                pa.Lxy[:] = 0.0
                pa.Lyx[:] = 0.0
                pa.Lyy[:] = 1.0
                pa.add_property("is_filler")
                pa.is_filler[:] = 0.0
                pa.add_property("is_wake")
                pa.is_wake[:] = 0.0
                pa.add_property("is_matrix")
                pa.is_matrix[:] = 0.0
                pa.add_property("is_conv")
                pa.is_conv[:] = 0.0
                pa.add_property("is_env")
                pa.is_env[:] = 0.0
                pa.add_property("rho_b_pre")
                pa.rho_b_pre[:] = 0.0
                pa.add_property("x_dep")
                pa.add_property("y_dep")
                pa.x_dep[:] = pa.x[:]
                pa.y_dep[:] = pa.y[:]
                pa.add_property("gen")
                pa.gen[:] = 0.0
                pa.add_property("ax_drag")
                pa.add_property("ay_drag")
                pa.add_property("au_flag")
                pa.add_property("x_spawn_ref")
                pa.add_property("y_spawn_ref")
                pa.x_spawn_ref[:] = pa.x[:]
                pa.y_spawn_ref[:] = pa.y[:]
                pa.add_output_arrays(
                    [
                        "rho_b_grown",
                        "cs",
                        "c_o",
                        "u",
                        "v",
                        "p",
                        "noise",
                        "au_flag",
                        "au_mar",
                        "sigma_a",
                        "is_wake",
                        "is_matrix",
                        "is_conv",
                        "is_env",
                        "rho_b_pre",
                        "is_filler",
                        "c_n",
                        "m",
                        "rho",
                    ]
                )
            elif pa.name == "solid":
                pa.add_property("p")

        return fluid_solid

    def create_scheme(self):
        return MyBiomassScheme(
            fluids=["fluid"],
            solids=["solid"],
            dim=2,
            mu=mu,
            gamma=gamma,
            beta=beta,
            sigma=sigma,
            D=D,
            D_ext=D_ext,
            lambda_=lambda_,
            r_growth=r_growth,
            rho_max=rho_max,
            c0=c0,
            alpha_mon=alpha_mon,
            D_o=D_o,
            k_o=k_o,
            lambda_o=lambda_o,
            Q0=Q0,
            D_n=D_n,
            D_n_int=D_n_int,
            k_n=k_n,
            k_src=k_src,
            cs_max=CS_CEILING,
            cs_sink_k=CS_SINK_K,
            col_cs_min=COL_CS_MIN,
            col_filler_donor=COL_FILLER_DONOR,
            col_target_frac=COL_TARGET_FRAC,
            growth_motile_boost=GROWTH_MOTILE_BOOST,
            growth_mass_gain=GROWTH_MASS_GAIN,
            growth_gain_rho_b_min=GROWTH_GAIN_RHO_B_MIN,
            growth_v_sat=GROWTH_V_SAT,
            prod_cn=CS_PROD_CN,
            chi=chi,
            use_shift=use_shift,
            shift_coeff=SHIFT_COEFF,
            shift_cap=SHIFT_CAP,
            shift_rho_b_min=SHIFT_RHO_B_MIN,
            use_kgc=use_kgc,
            kgc_det_min=KGC_DET_MIN,
            filler_nutrient_transparent=FILLER_NUTRIENT_TRANSPARENT,
            D_b=D_b,
            k_col=k_col,
            hill_k=HILL_K,
            flag_gate_lo=FLAG_GATE_LO,
            flag_gate_hi=FLAG_GATE_HI,
            lambda_bio_ratio=LAMBDA_BIO_RATIO,
            agar_drag_ratio=AGAR_DRAG_RATIO,
            agar_fade=AGAR_FADE,
        )

    def create_solver(self):
        kernel = CubicSpline(dim=2)
        solver = Solver(
            dim=2,
            integrator=self.scheme.get_integrator(),
            kernel=kernel,
            dt=5e-5,
            adaptive_timestep=True,
            cfl=0.4,
        )
        solver.tf = total_sim_time
        solver.set_print_freq(print_freq)
        return solver

    def post_step(self, solver):
        if self._m_initial is None:
            self._m_initial = float(np.sum(self.particles[0].m))

        # P4 — P1 (promocao do limbo) + `k_src`. O P1 e o unico dos tres que ENCHEU a
        # juncao: `rho_b_dip` 0.0000 -> 0.1398 e `n_junc_bio` 37 -> 1391, sustentado os
        # 50 s. Ele falhou pelo custo, nao pelo mecanismo: as 2400 promovidas consomem
        # nutriente (`OxigenConsumption` e proporcional a `rho_b`), `c_n` da juncao caiu
        # de 0.997 para 0.186, o gate do ParticleShift (`c_n>=0.6`) fechou de 3030 para
        # 263 e o `frac_clump` foi de 0.00 a 0.54.
        #
        # `k_src=0.3` repoe exatamente o que a promocao consome. Nao e alavanca nova:
        # medida em N1/D3a com `min_c_n` nunca abaixo de 0.375 e motor na frente +23%.
        #
        # P2 (limbo -> filler 0.45) e P3 (envelope 1.5h) REPROVADOS — ver runs/. O P3
        # e o mais informativo: converteu so 1941 particulas e ja bastou para `R99`
        # ficar em 1.93 contra 4.62. A colonia nao absorve nem 2000 particulas novas,
        # o que aponta capacidade de absorcao, nao mecanismo de conversao.
        if use_floor and solver.count % FLOOR_FREQ == 0:
            fl = self.particles[0]
            # O piso NAO e portadora de si mesmo, e e recalculado do zero: sem isso a
            # aureola avanca um raio de busca por chamada e vira flood-fill (medido:
            # R99 do piso 0.78 -> 2.13 em 23 s, sempre a frente do R99 real). Mesmo
            # papel do `is_conv` da serie M.
            # O reset devolve a particula ao estado que ela teria SEM o piso:
            # `rho_b_pre` (o valor no momento da marcacao) mais o que a colonizacao
            # acrescentou por cima do piso. Zerar direto destruia 102 recrutamentos
            # no E8 — e era isso, nao o arrasto, que deslocava os dedos.
            _old = fl.is_env > 0.5
            fl.rho_b_grown[_old] = fl.rho_b_pre[_old] + np.maximum(
                0.0, fl.rho_b_grown[_old] - RHO_B_FLOOR
            )
            fl.is_env[_old] = 0.0
            _corpo = (fl.rho_b_grown >= 0.1) | (fl.is_filler > 0.5)
            if int(np.count_nonzero(_corpo)) > 20:
                _dx = FLOOR_CELL * float(fl.h[0]) / 1.8
                _rc = 1.2 * float(
                    np.percentile(np.hypot(fl.x[_corpo], fl.y[_corpo]), 99)
                )
                _nc = int(2.0 * _rc / _dx) + 1
                _ix = np.clip(((fl.x + _rc) / _dx).astype(int), 0, _nc - 1)
                _iy = np.clip(((fl.y + _rc) / _dx).astype(int), 0, _nc - 1)
                _A = np.zeros((_nc, _nc), dtype=bool)
                _A[_ix[_corpo], _iy[_corpo]] = True
                _k = int(np.ceil(FLOOR_FECHA / FLOOR_CELL))
                _st = np.ones((2 * _k + 1, 2 * _k + 1), dtype=bool)
                _buraco = ndi.binary_fill_holes(
                    ndi.binary_closing(_A, structure=_st)
                ) & (~_A)
                _pv = (
                    (fl.rho_b_grown < RHO_B_FLOOR)
                    & (fl.is_filler < 0.5)
                    & _buraco[_ix, _iy]
                )
                if np.any(_pv):
                    fl.rho_b_pre[_pv] = fl.rho_b_grown[_pv]
                    fl.rho_b_grown[_pv] = RHO_B_FLOOR
                    fl.is_env[_pv] = 1.0
                # marcador que a colonizacao promoveu de verdade deixa de ser marcador
                fl.is_env[fl.rho_b_grown > 1.5 * RHO_B_FLOOR] = 0.0

        if use_bridge and solver.count % BRIDGE_FREQ == 0:
            fl = self.particles[0]
            _rb = fl.rho_b_grown
            _r = np.hypot(fl.x, fl.y)
            _body_all = _rb > 0.1
            if int(_body_all.sum()) > 40:
                _R = float(np.percentile(_r[_body_all], 99))
                _sel = np.where(_r < 1.15 * _R)[0]
                _X, _Y = fl.x[_sel], fl.y[_sel]
                _B = _rb[_sel] > 0.1
                _n = len(_sel)
                _h0 = float(fl.h[0])
                _p = cKDTree(np.column_stack([_X, _Y])).query_pairs(
                    BRIDGE_LINK * _h0, output_type="ndarray"
                )
                if len(_p) > 0:
                    _pb = _p[_B[_p[:, 0]] & _B[_p[:, 1]]]
                    _nc, _lab = connected_components(
                        coo_matrix(
                            (np.ones(len(_pb)), (_pb[:, 0], _pb[:, 1])), shape=(_n, _n)
                        ),
                        directed=False,
                    )
                    _lb = np.where(_B, _lab, -1)
                    _tam = np.bincount(_lb[_B], minlength=_nc)
                    _ic = int(np.argmin(np.hypot(_X, _Y) + 1e9 * (~_B)))
                    _nuc = _lb[_ic]
                    _alvos = [
                        i for i in range(_nc)
                        if _tam[i] >= BRIDGE_MIN_COMP and i != _nuc
                    ]
                    if _alvos:
                        # custo de ENTRAR num no: 0 se ja e corpo, 1 se precisa recrutar
                        _w = np.where(_B[_p[:, 1]], 0.0, 1.0)
                        _w2 = np.where(_B[_p[:, 0]], 0.0, 1.0)
                        _g = coo_matrix(
                            (
                                np.concatenate([_w, _w2]),
                                (
                                    np.concatenate([_p[:, 0], _p[:, 1]]),
                                    np.concatenate([_p[:, 1], _p[:, 0]]),
                                ),
                            ),
                            shape=(_n, _n),
                        ).tocsr()
                        _d, _pred, _ = dijkstra(
                            _g,
                            indices=np.where(_lb == _nuc)[0],
                            return_predecessors=True,
                            min_only=True,
                        )
                        _ponte = set()
                        for _i in _alvos:
                            _ix = np.where(_lb == _i)[0]
                            _j = _ix[int(np.argmin(_d[_ix]))]
                            if not np.isfinite(_d[_j]) or _d[_j] > BRIDGE_MAX_COST:
                                continue
                            _c = _j
                            while _c >= 0 and _pred[_c] >= 0:
                                if not _B[_c]:
                                    _ponte.add(int(_c))
                                _c = _pred[_c]
                        if _ponte:
                            _gi = _sel[np.array(sorted(_ponte), dtype=int)]
                            fl.rho_b_grown[_gi] = np.maximum(
                                fl.rho_b_grown[_gi], BRIDGE_VALUE
                            )
                            self._bridge_total += len(_gi)
                            print(
                                f"ponte t={solver.t:.1f}s: {len(_gi)} recrutadas "
                                f"({len(_alvos)} bracos soltos, total={self._bridge_total})"
                            )

        if use_matrix and solver.count % MATRIX_FREQ == 0:
            fl = self.particles[0]
            _seed = np.where((fl.rho_b_grown > 0.1) & (fl.is_conv < 0.5))[0]
            if len(_seed) > 10:
                _h0 = float(fl.h[0])
                _tb = cKDTree(np.column_stack([fl.x[_seed], fl.y[_seed]]))
                _dmin, _ = _tb.query(np.column_stack([fl.x, fl.y]))
                _conv = (
                    (fl.is_filler < 0.5)
                    & (fl.rho_b_grown < 1e-12)
                    & (_dmin < MATRIX_ENV * _h0)
                )
                if np.any(_conv):
                    fl.rho_b_grown[_conv] = MATRIX_VALUE
                    fl.is_filler[_conv] = 1.0
                    fl.is_wake[_conv] = 1.0
                    fl.is_matrix[_conv] = 1.0
                    fl.is_conv[_conv] = 1.0
                    self._matrix_total += int(_conv.sum())

        if use_promo:
            fl = self.particles[0]
            _r = np.hypot(fl.x, fl.y)
            _m = (
                (fl.is_filler < 0.5)
                & (fl.rho_b_grown > 0.0)
                & (fl.rho_b_grown <= PROMO_RHO_B_MAX)
                & (_r < PROMO_R_MAX)
            )
            if np.any(_m):
                fl.rho_b_grown[_m] = PROMO_VALUE
                self._promo_total += int(_m.sum())

        # Print stats
        if solver.count % print_freq == 0:
            fluid = self.particles[0]
            v_mag = np.sqrt(fluid.u**2 + fluid.v**2)
            max_v = np.max(v_mag)
            mean_v = np.mean(v_mag)
            n_fast = int(np.sum(v_mag > 0.1))

            a_mar = np.max(np.abs(fluid.au_mar))
            a_drag = np.max(np.abs(fluid.au_drag))
            a_total = np.max(np.sqrt(fluid.au**2 + fluid.av**2))
            ax_p = fluid.au - fluid.ax_mar - fluid.ax_drag
            ay_p = fluid.av - fluid.ay_mar - fluid.ay_drag
            a_pressure = np.max(np.sqrt(ax_p**2 + ay_p**2))
            a_flag = np.max(np.abs(fluid.au_flag))

            # Filler e quimicamente transparente: seu cs fica CONGELADO no valor
            # herdado e nao representa o campo. Incluí-lo infla mean_cs e deprime
            # contrast_cs artificialmente.
            # o piso tambem: ele CONDUZ cs mas nao produz, entao seu cs e baixo
            # (mediana 0.170 contra 0.481 nas vivas). Incluí-lo derruba mean_cs pela
            # metade e DOBRA contrast_cs — artefato medido no E8.
            _chem = (fluid.is_filler < 0.5) & (fluid.is_env < 0.5)
            _cs_real = fluid.cs[_chem] if np.any(_chem) else fluid.cs
            min_cs = np.min(_cs_real)
            max_cs = np.max(_cs_real)
            mean_cs = np.mean(_cs_real)
            contrast_cs = (max_cs - min_cs) / (mean_cs + 1e-9)

            min_c_n = np.min(fluid.c_n)
            max_c_n = np.max(fluid.c_n)
            mean_c_n = np.mean(fluid.c_n)
            contrast_c_n = (max_c_n - min_c_n) / (mean_c_n + 1e-9)

            mass_total = float(np.sum(fluid.m))

            # C2 (§2.5) exige a colonia INTEIRA, inseridas incluidas: a banda dos
            # braços e populacao GEOMETRICA (phi_s), enquanto bio_mask e biologica
            # (rho_b puro). Pos-D2 o filler tem rho_b=0 e sairia da banda.
            _phi_s = fluid.rho_b_grown
            # o piso e marcador visual, nao biomassa: em `biomass_total` ele valia
            # 0.72 de 1.05 no E8, escondendo que a biomassa real CAIU 28%.
            _real = (fluid.is_filler < 0.5) & (fluid.is_env < 0.5)
            arms_mask = (_phi_s >= 0.1) & (_phi_s < 0.5)
            bio_mask = (
                (fluid.rho_b_grown >= 0.1)
                & (fluid.rho_b_grown < 0.5)
                & (fluid.is_filler < 0.5)
            )
            n_bio_arms = int(np.sum(bio_mask))
            n_ins_arms = int(np.sum(arms_mask & (fluid.is_filler > 0.5)))

            if n_bio_arms > 0:
                sig_bio = fluid.sigma_a[bio_mask]
                min_sig_bio = float(np.min(sig_bio))
                mean_sig_bio = float(np.mean(sig_bio))
                frac_lowsig_bio = float(np.mean(sig_bio < 0.85))
            else:
                min_sig_bio = mean_sig_bio = 1.0
                frac_lowsig_bio = 0.0

            if int(np.sum(arms_mask)) > 0:
                sig_all = fluid.sigma_a[arms_mask]
                mean_sig_all = float(np.mean(sig_all))
                frac_lowsig_all = float(np.mean(sig_all < 0.85))
            else:
                mean_sig_all = 1.0
                frac_lowsig_all = 0.0

            # Envelope da colonia = suporte mecanico (celula viva OU matriz), senao o
            # filler sai do envelope ao perder rho_b e o vazio vira artefato.
            vf = void_fraction(fluid.x, fluid.y, fluid.rho_b_grown, dx)

            # Motor medido SO na biomassa real: 'a_marangoni' acima e um MAXIMO
            # (invariante #1 do §9 — max e cego ao tipico).
            # Biomassa REAL (exclui filler) e contagem de pinados: n_pinned e o
            # teste direto da maturacao dos braços em nucleo (licao #18/#41).
            # Clumping (Liu §6.4): resolvido o vacuo, o defeito remanescente e
            # sobreposicao, nao falta de particula. Mede-se pelo vizinho mais proximo.
            _rr = np.hypot(fluid.x, fluid.y)
            if np.any(_phi_s > 0.1):
                _Rc = float(np.percentile(_rr[_phi_s > 0.1], 99))
            else:
                _Rc = 1.0
            _colony = (_phi_s > 0.05) & (_rr > 0.3 * _Rc)
            if int(np.sum(_colony)) > 10:
                _d, _ = cKDTree(np.column_stack([fluid.x, fluid.y])).query(
                    np.column_stack([fluid.x[_colony], fluid.y[_colony]]), k=2
                )
                _nn = _d[:, 1] / dx
                frac_clump = float(np.mean(_nn < 0.5))
                nn_median = float(np.median(_nn))
            else:
                frac_clump = 0.0
                nn_median = 1.0

            # Junção nucleo-braço (r 0.4-1.2): onde o degrau de rho_b aparece.
            # rho_b_dip = minimo do perfil AO LONGO dos braços (nao media azimutal,
            # que mistura braço com baia e da a falsa impressao de desconexao).
            # Pos-D2 estas metricas sao medidas SO nas PORTADORAS vivas. O anel e
            # dominado por agar invadido (rho_b=0 e is_filler=0 — licao #58), entao
            # um percentil sobre a populacao inteira mede composicao, nao biologia:
            # dava 0.0004 no D2 e 0.31 no C4 sem que a biomassa viva mudasse.
            _carrier = _real & (fluid.rho_b_grown > 1e-6)
            _ju = (_rr >= 0.4) & (_rr < 1.2)
            c_n_junc = float(np.mean(fluid.c_n[_ju])) if np.any(_ju) else 0.0
            _juc = _ju & _carrier
            rho_b_junc = (
                float(np.percentile(fluid.rho_b_grown[_juc], 90))
                if np.any(_juc)
                else 0.0
            )
            n_junc_bio = int(np.sum(_ju & _real & (fluid.rho_b_grown > 0.1)))
            _th = np.arctan2(fluid.y, fluid.x)
            # Braço detectado por GEOMETRIA, com a banda ESCALADA por _Rc e o limiar no
            # quorum. As duas correcoes vem de falhas medidas: `rho_b>0.3` cego quando a
            # promocao poe todos em 0.12, e `r in [2,3]` fixo cego quando a colonia
            # encolhe (no P1, R99=3.05 e a banda tinha ZERO particulas -> dip nunca
            # calculado e reportado como 1.0, "sem degrau" sem ninguem medir).
            # `rho_b_dip_r` diz ONDE o vale esta — sem isso um dip=0 nao localiza nada.
            _arm = (_rr > 0.55 * _Rc) & (_rr < 0.75 * _Rc) & (_phi_s > 0.1)
            rho_b_dip = 1.0
            rho_b_dip_r = 0.0
            if int(np.sum(_arm)) > 20:
                _h, _e = np.histogram(_th[_arm], bins=72, range=(-np.pi, np.pi))
                _step = max(0.1, 0.05 * _Rc)
                for _b in np.argsort(_h)[-4:]:
                    _c = 0.5 * (_e[_b] + _e[_b + 1])
                    _in = np.abs(
                        ((_th - _c + np.pi) % (2 * np.pi)) - np.pi
                    ) < np.deg2rad(9)
                    for _r0 in np.arange(0.15 * _Rc, 0.75 * _Rc, _step):
                        _m = _in & (np.abs(_rr - _r0) < _step)
                        if not np.any(_m):
                            continue
                        _mc = _m & _carrier
                        _v = (
                            float(np.max(fluid.rho_b_grown[_mc]))
                            if np.any(_mc)
                            else 0.0
                        )
                        if _v < rho_b_dip:
                            rho_b_dip = _v
                            rho_b_dip_r = float(_r0)

            # quantas particulas o ParticleShift efetivamente processa (gate de c_n)
            n_shift_gate = int(
                np.sum(
                    (_phi_s >= SHIFT_RHO_B_MIN) & (_phi_s < 0.8) & (fluid.c_n >= 0.6)
                )
            )

            _vol = fluid.m[_real] / np.maximum(fluid.rho[_real], 1e-9)
            biomass_total = float(np.sum(fluid.rho_b_grown[_real] * _vol))
            n_pinned = int(np.sum(fluid.rho_b_grown >= 0.8))

            _bio = (fluid.rho_b_grown > 0.1) & (fluid.is_filler < 0.5)
            if int(np.sum(_bio)) > 0:
                _r = np.hypot(fluid.x, fluid.y)
                _R = float(np.percentile(_r[_phi_s > 0.1], 99))
                _arm = _bio & (_r > 0.3 * _R)
                a_mar_bio_med = float(np.median(fluid.au_mar[_bio]))
                a_mar_bio_p95 = float(np.percentile(fluid.au_mar[_bio], 95))
                if int(np.sum(_arm)) > 0:
                    cs_bio_arms = float(np.mean(fluid.cs[_arm]))
                    c_n_bio_arms = float(np.mean(fluid.c_n[_arm]))
                    biomass_arms = float(
                        np.sum(
                            fluid.rho_b_grown[_arm]
                            * fluid.m[_arm]
                            / np.maximum(fluid.rho[_arm], 1e-9)
                        )
                    )
                else:
                    cs_bio_arms = c_n_bio_arms = biomass_arms = 0.0
            else:
                a_mar_bio_med = a_mar_bio_p95 = cs_bio_arms = c_n_bio_arms = 0.0
                biomass_arms = 0.0

            print("-" * 50)
            print(f"Tempo: {solver.t:.2f}s | Iteração: {solver.count}")
            print(f"Velocidade Máx: {max_v:.4f}")
            print(f"Contraste CS: {contrast_cs:.4f}")
            print(
                f"c_n: mean={mean_c_n:.4f} max={max_c_n:.4f} contrast={contrast_c_n:.2f} | massa: {mass_total:.2f}"
            )
            print(
                f"sigma_a braços biomassa (n={n_bio_arms}): min={min_sig_bio:.3f} "
                f"mean={mean_sig_bio:.3f} frac<0.85={frac_lowsig_bio:.2%} | "
                f"c/ inseridas (n={n_ins_arms}): mean={mean_sig_all:.3f} "
                f"frac<0.85={frac_lowsig_all:.2%}"
            )
            print(
                f"vazio areal: >0.7dx={vf[0.7]:.2%} >1.0dx={vf[1.0]:.2%} "
                f">1.5dx={vf[1.5]:.2%}  (alvo: >1.5dx -> 0)"
            )
            print(
                f"motor na biomassa: a_mar med={a_mar_bio_med:.3f} "
                f"p95={a_mar_bio_p95:.2f} | bracos: cs={cs_bio_arms:.4f} "
                f"c_n={c_n_bio_arms:.3f}   (S0: med=3.38 p95=9.46 cs=0.353 c_n=0.552)"
            )
            print(
                f"biomassa real: total={biomass_total:.4f} bracos={biomass_arms:.4f}"
                f" | n_pinned(rho_b>=0.8)={n_pinned}"
            )
            print(
                f"empacotamento: clump(<0.5dx)={frac_clump:.1%} "
                f"nn_mediana={nn_median:.3f}dx | shift processa {n_shift_gate} part."
            )
            print(
                f"junção (so vivas): c_n={c_n_junc:.3f} "
                f"rho_b_p90={rho_b_junc:.3f} n_bio={n_junc_bio} | "
                f"vale ao longo do braço={rho_b_dip:.3f} em r={rho_b_dip_r:.2f}   (alvo: > 0.5)"
            )
            print("Acelerações:")
            print(f"  > Marangoni (líq): {a_mar:.2f}")
            print(f"  > Drag:            {a_drag:.2f} (Freio)")
            print(f"  > Pressão (est):   {a_pressure:.2f}")
            print(f"  > Total |a|:       {a_total:.2f}")

            with open(LOG_FILE, "a", newline="") as f:
                csv.writer(f).writerow(
                    [
                        f"{solver.t:.4f}",
                        solver.count,
                        f"{max_v:.6f}",
                        f"{mean_v:.6f}",
                        n_fast,
                        f"{a_mar:.4f}",
                        f"{a_drag:.4f}",
                        f"{a_pressure:.4f}",
                        f"{a_flag:.4f}",
                        f"{a_total:.4f}",
                        f"{min_cs:.4f}",
                        f"{max_cs:.4f}",
                        f"{mean_cs:.4f}",
                        f"{contrast_cs:.4f}",
                        f"{min_c_n:.4f}",
                        f"{max_c_n:.4f}",
                        f"{mean_c_n:.4f}",
                        f"{contrast_c_n:.4f}",
                        f"{mass_total:.6e}",
                        self._pass_n_spawned_since_log,
                        f"{min_sig_bio:.4f}",
                        f"{mean_sig_bio:.4f}",
                        f"{frac_lowsig_bio:.4f}",
                        f"{mean_sig_all:.4f}",
                        f"{frac_lowsig_all:.4f}",
                        n_bio_arms,
                        n_ins_arms,
                        f"{vf[0.7]:.4f}",
                        f"{vf[1.0]:.4f}",
                        f"{vf[1.5]:.4f}",
                        f"{a_mar_bio_med:.4f}",
                        f"{a_mar_bio_p95:.4f}",
                        f"{cs_bio_arms:.5f}",
                        f"{c_n_bio_arms:.4f}",
                        f"{biomass_total:.5f}",
                        f"{biomass_arms:.5f}",
                        n_pinned,
                        f"{frac_clump:.4f}",
                        f"{nn_median:.4f}",
                        n_shift_gate,
                        f"{c_n_junc:.4f}",
                        f"{rho_b_junc:.4f}",
                        f"{rho_b_dip:.4f}",
                        f"{rho_b_dip_r:.3f}",
                        n_junc_bio,
                    ]
                )
                self._pass_n_spawned_since_log = 0  # reseta após registrar

        if use_mitose and solver.count > 0 and solver.count % MITOSE_FREQ == 0:
            fluid = self.particles[0]
            m_nom = dx * dx

            madura = (
                (fluid.m > MITOSE_M_RATIO * m_nom)
                & (fluid.rho_b_grown > MITOSE_RHO_B_MIN)
                & (fluid.is_filler < 0.5)
                & (fluid.is_env < 0.5)
                & (fluid.rho < MITOSE_RHO_MAX * RHO0_NOMINAL)
            )
            idx = np.where(madura)[0]

            if len(idx) > MITOSE_MAX:
                # as mais maduras primeiro — divide quem ja dobrou mais
                idx = idx[np.argsort(-fluid.m[idx])][:MITOSE_MAX]

            if len(idx) > 0:
                # TODAS as propriedades persistentes. Copiar so um subconjunto foi o
                # defeito 3 do bloco antigo: a filha herdava lixo do realloc em
                # `is_filler`, `c_n`, `noise`, `x_dep`. Acumuladores (a_*, au, av, grad_*,
                # L*, M*, shift_*) sao recomputados a cada passo e vao a zero.
                herda = [
                    "rho_b_grown", "cs", "c_n", "c_o", "noise", "h", "rho", "m0",
                    "u", "v", "phi_m" if "phi_m" in fluid.properties else "rho_b_grown",
                    "is_filler", "is_wake", "is_matrix", "is_conv", "is_env", "gen",
                ]
                herda = list(dict.fromkeys(p for p in herda if p in fluid.properties))

                ang = np.random.rand(len(idx)) * 2.0 * np.pi
                off = MITOSE_EPS * fluid.h[idx]
                # filhas em direcoes OPOSTAS: centro de massa preservado
                nx = np.concatenate([fluid.x[idx] + off * np.cos(ang),
                                     fluid.x[idx] - off * np.cos(ang)])
                ny = np.concatenate([fluid.y[idx] + off * np.sin(ang),
                                     fluid.y[idx] - off * np.sin(ang)])

                data = {"x": nx, "y": ny}
                for p in herda:
                    data[p] = np.concatenate([fluid.get(p)[idx], fluid.get(p)[idx]])
                # MASSA CONSERVADA: cada filha leva metade. Era o defeito 1.
                data["m"] = np.concatenate([fluid.m[idx] * 0.5, fluid.m[idx] * 0.5])
                # o rastro do wake recomeca na posicao da filha
                for p in ("x_dep", "y_dep"):
                    if p in fluid.properties:
                        data[p] = nx.copy() if p == "x_dep" else ny.copy()
                # nasce viva por construcao
                if "is_filler" in data:
                    data["is_filler"] = np.zeros(2 * len(idx))

                filhas = fluid.empty_clone()
                filhas.add_particles(**{k: np.asarray(v) for k, v in data.items()})
                fluid.append_parray(filhas)
                fluid.remove_particles(idx)
                solver.nnps.update()

                self._mitose_since_log = getattr(self, "_mitose_since_log", 0) + len(idx)
                print(f"mitose t={solver.t:.1f}s: {len(idx)} divisoes")

        if use_splitting:
            if solver.count > 0 and solver.count % 500 == 0:
                fluid = self.particles[0]
                min_rho_b_for_division = 0.05
                mature_by_mass_indices = np.where(fluid.m > 1.99 * fluid.m0)[0]
                mature_by_rho_b_indices = np.where(
                    fluid.rho_b_grown[mature_by_mass_indices] > min_rho_b_for_division
                )[0]
                indices_to_split = []
                for idx in mature_by_rho_b_indices:
                    original_idx = mature_by_mass_indices[idx]
                    if np.random.rand() < prob_of_splitting:
                        indices_to_split.append(original_idx)

                if len(indices_to_split) > 0:
                    print(
                        f"\n--- Divisão Celular em t={solver.t:.2f}: {len(indices_to_split)} partículas se dividindo. ---"
                    )

                    daughters = fluid.empty_clone()

                    props_to_copy = [
                        "x",
                        "y",
                        "m",
                        "h",
                        "rho",
                        "rho_b_grown",
                        "cs",
                        "u",
                        "v",
                        "au",
                        "av",
                        "a_rho_b_grown",
                        "a_c_s",
                        "m0",
                    ]

                    for parent_idx in indices_to_split:
                        parent_props = {
                            prop: getattr(fluid, prop)[parent_idx]
                            for prop in props_to_copy
                        }

                        for i in range(2):
                            daughter_data = parent_props.copy()
                            daughter_data["m"] = parent_props["m"] / 2.0
                            daughter_data["m0"] = parent_props["m0"]
                            # daughter_data["rho_b_grown"] = parent_props["rho_b_grown"] / 2.0
                            dx_local = parent_props["h"] * 0.5
                            offset_x = dx_local * (np.random.rand() - 0.5)
                            offset_y = dx_local * (np.random.rand() - 0.5)

                            daughter_data = parent_props.copy()
                            daughter_data["x"] = parent_props["x"] + offset_x
                            daughter_data["y"] = parent_props["y"] + offset_y

                            data_to_add = {
                                key: [value] for key, value in daughter_data.items()
                            }

                            daughters.add_particles(**data_to_add)

                    fluid.append_parray(daughters)

                    fluid.remove_particles(indices_to_split)

                    solver.nnps.update()

        if use_pass_n and solver.count > 0 and solver.count % PASS_N_FREQ == 0:
            fluid = self.particles[0]

            colony_mask = (fluid.rho_b_grown > PASS_N_RHO_B_MIN) & (
                fluid.rho_b_grown < PASS_N_RHO_B_MAX
            )

            splittable_gen_mask = fluid.gen < PASS_N_MAX_GEN

            sigma_a = fluid.sigma_a
            split_mask = (
                colony_mask & (sigma_a < PASS_N_SIGMA_TRIG) & splittable_gen_mask
            )
            split_idx = np.where(split_mask)[0]

            if len(split_idx) > 0:
                if len(split_idx) > PASS_N_MAX_PARENTS:
                    order = np.argsort(sigma_a[split_idx])
                    split_idx = split_idx[order][:PASS_N_MAX_PARENTS]

                eps = PASS_N_EPSILON
                alpha = PASS_N_ALPHA
                prox_min = PASS_N_PROXIMITY_MIN * dx
                prox_min_sq = prox_min * prox_min

                # Padrão hexagonal: 6 vértices a 60°
                angles_hex = np.arange(6) * (np.pi / 3.0)
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                tree_mask = np.ones(len(fluid.x), dtype=bool)
                tree_mask[split_idx] = False
                positions = np.column_stack([fluid.x[tree_mask], fluid.y[tree_mask]])
                tree = cKDTree(positions)

                daughters = fluid.empty_clone()
                new_positions = []  # filhas já adicionadas nesta call
                mothers_used = []  # índices de mães que efetivamente splitaram
                n_daughters_total = 0

                # Snapshot dos campos (índices estáveis até o remove)
                x_arr = fluid.x[split_idx]
                y_arr = fluid.y[split_idx]
                h_arr = fluid.h[split_idx]
                m_arr = fluid.m[split_idx]
                u_arr = fluid.u[split_idx]
                v_arr = fluid.v[split_idx]
                rho_arr = fluid.rho[split_idx]
                rhob_arr = fluid.rho_b_grown[split_idx]
                cs_arr = fluid.cs[split_idx]
                co_arr = fluid.c_o[split_idx]
                cn_arr = fluid.c_n[split_idx]
                noise_arr = fluid.noise[split_idx]
                sigma_arr = fluid.sigma_a[split_idx]
                gen_arr = fluid.gen[split_idx]

                for k in range(len(split_idx)):
                    x_m, y_m = float(x_arr[k]), float(y_arr[k])
                    h_m = float(h_arr[k])
                    m_m = float(m_arr[k])
                    r_off = eps * h_m

                    valid_vertices = []
                    for j in range(6):
                        vx = x_m + r_off * cos_a[j]
                        vy = y_m + r_off * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox_min:
                            continue
                        too_close = False
                        for nx, ny in new_positions:
                            if (vx - nx) ** 2 + (vy - ny) ** 2 < prox_min_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        valid_vertices.append((vx, vy))

                    n_d = 1 + len(valid_vertices)
                    if n_d != 7:
                        continue

                    m_d = m_m / 7.0
                    h_d = alpha * h_m
                    u_d = float(u_arr[k])
                    v_d = float(v_arr[k])

                    xs = [x_m] + [p[0] for p in valid_vertices]
                    ys = [y_m] + [p[1] for p in valid_vertices]

                    new_positions.append((x_m, y_m))
                    new_positions.extend(valid_vertices)

                    data = {
                        "x": xs,
                        "y": ys,
                        "m": [m_d] * 7,
                        "h": [h_d] * 7,
                        "rho": [float(rho_arr[k])] * 7,
                        "rho_b_grown": [float(rhob_arr[k])] * 7,
                        "cs": [float(cs_arr[k])] * 7,
                        "c_o": [float(co_arr[k])] * 7,
                        "c_n": [float(cn_arr[k])] * 7,
                        "u": [u_d] * 7,
                        "v": [v_d] * 7,
                        "noise": [float(noise_arr[k])] * 7,
                        "sigma_a": [float(sigma_arr[k])] * 7,
                        "is_filler": [0.0] * 7,
                        "gen": [float(gen_arr[k]) + 1.0] * 7,
                        "is_wake": [0.0] * 7,
                        "is_matrix": [0.0] * 7,
                        "is_conv": [0.0] * 7,
                        "is_env": [0.0] * 7,
                        "rho_b_pre": [0.0] * 7,
                        "x_dep": xs,
                        "y_dep": ys,
                    }
                    daughters.add_particles(**data)
                    mothers_used.append(int(split_idx[k]))
                    n_daughters_total += 7

                if mothers_used:
                    fluid.append_parray(daughters)
                    fluid.remove_particles(np.asarray(mothers_used, dtype=np.uint32))
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_daughters_total
                    print(
                        f"Pass N v2.4 t={solver.t:.1f}s: "
                        f"{len(mothers_used)} mães → {n_daughters_total} filhas "
                        f"(n_d=7 estrito, ε={eps}, α={alpha}, σ_trig={PASS_N_SIGMA_TRIG}, gen≤1)"
                    )

        if use_insert and solver.count > 0 and solver.count % INSERT_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx

            _phi_s = fluid.rho_b_grown
            void_mask = (_phi_s > INSERT_RHO_B_MIN) & (
                fluid.sigma_a < INSERT_SIGMA_TRIG
            )
            void_idx = np.where(void_mask)[0]

            if len(void_idx) > 0:
                positions = np.column_stack([fluid.x, fluid.y])
                tree = cKDTree(positions)
                prox = INSERT_PROX * dx  # distância de segurança entre duas partículas
                prox_sq = prox * prox  # o quadrado da distância de segurança

                angles_hex = np.arange(6) * (
                    np.pi / 3.0
                )  # angulos hexagonais 0°, 60°, 120°, 180°, 240°, 300°
                cos_a = np.cos(angles_hex)
                sin_a = np.sin(angles_hex)

                new_x = []
                new_y = []
                parent_idx = []
                added_pts = []  # dedupe entre candidatos da mesma call

                for k in void_idx:
                    xk = float(fluid.x[k])
                    yk = float(fluid.y[k])
                    for j in range(6):
                        # coordenadas polares
                        vx = xk + dx * cos_a[j]
                        vy = yk + dx * sin_a[j]
                        d_existing, _ = tree.query([vx, vy])
                        if d_existing < prox:
                            continue
                        too_close = False
                        for ax_, ay_ in added_pts:
                            if (vx - ax_) ** 2 + (vy - ay_) ** 2 < prox_sq:
                                too_close = True
                                break
                        if too_close:
                            continue
                        new_x.append(vx)
                        new_y.append(vy)
                        parent_idx.append(int(k))
                        added_pts.append((vx, vy))
                    if len(new_x) >= INSERT_MAX:
                        break

                if len(new_x) > 0:
                    n_ins = len(new_x)
                    p = np.asarray(parent_idx, dtype=int)
                    inserted = fluid.empty_clone()
                    data = {
                        "x": new_x,
                        "y": new_y,
                        "m": [m_target] * n_ins,
                        "h": list(fluid.h[p]),
                        "rho": list(fluid.rho[p]),
                        "rho_b_grown": list(np.maximum(fluid.rho_b_grown[p], FILLER_RHO_B_FLOOR)),
                        "cs": list(fluid.cs[p]),
                        "c_o": list(fluid.c_o[p]),
                        "c_n": list(fluid.c_n[p]),
                        "u": [0.0] * n_ins,
                        "v": [0.0] * n_ins,
                        "noise": list(fluid.noise[p]),
                        "is_filler": [1.0] * n_ins,
                        "is_wake": [0.0] * n_ins,
                        "is_matrix": [0.0] * n_ins,
                        "is_conv": [0.0] * n_ins,
                            "is_env": [0.0] * n_ins,
                            "rho_b_pre": [0.0] * n_ins,
                        "x_dep": new_x,
                        "y_dep": new_y,
                    }
                    inserted.add_particles(**data)
                    fluid.append_parray(inserted)
                    solver.nnps.update()

                    self._pass_n_spawned_since_log += n_ins
                    print(
                        f"C3 inserção t={solver.t:.1f}s: {n_ins} partículas inseridas"
                    )

        if use_wake and solver.count > 0 and solver.count % WAKE_FREQ == 0:
            fluid = self.particles[0]
            m_target = dx * dx

            disp = np.hypot(fluid.x - fluid.x_dep, fluid.y - fluid.y_dep)
            _phi_s = fluid.rho_b_grown
            wake_idx = np.where((_phi_s > WAKE_RHO_B_MIN) & (disp >= WAKE_DISP * dx))[0]

            # Conta so a massa que o WAKE adicionou. Antes comparava a massa TOTAL,
            # entao com a biologia ativa o BiomassGrowth sozinho estouraria o teto e
            # o wake morreria por um motivo alheio a ele.
            budget_ok = (
                self._m_initial is None
                or self._wake_mass_added < WAKE_MASS_BUDGET * self._m_initial
            )
            if not budget_ok:
                wake_idx = np.array([], dtype=int)

            if len(wake_idx) > 0:
                tree = cKDTree(np.column_stack([fluid.x, fluid.y]))
                prox = WAKE_PROX * dx
                prox_sq = prox * prox
                h0 = float(fluid.h[0])
                w_self = float(cubic_spline_w(np.array([0.0]), h0)[0])
                r_ring = WAKE_RING_RATIO * dx
                w_ring = float(cubic_spline_w(np.array([r_ring]), h0)[0])

                new_x = []
                new_y = []
                parent_idx = []
                added_pts = []

                def _livre(px, py):
                    d_ex, _ = tree.query([px, py])
                    if d_ex < prox:
                        return False
                    d_min = d_ex
                    for ax_, ay_ in added_pts:
                        dd = (px - ax_) ** 2 + (py - ay_) ** 2
                        if dd < prox_sq:
                            return False
                        if dd < d_min * d_min:
                            d_min = dd ** 0.5
                    if WAKE_ATTACH > 0.0 and d_min > WAKE_ATTACH * dx:
                        return False
                    return True

                for k in wake_idx:
                    sx = float(fluid.x_dep[k])
                    sy = float(fluid.y_dep[k])
                    ex = float(fluid.x[k])
                    ey = float(fluid.y[k])
                    fluid.x_dep[k] = ex  # reset: deslocamento ja consumido
                    fluid.y_dep[k] = ey

                    if WAKE_SEG:
                        seg = np.hypot(ex - sx, ey - sy)
                        n_seg = min(int(seg / dx), WAKE_SEG_MAX)
                        for j in range(1, n_seg + 1):
                            f_ = j / (n_seg + 1.0)
                            px = sx + f_ * (ex - sx)
                            py = sy + f_ * (ey - sy)
                            if not _livre(px, py):
                                continue
                            new_x.append(px)
                            new_y.append(py)
                            parent_idx.append(int(k))
                            added_pts.append((px, py))

                    if not _livre(sx, sy):  # rastro ja refluido: nao ha vazio
                        continue

                    n_extra = 0
                    if WAKE_MODE == 2:
                        nn = tree.query_ball_point([sx, sy], 2.0 * h0)
                        if len(nn) > 0:
                            nn = np.asarray(nn, dtype=int)
                            d_nn = np.hypot(fluid.x[nn] - sx, fluid.y[nn] - sy)
                            rho_void = float(
                                np.sum(fluid.m[nn] * cubic_spline_w(d_nn, h0))
                            )
                            rho_local = float(np.mean(fluid.rho[nn]))
                            deficit = rho_local - (rho_void + m_target * w_self)
                            if deficit > 0.0:
                                n_extra = int(round(deficit / (m_target * w_ring)))
                            n_extra = max(0, min(n_extra, WAKE_CLUSTER_MAX - 1))

                    if WAKE_SPREAD:
                        n_tot = 1 + n_extra
                        for j in range(n_tot):
                            f_ = j / float(n_tot)
                            px = sx + f_ * (ex - sx)
                            py = sy + f_ * (ey - sy)
                            if not _livre(px, py):
                                continue
                            new_x.append(px)
                            new_y.append(py)
                            parent_idx.append(int(k))
                            added_pts.append((px, py))
                    else:
                        new_x.append(sx)
                        new_y.append(sy)
                        parent_idx.append(int(k))
                        added_pts.append((sx, sy))

                        for j in range(n_extra):
                            ang = 2.0 * np.pi * j / max(n_extra, 1)
                            vx = sx + r_ring * np.cos(ang)
                            vy = sy + r_ring * np.sin(ang)
                            d_ex, _ = tree.query([vx, vy])
                            if d_ex < prox:
                                continue
                            bad = False
                            for ax_, ay_ in added_pts:
                                if (vx - ax_) ** 2 + (vy - ay_) ** 2 < prox_sq:
                                    bad = True
                                    break
                            if bad:
                                continue
                            new_x.append(vx)
                            new_y.append(vy)
                            parent_idx.append(int(k))
                            added_pts.append((vx, vy))

                    if len(new_x) >= WAKE_MAX:
                        break

                if len(new_x) > 0:
                    n_ins = len(new_x)
                    p = np.asarray(parent_idx, dtype=int)
                    if WAKE_RECYCLE:
                        # Rota C: em vez de CRIAR particula (que adiciona massa e
                        # esbarra no teto), REALOCA agar ocioso do campo distante.
                        # Massa exatamente conservada; o buraco deixado fica fora do
                        # disco da colonia, onde C1/C2 nao sao medidos e sigma_a~1.
                        r_all = np.hypot(fluid.x, fluid.y)
                        r_col = float(np.percentile(r_all[fluid.rho_b_grown > 0.1], 99))
                        donor_pool = np.where(
                            (r_all > r_col + WAKE_DONOR_MARGIN)
                            & (fluid.rho_b_grown < 0.05)
                            & (fluid.is_filler < 0.5)
                        )[0]
                        if len(donor_pool) < n_ins:
                            n_ins = len(donor_pool)
                        if n_ins > 0:
                            # amostra ESPALHADA pelo campo distante: tirar sempre os
                            # mais distantes concentraria a depleção nos 4 cantos e
                            # criaria rarefacao (pressao tensil) na fronteira.
                            donors = np.random.choice(
                                donor_pool, size=n_ins, replace=False
                            )
                            p = p[:n_ins]
                            fluid.x[donors] = np.asarray(new_x[:n_ins])
                            fluid.y[donors] = np.asarray(new_y[:n_ins])
                            fluid.u[donors] = 0.0
                            fluid.v[donors] = 0.0
                            fluid.rho_b_grown[donors] = fluid.rho_b_grown[p]
                            fluid.cs[donors] = fluid.cs[p]
                            fluid.c_o[donors] = fluid.c_o[p]
                            fluid.c_n[donors] = fluid.c_n[p]
                            fluid.noise[donors] = fluid.noise[p]
                            fluid.is_filler[donors] = 1.0
                            fluid.is_wake[donors] = 1.0
                            fluid.x_dep[donors] = fluid.x[donors]
                            fluid.y_dep[donors] = fluid.y[donors]
                            solver.nnps.update()
                            self._pass_n_spawned_since_log += n_ins
                            print(
                                f"wake reciclagem t={solver.t:.1f}s: {n_ins} "
                                f"agar realocado (pool={len(donor_pool)})"
                            )
                    else:
                        inserted = fluid.empty_clone()
                        data = {
                            "x": new_x,
                            "y": new_y,
                            "m": [m_target] * n_ins,
                            "h": list(fluid.h[p]),
                            "rho": list(fluid.rho[p]),
                            "rho_b_grown": list(np.maximum(fluid.rho_b_grown[p], FILLER_RHO_B_FLOOR)),
                            "cs": list(fluid.cs[p]),
                            "c_o": list(fluid.c_o[p]),
                            "c_n": list(fluid.c_n[p]),
                            "u": [0.0] * n_ins,
                            "v": [0.0] * n_ins,
                            "noise": list(fluid.noise[p]),
                            "is_filler": [1.0] * n_ins,
                            "is_wake": [1.0] * n_ins,
                            "is_matrix": [1.0 if use_matrix else 0.0] * n_ins,
                            "is_conv": [0.0] * n_ins,
                            "is_env": [0.0] * n_ins,
                            "rho_b_pre": [0.0] * n_ins,
                            "x_dep": new_x,
                            "y_dep": new_y,
                        }
                        inserted.add_particles(**data)
                        fluid.append_parray(inserted)
                        solver.nnps.update()

                        self._wake_mass_added += n_ins * m_target
                        self._pass_n_spawned_since_log += n_ins
                        print(f"wake inserção t={solver.t:.1f}s: {n_ins} inseridas")


if __name__ == "__main__":
    app = SwarmApp()
    app.run()
