# P2R18 — líder interno por VELOCIDADE (`RASTRO_V_MOV` = 0.02) no lugar do passo por chamada

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R16** (verificado: o `main.py` difere do
snapshot do P2R16 por UMA constante). Alavanca unica: `RASTRO_V_MOV` 0 -> 0.02. t=50, `SEED` fixo,
10 threads.

## O defeito, medido no proprio P2R16

O criterio atual de lider interno e `passo >= RASTRO_PASSO_MOV*dx` = **0.5 dx por chamada**, com a
chamada a cada `WAKE_FREQ`=100 iteracoes. Isso e um criterio de VELOCIDADE com limiar variavel:
`v_exigida = 0.5*dx / (100*dt)`. Medido no log do P2R16:

| janela | `dt` (p50) | `dt_call` | **`v` exigida** |
|---|---:|---:|---:|
| t 0-10 | 0.01583 | 1.58 s | 0.0170 |
| t 10-20 | 0.01811 | 1.81 s | 0.0154 |
| t 20-30 | 0.01951 | 1.95 s | **0.0138** |
| t 30-40 | 0.01924 | 1.92 s | 0.0140 |
| t 40-50 | 0.01692 | 1.69 s | 0.0159 |

Sobre o run: **min 0.0120, p50 0.0147, max 0.0263 — razao 2.2x**. O limiar oscila por um fator 2.2
sem que nada na fisica mude; e o `dt` adaptativo entrando no criterio. As rajadas de lideres saem
dai.

**A separacao entre as populacoes e limpa e justifica o corte:** velocidade radial das vivas,

| | frente (`r >= 0.7 R99`) | candidatos internos (`0.4 <= r/R99 < 0.7`) |
|---|---:|---:|
| `v_r` p50 | **0.057 a 0.072** | **0.0003 a 0.0069** |

uma ordem de grandeza. Com limiar 0.02 passam 0-12 por frame (p50 ~5); com o limiar efetivo atual
(0.0147) passam ~9.

## Metrica-alvo: onde os bracos NASCEM

Contando bracos por raio (componentes da mascara, EDT, `tools/ciclo_largura.py`), o P2R16 da
**11 (0.4R) -> 17 (0.5R) -> 19 (0.6R) -> 13 (0.7R)**. O acrescimo **+8 entre 0.4R e 0.6R** sao os
bracos que nascem no MEIO — exatamente o que a alavanca ataca.

## Predicao

- Lideres internos aprovados por chamada: p50 ~9 -> **~5** (-40%), **e a oscilacao de 2.2x some**.
- **Acrescimo de bracos 0.4R -> 0.6R: +8 -> +3 a +5.**
- Bracos em 0.5R: 17 -> **13-15**.
- **Largura por braco: 4.4 dx, INALTERADA** — a alavanca decide QUEM lidera, nao a geometria do
  rastro. E o teste de falseabilidade do mecanismo.
- Ciclo em 0.6R: 0.60 -> **0.48-0.55** (menos bracos, mesma largura — cai por construcao).
- `R99`/`Rmax` 4.68/4.97 -> iguais ou ate +5% (menos competicao por material na frente).
- Motor `a_mar_bio_p95` 0.44 -> **0.5-0.7** (menos agar convertido por lider interno deixa mais
  recrutavel).
- Nucleo 0.32 R inalterado (a alavanca nao toca a base).
- iter/t ~57, inalterado.

## Criterios (t=50 e janela t in [35, 50])

1. **Objetivo:** acrescimo de bracos 0.4R -> 0.6R **<= +5** (hoje +8) **e** bracos em 0.5R <= 15.
2. **Falseabilidade do mecanismo:** largura por braco em 4.4 ± 0.7 dx. **Se sair dessa faixa, o
   mecanismo NAO e o que suponho** (ele nao deveria tocar a largura) e a leitura do resultado muda.
3. **Guarda:** pedacos soltos <= 0.5; `R99` >= 4.4; nucleo <= 0.40 Rmax.
4. Campo/motor: `a_mar_bio_p95` (mediana t>30) >= 0.44 (nao piorar); `c_s` do corpo >= 0.15.
5. Numerica: `max_cs` <= 0.5; iter/t <= 90; massa sem runaway.
6. Visual (§11): classe (b); bracos mais uniformes em comprimento, baias atravessando da ponta ate
   perto do nucleo.

## Ambiguidade declarada

`RASTRO_V_MOV` = 0.02 faz **duas** coisas ao mesmo tempo: (a) remove a oscilacao de 2.2x e (b) aperta
36% sobre o limiar mediano atual (0.0147). **Se reprovar, nao sabera qual das duas causou.** A versao
isolante seria 0.015 (equivalente a mediana atual, so removendo a oscilacao) — mas ela nao reduziria
a contagem media de lideres, que e o objetivo declarado. Escolhi 0.02 por isso; registro a limitacao.

## Reguas usadas (todas persistidas antes da rodada)

`tools/ciclo_largura.py` (EDT sobre a mascara rasterizada — a v1 por fatia fina de particulas estava
errada, licao #101), `tools/rank_runs.py`, `tools/lit_rank.py`, `plots/fig_tese.py`.

## Resultado (2026-09-16) — NO-OP: bit-a-bit identico ao P2R16. Nenhum criterio e avaliavel.

Estado final (t=50, iteracao 2878): `x`, `y`, `rho_b_grown`, `is_filler`, `cs`, `u`, `v`, `m` com
**max|dif| = 0.000e+00 em 70 106 particulas**, e o `log.csv` identico linha a linha.

**NAO e erro de configuracao.** Verificado por instrumentacao (print no ramo do criterio, run curto
a t=12): o ramo de VELOCIDADE executa e **seleciona conjuntos diferentes** do criterio de passo —

| t | `dt_call` | limiar efetivo | candidatos | passam por VELOCIDADE | passariam por PASSO |
|---|---:|---:|---:|---:|---:|
| 4.49 | 2.195 | 0.82 dx | 7 | 3 | 3 |
| 5.19 | 0.706 | 0.26 dx | 7 | **4** | **2** |
| 6.53 | 1.341 | 0.50 dx | 8 | 5 | 5 |
| 8.13 | 1.597 | 0.59 dx | 19 | **6** | **8** |
| 9.70 | 1.569 | 0.58 dx | 33 | 6 | 6 |
| 11.79 | 2.095 | 0.78 dx | 62 | **7** | **17** |

(a tabela tambem confirma o defeito que motivou a alavanca: o limiar efetivo do criterio antigo
varia de 0.26 a 0.82 dx entre chamadas.)

**Em t=11.79 o criterio antigo elegeria 10 lideres a mais, e o estado nao mudou em nada.**

**Teste decisivo:** rodada com `RASTRO_R_MOV` = 0 (lideres internos DESLIGADOS por completo) —
tambem **bit-a-bit identica** ao P2R16 nas iteracoes 200/400/600, com a mesma contagem de filler
(146 / 478 / 731).

**Conclusao: na base P2R16 os lideres internos sao INERTES** — nao convertem nenhum agar que os
lideres externos ja nao convertam. Duas causas medidas: (a) com `RASTRO_W`=5 o rastro dos externos
cobre 10 dx de largura e a banda 0.4-0.7 R99 fica dentro do que eles ja converteram; (b) a guarda de
baia (`RASTRO_BAIA`=2.0) barra **100% do agar** vizinho aos candidatos internos ate t≈21 (e 47-52%
depois). O mecanismo validado no P2R8 (licao #94-F, onde `RASTRO_W` era 3 e resolveu os dendritos
tracejados) foi **neutralizado por duas alavancas adicionadas depois**.

**Consequencia:** o objetivo (reduzir os bracos que nascem em 0.4-0.7 R99) **nao e atingivel por esta
alavanca nesta base** — os bracos novos dali nao vem da conversao de agar pelos lideres internos.
Para atacar o defeito e preciso primeiro identificar o que de fato os cria (candidatos: o wake, cuja
mae precisa so de `rho_b` >= `WAKE_RHO_B_MIN`=0.05, ou a colonizacao).
