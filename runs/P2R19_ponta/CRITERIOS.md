# P2R19 — arredondar a PONTA do braço (`RASTRO_PONTA`): calota elíptica atrás de p1

Pre-registrado em 2026-09-16, ANTES de lancar (o P2R18 ainda estava rodando). Alavanca unica:
`RASTRO_PONTA` 0 -> **1.0**. **Base: P2R16** (o P2R18 saiu NO-OP bit-a-bit, licao #102). `RASTRO_V_MOV`
fica em 0.0 — inerte nesta base, mas voltaria a agir se `RASTRO_W` mudar. t=50, `SEED` fixo, 10 threads.

## A causa, no codigo

[main.py](../../main.py) `_rastro_segmento` selecionava o agar a converter com

    sel = cand[(s <= 1.0) & (dist <= w)]

`s` e a projecao normalizada ao longo do segmento `p0 -> p1` (posicao anterior -> atual do lider) e
`dist` a distancia ao segmento. Atras de `p0` o corte e um semicirculo de raio `w` (porque `t` e
clipado em 0); **na frente, `s <= 1.0` e um corte PLANO perpendicular ao segmento, de largura
`2w` = 10 dx** com `RASTRO_W`=5. A ponta do braco e literalmente a secao transversal do rastro.

## O defeito, medido

Perfil de largura recuando `u` a partir da extremidade do braco (mediana sobre os bracos que chegam
alem de 0.75 Rmax; mascara do painel (a), t=50):

| `u` (dx) | 0 | 1 | 2 | 3 | 4 | 5 | 7 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P2R15 (W=3) | 2.3 | 4.4 | 5.2 | 5.5 | 5.4 | 5.4 | 5.7 | 5.7 |
| **P2R16 (W=5)** | 2.5 | **6.6** | **9.0** | 9.3 | 9.5 | 9.7 | 9.6 | 9.8 |
| P2R17 (W=5 + conduz) | 1.9 | 7.0 | 8.8 | 9.7 | 9.9 | 9.7 | 9.9 | 9.8 |

**A largura plena e atingida em 2 dx de recuo.** Uma calota hemisferica de raio `w`=4.9 dx daria
5.9 / 8.0 / 9.2 em u = 1 / 2 / 3 — ou seja **a ponta de hoje e ~2x mais abrupta que hemisferica**.

## A correcao

Calota eliptica de semi-eixos `a = RASTRO_PONTA*w` (longitudinal) e `w` (transversal), atras de p1:

    u   = clip((1 - s)*|p1-p0|, 0, a)          # recuo longitudinal atras da ponta
    lim = w*sqrt(u*(2a - u))/a                 # largura permitida; lim(0)=0, lim(a)=w
    sel = cand[(s <= 1.0) & (dist <= lim)]

`RASTRO_PONTA` = 1.0 da `a = w`, ou seja **calota semicircular**. `RASTRO_PONTA` = 0 cai no ramo
`else` com `lim = w`, o codigo original.

**Por que NAO remover o `s <= 1.0` (a alternativa obvia):** isso converteria agar A FRENTE do lider,
que e exatamente o que fez o P2R2 serpentear (licao #94-A) — a ponta avanca sobre filler que ela
mesma criou, que nao entra no gradiente de `cs`, e perde o rumo. A calota afila por tras sem tocar
no que esta a frente.

## Verificacoes feitas antes de pre-registrar

1. **Teste-ouro da geometria** (grade de 801², segmento de 3 dx): com `RASTRO_PONTA`=1.0 a largura
   medida da regiao convertida e 0.0 / 5.9 / 8.0 / 9.2 em u = 0 / 1 / 2 / 3 dx, contra o semicirculo
   ideal 0.0 / 6.0 / 8.0 / 9.2 — **erro <= 0.3 dx**.
2. **`RASTRO_PONTA` = 0 reproduz o original bit-a-bit** (`np.array_equal` da mascara de selecao).
3. Sem `NaN`/warning: o `u` e clipado em `[0, a]` antes da raiz (a primeira versao usava `np.where`,
   que avalia os dois ramos e fazia `sqrt` de negativo para `u > 2a`).

## Predicao

- **A ponta passa a seguir o semicirculo:** u=1: 6.6 -> **5.5-6.2**; u=2: 9.0 -> **7.5-8.3**;
  u=3: 9.3 -> 8.8-9.4. Razao u1/u5: 0.68 -> **0.55-0.62**.
- **A largura do CORPO do braco (u >= 5 dx) nao muda: 9.7 dx.** O afilamento e transitorio — um ponto
  a `u` atras da ponta atual recebe `lim(u + j*p)` nas chamadas seguintes (`p` = passo por chamada,
  2.9 dx), e com `a`=5 dx basta **j=2 chamadas** para ele receber largura plena. So a ponta corrente
  fica afilada.
- Area convertida por chamada: -16% (teste-ouro), mas **sobrescrita nas duas chamadas seguintes** —
  o efeito liquido em massa/area total e desprezivel (~5 dx² por braco, so na ponta final).
- `R99`, `Rmax`, numero de bracos, nucleo, motor e campo: **inalterados**. O alcance e do LIDER, nao
  da conversao.

## Criterios (t=50 e janela t in [35, 50])

1. **Objetivo:** largura em u=2 dx **<= 8.3** (hoje 9.0) **e** razao u1/u5 **<= 0.62** (hoje 0.68).
2. **Falseabilidade:** largura do corpo (u >= 5 dx) em 9.7 ± 0.7 dx. **Se o corpo afinar junto, o
   afilamento NAO e transitorio como calculei** e a leitura muda.
3. **Guarda (nada mais pode mexer):** `R99` >= 4.4; bracos em 0.5R dentro de ±2 da base; nucleo
   <= 0.40 Rmax; pedacos soltos <= 0.5.
4. Motor/campo: `a_mar_bio_p95` >= o da base (nao piorar); `c_s` do corpo >= o da base.
5. Numerica: `max_cs` <= 0.5; iter/t <= 90; massa sem runaway.
6. Visual (§11): pontas visivelmente arredondadas, corpo do braco igual.

## Risco identificado

Se o lider MUDAR DE DIRECAO dentro da calota, o entalhe lateral que ela deixa pode nao ser preenchido
pela chamada seguinte e virar uma mordida no contorno. Medido nas rodadas da serie: tortuosidade
1.00-1.01 e desvio de direcao ~4 graus por chamada (licao #94-D/F) — **risco baixo**, mas e o que o
criterio 3 (pedacos soltos) e o 6 (visual) vigiam.

## Teste de efeito ANTES da rodada (regra da licao #102)

Run curto (t=12) com `RASTRO_PONTA`=1.0 contra o P2R16, mesmas iteracoes:

| iteracao | t | particulas | **filler** | max\|dx\| |
|---|---:|---:|---:|---:|
| 200 | 4.5 | 68267 / 68267 | 146 / 146 | 0 |
| 400 | 6.5 | 68380 / 68441 | **478 / 429 (-10%)** | 3.6e-1 |
| 600 | 9.7 | 68430 / 68483 | **731 / 915 (+25%)** | 1.28 |

**A alavanca TEM efeito** (diverge a partir da iteracao 400; a 200 e identica porque antes da primeira
conversao do rastro nao ha o que diferir). O `max|dx|` de 1.28 e assinatura de caos, nao de efeito
medio.

**DISCREPANCIA COM A PREDICAO, registrada antes de rodar:** previ area convertida liquida
desprezivel (-16% por chamada, sobrescrita nas duas seguintes). Em t=6.5 bate (-10%); **em t=9.7 o
filler esta +25%**, com +53 particulas. Pode ser flutuacao de uma realizacao ja caotica — ou a calota
esta mudando a trajetoria dos lideres de um jeito que converte MAIS. **O criterio 2 (largura do corpo
9.7 +- 0.7 dx) e o que distingue os dois**: se o filler extra for engrossamento do braco, ele reprova.

## Resultado (2026-09-16, t=50) — REPROVADO: a ponta afila demais, o braco inteiro afina, e a colonizacao MORRE

### Correcao da propria linha de base (antes do veredito)

O valor-base pre-registrado (largura em u=2 = 9.0 dx) veio do frame FINAL do P2R16, iteracao 2878,
que esta **fora de fase** com o ciclo de conversao do rastro (78 iteracoes depois de uma chamada; o do
P2R19 e a 2811, 11 depois). Nos frames multiplos de 200 (mesma fase nos dois runs) o mesmo P2R16 da
u=2 entre 2.5 e 9.0. **O perfil de ponta de um frame isolado e funcao da fase, nao so do run.** A
comparacao abaixo usa so frames em fase, iter 1600-2800, agregando bracos e frames (n ~56 amostras).

| em fase, iter 1600-2800 | u=1 | **u=2** | **u=3** | **corpo u=5** | u=7 | bracos 0.5R (t=46 / 50) |
|---|---:|---:|---:|---:|---:|---:|
| P2R16 | 2.8 | **7.1** [3.9-8.8] | 9.2 | **9.7** [9.4-10.0] | 9.6 | 17 / 17 |
| P2R19 | 2.2 | **3.3** [2.8-4.3] | **4.8** | **8.0** [7.4-8.8] | 8.6 | 19 / 23 |
| semicirculo previsto | 6.0 | 8.0 | 9.2 | 9.7 | 9.7 | — |

1. **Objetivo: atingido e ULTRAPASSADO.** u=2 7.1 -> 3.3 (limite 8.3). Mas o semicirculo previa 8.0 em
   u=2 e 9.2 em u=3; medido 3.3 e 4.8. **A ponta ficou PONTIAGUDA, nao arredondada** (figura
   `plots/cmp_P2R16_P2R19.png`).
2. **Falseabilidade: REPROVADO.** Corpo em u=5 9.7 -> **8.0** (faixa 9.0-10.4); largura por braco por
   EDT 4.3 -> 3.8 dx. **O afilamento NAO foi transitorio** — os bracos viraram fusos, afinando ao
   longo de todo o comprimento. Minha conta ("basta j=2 chamadas para a faixa receber largura
   plena") esta errada em algum ponto que nao identifiquei.
3. Guarda: R99 4.68 -> 4.57 ✅, nucleo 0.29 -> 0.32 R ✅, pedacos soltos 0 ✅ (ha fragmentos pequenos
   perto da base, abaixo do limiar da regua); **bracos em 0.5R 17 -> 23 em t=50 ❌** (limite ±2).
4. **Motor: REPROVADO**, `a_mar_bio_p95` 0.44 -> 0.35. Campo do corpo 0.19 -> 0.23 ✅.
5. Numerica: ✅ — iter/t 56, `max_cs` 0.494, massa 203.6.
6. Visual: pontas pontiagudas, bracos em fuso, nao arredondados.

### Achado inesperado: a COLONIZACAO PARA em t≈32

| t≈ | 25 | 28 | 32 | 36 | 39 | 43 | 46 | 50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| vivas P2R16 | 503 | 512 | 519 | 554 | 593 | 619 | 667 | 696 |
| **vivas P2R19** | 506 | 523 | **525** | **525** | **525** | **525** | **525** | **525** |
| limbo P2R16 | 64 | 91 | 135 | 143 | 138 | 145 | 134 | 146 |
| **limbo P2R19** | 19 | 2 | **0** | **0** | **0** | **0** | **0** | **0** |

O limbo vai a zero e FICA em zero: **nenhuma particula de agar e colonizada nos ultimos 18 s**. No P2R16
nascem 30-50 limbos novos por frame ate o fim. A simulacao segue evoluindo (R99, filler e biomassa
crescem), so o recrutamento morre. **E o que explica o motor mais fraco.**

### Mecanismo: NAO ESTABELECIDO

Testei que a guarda de baia seria a causa do afilamento nao-transitorio (a faixa lateral deixada pela
calota fica a <= 2 dx do filler do eixo e nao poderia mais ser convertida). **Refutado:** a guarda
bloqueia fracao parecida do agar dos bracos nos dois runs (43-49% no P2R16, 49-55% no P2R19). Nao tenho
mecanismo medido nem para o afilamento continuo nem para a morte da colonizacao, e nao afirmo um.

## Diagnostico pos-reprovacao (2026-09-17) — onde a colonizacao morre

Runs instrumentados ate t=36 (so prints; **bit-a-bit identicos** ao P2R16 e ao P2R19 nas iteracoes
600/1200/1800). Por chamada do rastro: lideres externos/internos, agar convertido e conversao
EXCLUSIVA dos internos (externos vem primeiro e o `setdefault` fica com o primeiro).

| | P2R16 (`PONTA`=0) | P2R19 (`PONTA`=1) |
|---|---:|---:|
| agar convertido ate t=36 | 6626 | 6973 |
| conversao exclusiva dos lideres internos | 58 (1%) | **0** |
| convertido abaixo de 0.7 R99 | 101 (2%) | 4 |
| **lideres externos (t >= 16)** | **19-20** | **25-27** |

**Hipotese 1 REFUTADA** (lideres internos convertendo os bolsoes do meio): eles sao inertes nos dois
runs, e o rastro quase nao converte abaixo de 0.7 R99 em nenhum.

**O que mudou foi o numero de pontas ativas (+30%)**, divergindo entre t≈10 e t≈17. Consequencia medida
na faixa 0.5-0.7 R99 (frames em fase):

| iteracao | 1000 | 1200 | 1400 | 1600 |
|---|---:|---:|---:|---:|
| agar a <= 2h de material, P2R16 | 28 | 88 | 143 | 212 |
| agar a <= 2h de material, P2R19 | **0** | **0** | **1** | **21** |
| `cs` p99 desse agar, P2R16 | 0.274 | 0.284 | 0.292 | 0.293 |

**No P2R16 a colonizacao vive nesses bolsoes do meio**, encostados em 50-70 celulas vivas e no filler,
onde o `cs` do agar sobe ate 0.29-0.30 — **encostado no gate `COL_CS_MIN`=0.3**. No P2R19 as baias ja
fecham enquanto a faixa ainda e frente (mais rastros vizinhos se encostam); sobra so agar junto a
frente, com `cs` p99 <= 0.25, e o gate nunca mais abre.

**Hipotese 2 REFUTADA** (lider com corte plano fica cercado do proprio filler, desacelera e sai da
frente): com a calota as celulas da frente tem mais agar vizinho (40% contra 15% em t≈7), mas a
velocidade radial e parecida, e o P2R16 ate e mais rapido apos t≈17. **Por que a calota mantem mais
pontas ativas continua sem mecanismo medido.**

## Mecanismo do +30% de pontas (2026-09-17, rastreio de identidade dos lideres nos frames)

Lider nao vira filler nem perde `rho_b`, e so 0-3 foram pinados; quem deixa de ser lider PARA (avanca
0-18% do crescimento do R99). A partir de t≈17 o conjunto e fechado (nenhuma entrada). A diferenca
inteira esta nas saidas entre t≈10 e 17: **P2R16 perde 33 de 53, P2R19 perde 18 de 47**.

| t≈10 -> 17 | no grupo da ponta | atras da mais externa | agar vizinho | `rho_b` | `v_r` |
|---|---|---:|---:|---:|---:|
| SAEM, P2R16 (33) | nunca sozinha, nunca a mais externa | **4.6 dx** | **0%** | 0.19 | 0.028 |
| FICAM, P2R16 (20) | mais externa | 0.0 dx | 67% | 0.35 | 0.069 |
| SAEM, P2R19 (18) | nunca sozinha, nunca a mais externa | **7.0 dx** | **0%** | 0.21 | 0.010 |
| FICAM, P2R19 (29) | mais externa ou logo atras | 1.3 dx | 75% | 0.37 | 0.063 |

**Quem para e sempre co-lider atras da ponta do proprio grupo, sem agar em volta.** Hipotese de perda
do flagelo por `rho_b` > 0.6 refutada (nenhuma cruza).

**A calota controla exatamente isso.** Com corte plano todo agar a <= 5 dx do segmento vira filler ate
a posicao do lider, e uma co-lider 4.6 dx atras fica cercada. A calota poupa o agar lateral nos
primeiros `a` = 5 dx atras da ponta; co-lideres nessa zona mantem contato com agar e seguem; so param
as que estao alem de 5 dx (medido: 7.0 dx).

**Cadeia completa, toda medida:** calota poupa agar lateral a < `a` atras da ponta -> co-lideres
nao sao cercadas -> +30% de pontas -> baias a meia altura fecham ainda na frente -> somem os bolsoes
de agar onde o `cs` chega ao gate -> colonizacao morre.
