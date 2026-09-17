# P2R21 — calota na ponta do grupo + recobrimento atras de p0 (`RASTRO_PONTA_RECOBRE`)

Pre-registrado em 2026-09-17, ANTES de rodar. Base: **P2R20** (P2R16 + `RASTRO_PONTA`=1.0 +
`RASTRO_PONTA_LINK`=6.0). Alavanca unica: `RASTRO_PONTA_RECOBRE` False -> True. t=50, `SEED` fixo,
10 threads.

## A correcao (licao #104-F)
O segmento `p0 -> p1` so converte largura plena entre `p0` e `p1`; atras de `p0` cobre apenas um disco
de raio `w`, e a zona poupada pela calota da chamada anterior (que fica logo atras do novo `p0`) nunca
e reposta. Com `RASTRO_PONTA_RECOBRE`, o lider que recebe calota passa a converter largura plena ate
`a = RASTRO_PONTA*w` atras de `p0` (clip de `t` ate `-a/L`) e o raio de busca de candidatos cresce `a`.

A guarda de baia NAO bloqueia a faixa reconvertida: ela ignora material a <= `w + 0.5 dx` da trajetoria
do proprio lider (ultimas `RASTRO_HIST`=10 posicoes).

## Verificacao antes de rodar (teste-ouro sobre o caminho exato do codigo, passos 1-4 dx)
| | corpo media / min | ponta u=1/2/3/5/7 | desligado == original |
|---|---:|---|---|
| corte plano | 10.00 / 10.00 | 10/10/10/10/10 | bit-a-bit |
| calota | 9.15 / 8.30 | 5.9/7.9/9.1/9.5/9.3 | bit-a-bit |
| **calota + RECOBRE** | **10.00 / 10.00** | **5.9/7.9/9.1/9.9/10.0** | — |

## Predicao
- **Corpo do braco volta a ~9.7 dx** (P2R16 9.7; P2R20 8.3); largura EDT 0.75R ~5.0 (P2R20 4.1).
- **Ponta semicircular**: em fase, u=2 **5.5-8.0**, u=3 **8.0-9.5** (P2R20 3.6 / 6.4; P2R16 7.1 / 9.2).
  A ponta medida fica mais cheia que no P2R20 porque o corpo nao afina mais.
- **Colonizacao viva**: a faixa reconvertida fica exatamente onde estao as co-lideres 3-5 dx atras da
  ponta, entao elas voltam a ser cercadas um passo depois -> lideres externos **19-22**, limbo > 0 apos
  t=32, vivas crescendo.
- Risco: braco de largura plena ate perto da ponta volta a fechar mais baia a meia altura que o P2R20
  (agar 0.5-0.7 R99 em iter 1600: 102 no P2R20, 212 no P2R16). Espero **100-212**.
- R99 volta para perto do P2R16 (4.39 -> 4.5-4.7).

## Criterios (frames em fase, iter 1600-2800; mais t=50)
1. **Objetivo:** corpo u=5 **>= 9.0** dx **e** ponta mais afilada que o P2R16 em u=2 (**< 7.1**).
2. **Colonizacao:** limbo > 0 em todos os frames com t > 32 e vivas crescendo entre t=32 e t=50.
3. Guarda: lideres externos (t >= 20) <= 24; bracos em 0.5R 17 ± 2; R99 >= 4.4; nucleo <= 0.40 Rmax;
   pedacos soltos <= 0.5; agar 0.5-0.7 R99 (iter 1600) >= 100.
4. Motor: `a_mar_bio_p95` (mediana t>30) >= 0.44 (P2R16).
5. Numerica: `max_cs` <= 0.5; iter/t <= 90; massa sem runaway.
6. Visual (§11): pontas arredondadas (nao pontiagudas), braco de largura constante.

## Resultado (2026-09-17, t=50) — REPROVADO: a largura do braco voltou, a ponta continua pontiaguda e o motor caiu

| | P2R16 | P2R20 | **P2R21** | criterio |
|---|---:|---:|---:|---|
| corpo u=5 / u=7 (em fase) | 9.7 / 9.6 | 8.3 / 8.7 | **8.8 / 9.6** | u5 >= 9.0 ❌ |
| largura por braco (EDT, iter 2600) | 4.3 | 3.8 | **4.4** | — |
| largura EDT 0.75R | 5.0 | 4.1 | **5.0** | — |
| ponta u=2 / u=3 | 7.1 / 9.2 | 3.6 / 6.4 | **3.6 / 5.9** | u2 < 7.1 ✅ (previ u3 8.0-9.5 ❌) |
| limbo apos t=32 | 142-145 | 36-100 | **39-163** | > 0 ✅ |
| vivas (iter 1600 -> 2800) | 512 -> 686 | 508 -> 612 | **501 -> 584** | crescendo ✅ |
| lideres externos | 16-20 | 21-22 | **20-21** | <= 24 ✅ |
| bracos em 0.5R (iter 2400/2600/2800) | 16/17/17 | 16/16/19 | **12/14/18** | 17 ± 2 ❌ |
| agar 0.5-0.7 R99 (iter 1600) | 212 (t=28.0) | 102 (t=27.7) | **74 (t=25.2)** | >= 100 ❌ |
| R99 / nucleo / soltos | 4.68 / 0.29 / 0 | 4.39 / 0.32 / 0 | **4.45 / 0.32 / 0** | ✅ |
| `a_mar_bio_p95` (t>30) | 0.44 | 0.79 | **0.37** | >= 0.44 ❌ |
| AR / dedos (lit_rank) | 12.0 / 18 | 12.6 / 19 | **11.4 / 20** | — |
| iter/t / `max_cs` / massa | 57 / 0.494 / 203.2 | 58 / 0.494 / 203.0 | **60 / 0.494 / 202.4** | ✅ |

1. **Largura do braco: restaurada** na medida por EDT (4.4 dx por braco, 5.0 em 0.75R, igual ao P2R16) e em
   u=7 (9.6). O corpo em u=5 (8.8) ainda fica abaixo de 9.0.
2. **Ponta: NAO ficou semicircular.** u=2 3.6 e u=3 5.9, contra 7.9 / 9.1 do teste-ouro. Visualmente
   (`plots/cmp_P2R16_P2R20_P2R21.png`, iter 2600) as pontas seguem pontiagudas, como no P2R20. **O teste-ouro
   descreve so a regiao CONVERTIDA; a ponta real da mascara e mais fina que isso por um motivo que nao
   medi.**
3. Colonizacao viva ✅; lideres 20-21 ✅.
4. **Reprovado:** motor 0.37 (o P2R20 tinha 0.79), bracos em 0.5R 12/14 em dois dos tres frames, agar a meia
   altura 74. Ressalva: a iter 1600 cai em t=25.2 aqui contra t≈28 nos outros (dt diferente); em t=32.7 o
   agar ja e 206.

## O que forma a ponta visivel (2026-09-17, analise offline dos frames)

1. **A ponta visivel e filler do rastro, nao o lider.** Celulas vivas nos ultimos 8 dx atras da extremidade:
   0.0-0.1 por ponta nos tres runs. O lider fica a 0.9-1.3 dx da extremidade, ligeiramente a frente, cercado
   de agar (38-40 particulas a <= 2h) — por isso nao aparece no painel (a), que so marca pontos a > 0.8 dx
   de agar.
2. **Descartados:** (a) discretizacao da rede de particulas — a conversao aplicada a uma rede dx da ponta ate
   MAIS larga (7.2/8.8/9.6); (b) renderizador — bracos sinteticos com a geometria exata passam por
   `FT.campo` com ponta cheia (P2R21 6.5/8.1/8.7) e corpo que bate com a simulacao; (c) movimento do filler —
   ele se desloca 0.5-1 dx para a frente e para FORA do eixo (alarga); (d) guarda de baia — bloquearia so
   2-4% do agar perto do lider.
3. **Sobram 16.7 (P2R16) e 18.7 (P2R21) particulas de AGAR nos 2 dx logo atras de cada lider de ponta**, dentro
   da zona que o corte deveria ter convertido.
4. **Causa: ordem de gravacao.** No `Solver.solve` do PySPH cada iteracao faz `integrator.step` ->
   `post_step` (com `count` ainda antigo) -> `count += 1` -> dump. O frame `main_02600` e gravado ANTES da
   conversao da chamada 2600 e 99 passos DEPOIS da anterior. **Todo frame multiplo de `WAKE_FREQ` mostra o
   rastro no atraso maximo**, com o trecho recem-percorrido ainda como agar. Toda medida de forma de ponta
   feita nesses frames (licoes #103, #104, P2R19-P2R21) mede esse estado.

**Nao estabelecido:** quanto do aspecto pontiagudo e so esse atraso. Separar exige um frame gravado logo
DEPOIS da conversao, que os runs existentes nao tem.

## Diagnostico do atraso de gravacao (2026-09-17) — o atraso ADIANTA a ponta, nao a afina

Run curto do P2R21 ate t=25 (copia do snapshot em scratchpad, sem tocar no repositorio) com uma copia
extra do estado gravada logo DEPOIS de cada conversao (`post_step` em `count % 200 == 0`). Frames
normais **bit-a-bit identicos** ao P2R21 (62 arrays, iteracoes 600 e 1200): a copia extra nao muda a
dinamica.

A regua de perfil por braco (corte em 0.5 R, braco >= 14 dx) nao funciona antes de t≈25: os bracos
ainda nao se separam do corpo e ela mede o disco. Trocada por uma regua centrada em cada lider de ponta
(mais externo do grupo, ligacao 6 dx): extremidade da mascara ao longo do eixo radial do lider e
largura do trecho contiguo que cruza o eixo a `u` dx atras dela. Calibrada nos frames tardios:

| mediana (IQR) | extremidade - lider | u1 | u2 | u3 | u5 | u7 |
|---|---:|---:|---:|---:|---:|---:|
| P2R16, iter 1600-2800 (n=130) | -1.2 dx | 4.4 | 8.0 | 9.3 | 9.6 | 9.6 |
| P2R21, iter 1600-2800 (n=146) | -1.4 dx | 3.3 | 5.9 | 7.2 | 8.8 | 9.6 |
| **P2R21 PRE-conversao, iter 800-1400 (n=74)** | **-2.2 dx** | 5.3 (3.3-11.1) | 9.2 (6.6-13.0) | 12.8 | 16.5 | 20.2 |
| **P2R21 POS-conversao, iter 800-1400 (n=72)** | **-0.4 dx** | 6.1 (5.3-11.3) | 8.5 (7.7-13.3) | 12.9 | 13.8 | 17.6 |

Nos frames tardios o corpo bate com a regua antiga (u5/u7 8.8/9.6), e a ordem P2R16 > P2R21 na ponta se
mantem. Em t≈13-23 as larguras em u5-u7 ainda incluem material vizinho (bracos nao separados); so
u0-u3 descrevem a ponta.

1. **A conversao ADIANTA a ponta 1.8 dx** (a extremidade vai de 2.2 para 0.4 dx atras do lider) e
   limpa o agar atras do lider (mediana 3-4 -> 0-1 particulas a <= 2 dx).
2. **A FORMA da ponta nao muda:** u1/u2/u3 = 5.3/9.2/12.8 antes e 6.1/8.5/12.9 depois, dentro do IQR.
   A ponta anda inteira, ela nao arredonda.
3. **Conclusao:** o frame com atraso maximo desloca a ponta mas nao explica o aspecto pontiagudo. O
   afilamento do P2R21 esta na regiao convertida, nao no atraso de gravacao.

**Ressalva:** so medido em t <= 25, antes de os bracos se separarem. Nos frames tardios, onde a ponta
pontiaguda foi medida (u2 3.6 na regua antiga), o estado pos-conversao nao foi gravado.

### CORRECAO (mesmo dia, run estendido a t=50) — a conclusao acima ESTA ERRADA: a ponta pontiaguda E o atraso

O diagnostico foi repetido ate t=50 (frames normais de novo bit-a-bit identicos ao P2R21 nas iteracoes
600 e 1200). Nos frames em fase onde a ponta foi medida (iter 2200-2800, bracos ja finos e separados),
o estado logo apos a conversao e OUTRO:

| iter 2200-2800, mediana | u1 | u2 | u3 | u5 | u7 |
|---|---:|---:|---:|---:|---:|
| regua por braco — PRE (o que os frames mostram) | 2.1-2.7 | 3.2-4.0 | 5.1-7.0 | 8.4-8.8 | 9.4-9.7 |
| regua por braco — **POS-conversao** | **5.7-6.2** | **7.5-8.0** | **8.7-8.9** | **9.4-9.6** | **9.7-10.5** |
| regua por lider — PRE (n=83) | 2.9 | 4.5 | 6.6 | 8.8 | 9.6 |
| regua por lider — **POS** (n=75) | **5.3** | **7.4** | **8.5** | **9.6** | **9.8** |
| **teste-ouro calota+recobre (geometria pura)** | **5.9** | **7.9** | **9.1** | **9.9** | **10.0** |

**O estado pos-conversao reproduz o teste-ouro.** A calota + recobrimento produzem a ponta semicircular
que a geometria previa; os frames gravados nunca a mostram porque sao gravados no atraso maximo (99
passos depois da conversao anterior). A extremidade tambem sai de 0.9 dx atras do lider para 0.1.

**Isso corrige, alem desta entrada, o item 2 do Resultado acima e as medidas de forma de ponta das
licoes #103 e #104:** todas foram feitas em frames de atraso maximo e descrevem o rastro antes da
reposicao, nao a forma que o mecanismo produz. **O corpo tambem e subestimado:** u5 8.8 (pre) contra
**9.6** (pos), ou seja o criterio 1 pre-registrado (u5 >= 9.0) PASSA no estado pos-conversao.

**Por que a primeira leitura (t <= 25) deu o contrario:** ali os bracos ainda nao se separam do corpo,
a regua mede o disco (larguras de 100-230 dx) e a unica parte valida — u0-u3 — estava saturada pelo
material vizinho. **Regra: antes de concluir a partir de uma regua, checar se ela esta no regime em que
foi calibrada; largura de 200 dx num braco de 5 dx e o sinal de que nao esta.**

Os criterios morfologicos do P2R21 (bracos em 0.5R, agar a meia altura, motor) NAO mudam — sao medidos
sobre a colonia inteira e nao dependem do estado da ponta. A reprovacao do P2R21 se mantem pelo motor
(0.37) e pelos bracos (12/14); o que muda e o diagnostico da PONTA, que nunca foi pontiaguda.
