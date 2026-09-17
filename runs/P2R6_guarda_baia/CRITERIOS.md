# P2R6 — rastro por segmento W = 3 dx + guarda de baia (sobre o P2R5)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R5: `RASTRO_BAIA`
0 -> 2.0 dx (`RASTRO_HIST`=10). t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (P2R5, t=50; regioes marcadas pela usuaria + viewer do PySPH)
Com W = 3 a base dos bracos se fundiu num corpo largo ate r ≈ 2 (0.4 R99); as baias da base
viraram frestas de agar de 1-2 dx (B, D) e bolsas fechadas (E). O agar de B e D foi convertido
em t = 18-33 (mediana 24.5), quando aquele raio era perto da frente e os bracos estavam a ~9 dx
um do outro: as faixas de ±3 dx de dois vizinhos se sobrepoem. Nas frestas A e C, lideres
andando dentro da baia deixam rastro tracejado do wake (abaixo de 0.7 R99, sem conversao).

## Regra
Agar a <= W do segmento percorrido so e convertido se nao houver material nao-agar VIZINHO a
<= 2 dx dele. Vizinho = qualquer particula nao-agar a mais de W + 0.5 dx da trilha do proprio
lider (ultimas 10 posicoes de chamada + posicao atual). Teste sintetico: vizinho a 4 dx -> a
faixa para em 1 dx e sobram 2 fileiras de agar; o proprio braco atras nao bloqueia.

## Predicao
- Baias da base (r/R99 = 0.4 e 0.5): frestas < 2 dx caem para menos da metade do P2R5; baia p50 sobe.
- Bolsas de agar fechadas no corpo: menos que no P2R5.
- Largura p50 em 0.8 R99: 5.5 -> 5.2-5.6 (la as baias tem 10-18 dx, a guarda quase nao age).
- Amplitude (R99=2.4): 0.168 -> 0.18-0.21; agar limpo nas baias 42% -> 45-55%.
- Filler t≈41: 7564 -> 10-20% menos. Bracos: numero igual ±2.

## Criterios (janela t in [35, 50], contra P2R5)
1. **Objetivo:** numero de baias < 2 dx em r/R99 = 0.4 e 0.5 (soma, media na janela) <= 0.5x
   o do P2R5, E bolsas de agar fechadas no painel (a) <= 0.5x o do P2R5.
2. Manter a largura: p50 em 0.8 R99 >= 4.5 dx.
3. Retos: tortuosidade p50 <= 1.08. Continuidade: pedacos soltos <= 2 por quadro.
4. `Rmin`(t=50) <= 2.0. Numerica: iter/t <= 1.5x o P2, massa sem runaway.
5. Visual (§11 + viewer): baias da base abertas (nao frestas), classe (b).

**Emenda ao criterio 1, feita com o run em andamento e ANTES de ver qualquer resultado dele**
(base medida nos runs existentes, [baias_base2.py]): o numero de frestas < 2 dx cai tambem
quando a baia FECHA DE VEZ (baia fechada nao e contada) — o P2R5 tem menos frestas que o P2R4
(7.2 contra 12.0) justamente por isso. Contrapeso: **baias abertas (>= 2 dx) em r/R99 = 0.4 e
0.5 (soma, media na janela) >= a do P2R5**. Bases: frestas < 2 dx — P2R4 12.0, P2R5 7.2
(alvo <= 3.6); baias abertas — P2R4 36.2, P2R5 30.3 (alvo >= 30.3); bolsas fechadas >= 4 dx² —
P2R4 1.6, P2R5 1.0 (alvo <= 0.5). P2: 2.8 / 4.2 / 0.0 (em 0.4-0.5 R99 o P2 ainda e o disco).

## Resultado (2026-09-15, t=50, 1551 s) — renderizador final, janela t in [35, 50]

| | P2R4 | P2R5 | **P2R6** |
|---|---:|---:|---:|
| frestas < 2 dx em 0.4+0.5 R99 | 12.0 | 7.2 | **8.0** |
| baias abertas >= 2 dx em 0.4+0.5 R99 | 36.2 | 30.3 | **36.7** |
| bolsas de agar fechadas >= 4 dx² | 1.6 | 1.0 | **1.8** |
| largura p50 0.65 / 0.8 R99 (dx) | 3.2 / 3.4 | 5.5 / 5.5 | 5.4 / 5.5 |
| bracos em 0.65 / 0.8 | 26.4 / 19.6 | 21.8 / 18.0 | 23.8 / 18.2 |
| pedacos soltos / tortuosidade | 0.4 / 1.01 | 1.7 / 1.01 | 1.3 / 1.01 |
| amplitude / agar limpo / Rmin | 0.241 / 65% / 1.25 | 0.168 / 42% / 1.17 | 0.156 / 46% / 1.28 |
| filler t≈41 / iter/t / massa | 4850 / 57.7 / 204.3 | 7564 / 60.3 / 204.3 | 7299 / 60.4 / 204.2 |

1. Objetivo: **REPROVADO** — baias abertas +21% (aprovado), mas frestas 8.0 (alvo <= 3.6) e
   bolsas 1.8 (alvo <= 0.5).
2-4. **APROVADOS** (largura 5.5, tortuosidade 1.01, 1.3 pedacos, Rmin 1.28, numerica igual).
5. Visual: baias da base ainda frestas — **REPROVADO**.
Predicoes: agar limpo 46% (ok); amplitude 0.156 (errou — previ 0.18-0.21); filler −3.5% (previ −10-20%).

**Causa (medida):** a borda das baias da base e 88% agar convertido pelo rastro, igual ao P2R5.
Dos pares convertidos em lados opostos de uma fresta (< 3 dx, base), 46% foram convertidos no
MESMO quadro (P2R5: 35%) e so 17% com >= 2 quadros de diferenca (P2R5: 27%). Dois defeitos:
(a) **falha de implementacao** — a arvore de material vizinho e montada no inicio da chamada,
entao dois lideres vizinhos que convertem na mesma chamada nao se enxergam; (b) **G = 2 dx deixa
uma fresta por construcao** — a guarda preserva ~2 fileiras de agar, que no painel (a) viram uma
fenda de 1.6-2.6 dx. Proxima rodada: incluir as conversoes da propria chamada na guarda e G = 4 dx.

Figura do diagnostico (classe das particulas nas regioes A-E marcadas pela usuaria no P2R5): `plots/zoom_regioes_P2R5.png`.
