# P2R4 — rastro largo por SEGMENTO percorrido (sobre o P2R3)

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica contra o P2R3: a regra de
conversao passa de "disco de 2 dx em volta do lider, menos um cone frontal de ±60°"
(`RASTRO_COS_FRENTE`=0.5) para "agar a <= 2 dx do SEGMENTO que o lider percorreu desde a
chamada anterior, e so o que esta atras da posicao atual" (`RASTRO_SEG`=True, `RASTRO_W`=2.0).
t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (medido no P2R3, t≈43)
Os bracos do P2R3 saem cortados no painel (a). O agar que sobra no corredor dos bracos
(30 particulas) fica a 0.73 dx a frente da posicao de chamada (p50), a 0.85 dx do eixo, 47%
dentro da cunha frontal: na chamada k ele esta no cone excluido; na k+1 o lider ja andou
~2.9 dx (p90 3.5) e ele ficou fora do raio de 2 dx. P2R2 (disco cheio): 3 no corredor, mas
sinuoso (tortuosidade 1.22), porque converter a frente abre o caminho do lider.

A posicao anterior e guardada pela propria rotina a cada chamada (indices estaveis). `x_dep`
nao serve: so e reiniciado para quem deposita wake e congela quando o orcamento acaba.
Geometria verificada offline num lider sintetico (24 de 24 convertidas corretas, zero a frente).

## Predicao
- Agar no corredor (mesma medida, t≈43): 30 -> **<= 5**.
- Tortuosidade p50 t 32->43: 1.02 -> 1.00-1.05 (nada a frente do lider e convertido; menos
  vies frontal que o cone de 60° do P2R3). Cos velocidade-radial >= 0.95.
- Largura p50 em r/R99=0.8: 3.2 -> 3.2-3.8 dx; em 0.65: 3.1 -> 3.2-3.8.
- Filler em t≈43: P2R3 4922, P2R4 dentro de ±25%.
- R99, amplitude, dedos: iguais ao P2R3 dentro do ruido.

## Criterios (janela t in [35, 50], contra P2R3 e P2)
1. **Objetivo — continuidade:** agar no corredor <= 8 em t≈43 E pedacos de colonia soltos
   no painel (a) (componentes alem do principal, area >= 4 dx²) <= o do P2R2 no mesmo t.
2. Bracos retos: tortuosidade p50 <= 1.08, mudanca de direcao <= 15°, cos >= 0.93.
3. Largura p50 em r/R99=0.8 >= 2.9 dx (P2R3 3.2 menos o ruido entre realizacoes, ~0.3).
4. Forma: amplitude em R99=2.4 >= 0.8x o P2; sem laços.
5. Numerica: iter/t ate t=50 <= 1.5x o P2, massa sem runaway, `a_press_med` <= 2x o P2.
6. Visual (§11): classe (b) Fingering, bracos continuos no painel (a), baias abertas.

**Base do criterio 1, medida DEPOIS de lancar (texto do criterio inalterado)** — pedacos
>= 4 dx² fora do corpo principal, renderizador final, janela t in [35, 50]: P2 35.0, P2R 25.7,
P2R2 **0.5** (0 em t=44.6), P2R3 14.4 (12 em t=43.4). Largura com o renderizador final
(`fecha`=0): P2R3 3.2 dx em 0.65 e 3.3 em 0.8 (a tabela acima usa os numeros de t≈43).

## Resultado (2026-09-14, t=50, 1526 s) — renderizador final, janela t in [35, 50]

| | P2 | P2R2 | P2R3 | **P2R4** |
|---|---:|---:|---:|---:|
| pedacos soltos >= 4 dx² / > 20 dx² | 35.0 / 8.7 | 0.5 / 0.0 | 14.4 / 10.6 | **0.4 / 0.0** |
| agar no corredor (medida do P2R3, t≈43) | — | 3 | 30 | **13** |
| tortuosidade p50 / rumo / cos velocidade-raio | ~1.02 / 5-10° / 0.98 | 1.22 / 25° / 0.82 | 1.03 / 8° / 0.97 | **1.01 / 8° / 0.97** |
| largura p50 r/R99=0.65 / 0.8 (dx) | 1.6 / 2.1 | 4.3 / 3.9 | 3.2 / 3.3 | **3.2 / 3.4** |
| bracos em 0.65 / 0.8 | 10.0 / 11.3 | 18.0 / 18.8 | 23.4 / 18.0 | 26.4 / 19.6 |
| amplitude (R99=2.4) / dedos / agar limpo | 0.312 / 36 / 78% | 0.163 / 26 / 55% | 0.255 / 31 / 67% | **0.241** / 31 / 65% |
| R99 / Rmin (t=50) / dR/dt | 3.99 / 1.77 / 0.078 | 3.98 / 1.41 / 0.073 | 4.71 / 1.26 / 0.091 | 4.52 / 1.25 / 0.091 |
| C5a / C5b (quorum) | 0.45 / 17.0 | 0.35 / 10.9 | 0.31 / 10.0 | 0.39 / 12.7 |
| filler t≈43 / vivas t=50 | — / 1206 | — / 581 | 4922 / 616 | 4850 / 595 |
| iter/t / massa t=50 / `a_press_med` | 56.3 / 208.0 / — | 69.0 / 202.3 / 2.1e-3 | 63.8 / 203.8 / 1.7e-3 | **57.7 / 204.3 / 1.6e-3** |

1. Continuidade: **pedacos APROVADO** (0.4 contra 0.5 do P2R2; zero acima de 20 dx²).
   **Agar no corredor REPROVADO** — 13 contra o limite de 8 (predicao <= 5), queda de 57% sobre
   o P2R3. O agar restante nao corta os bracos no painel (a). Uma medida refinada (separando o
   que esta a frente da posicao atual) deu numeros ruidosos — ate o P2R2, coberto por
   construcao, aparece com 10 "buracos" — e nao foi usada.
2. Retos: **APROVADO** — 1.01 / 8° / 0.97, o mais reto da serie.
3. Largura: **APROVADO** — 3.4 dx em 0.8 (predicao 3.2-3.8).
4. Forma: sem lacos; **amplitude 0.77x o P2 — REPROVADO no limite** (limite 0.8x; P2R3 0.82x).
5. Numerica: **APROVADO** — iter/t 1.02x o P2, massa abaixo do P2, `a_press_med` como P2R/P2R3.
6. Visual (§11): **APROVADO** — ~22 bracos continuos, retos e de largura uniforme, baias
   abertas, nucleo compacto; classe (b). Figura: `plots/cmp_colonia_P2R4_t50.png`.

Filler em t≈43 dentro da predicao (4850 contra 4922, −1.5%). Conversoes por chamada 89-185.

Figuras: `plots/cmp_colonia_P2R4_t50.png` (P2R2/P2R3/P2R4), `plots/fig_P2R4_t50.png` (figura de tese, `--contorno`).
