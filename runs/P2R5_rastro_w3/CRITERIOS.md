# P2R5 — rastro por segmento com W = 3 dx (sobre o P2R4)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R4: `RASTRO_W`
2.0 -> 3.0 (`RASTRO_SEG`=True). t=50, `SEED` fixo, 10 threads.

## Objetivo (pedido da usuaria)
Engrossar ainda mais os bracos; aceitavel reduzir o numero de bracos.

## Estimativa estatica (P2R4 em t=50, colonia do painel (a) dilatada 1 dx por lado)
| raio | bracos | baia p50 | +1 dx |
|---|---:|---:|---:|
| 0.4 R99 | 24 | 3.6 dx | 19 |
| 0.5 R99 | 26 | 5.7 dx | 23 |
| 0.65 R99 | 24 | 9.9 dx | 26 |
| 0.8 R99 | 20 | 17.9 dx | 20 |
Engrossa sem fundir no meio e fora; fecha baias so na base (nucleo cresce, licao #90).

## Predicao
- Largura p50 em r/R99 = 0.65 / 0.8: 3.2 / 3.4 -> **5.0-5.6 dx**.
- Bracos em 0.65 / 0.8: ~20-26 / 18-20 (sem fusao); em 0.4 R99: 24 -> ~19.
- `Rmin` em t=50: 1.25 -> 1.4-1.7. Amplitude (R99=2.4): 0.241 -> ~0.20.
- Filler em t≈42: 4850 -> ~6500-7500 (faixa convertida 4 -> 6 dx). Vivas ~600.
- Tortuosidade inalterada (~1.01); numerica inalterada.
- Reducao dinamica do numero de bracos (filler ao lado das pontas vizinhas): possivel, nao prevista.

## Criterios (janela t in [35, 50], contra P2R4 e P2)
1. **Objetivo:** largura p50 em r/R99=0.8 >= 4.5 dx.
2. Continuidade: pedacos soltos >= 4 dx² <= 1 por quadro (media).
3. Retos: tortuosidade p50 <= 1.08.
4. Nao virar disco: `Rmin`(t=50) <= 2.0; amplitude >= 0.6x o P2 (0.187).
5. Numerica: iter/t <= 1.5x o P2, massa sem runaway, `a_press_med` <= 2x o P2R4.
6. Visual (§11): classe (b) Fingering, baias abertas.

## Resultado (2026-09-15, t=50, 1670 s) — renderizador final, janela t in [35, 50]

| | P2 | P2R4 (W=2) | **P2R5 (W=3)** |
|---|---:|---:|---:|
| largura p50 r/R99=0.65 / 0.8 (dx) | 1.6 / 2.1 | 3.2 / 3.4 | **5.5 / 5.5** |
| bracos em 0.65 / 0.8 | 10.0 / 11.3 | 26.4 / 19.6 | 21.8 / 18.0 |
| pedacos soltos >= 4 dx² / > 20 dx² | 35.0 / 8.7 | 0.4 / 0.0 | **1.7** / 0.0 |
| tortuosidade p50 / rumo / cos | ~1.02 / 5-10° / 0.98 | 1.01 / 8° / 0.97 | **1.01 / 4° / 0.98** |
| amplitude (R99=2.4) / dedos / agar limpo | 0.312 / 36 / 78% | 0.241 / 31 / 65% | **0.168** / 30 / 42% |
| R99 / Rmin (t=50) / dR/dt | 3.99 / 1.77 / 0.078 | 4.52 / 1.25 / 0.091 | 4.87 / 1.17 / 0.098 |
| filler t≈41 / vivas t=50 | — / 1206 | 4850 / 595 | 7564 / 598 |
| iter/t / massa / `a_press_med` | 56.3 / 208.0 / — | 57.7 / 204.3 / 1.6e-3 | 60.3 / 204.3 / 1.4e-3 |

1. Objetivo (largura >= 4.5 dx em 0.8): **APROVADO** — 5.5 dx, dentro da predicao (5.0-5.6).
2. Continuidade (<= 1 pedaco por quadro): **REPROVADO** — 1.7 (0, 0, 2, 2, 3, 3), todos
   pequenos (nenhum > 20 dx²), fragmentos perto do nucleo.
3. Retos: **APROVADO** — 1.01 / 4° / 0.98.
4. Disco: `Rmin` 1.17, estavel desde t≈25 — **APROVADO** (predicao 1.4-1.7 ERROU: a base nao
   fechou; a estimativa estatica superestimou a fusao na base). Amplitude 0.168 = 0.54x o P2 —
   **REPROVADO** (limite 0.6x; predicao ~0.20 tambem otimista). Agar limpo nas baias 65% -> 42%.
5. Numerica: **APROVADO**.
6. Visual (§11): **APROVADO** — 18-20 bracos grossos, retos e continuos, baias abertas ate o
   nucleo; classe (b). Figuras: `plots/fig_P2R5_t50.png`, `plots/cmp_colonia_P2R5_t50.png`.

Numero de bracos caiu pouco (0.65 R99: 26.4 -> 21.8; 0.8: 19.6 -> 18.0). A queda de amplitude e
de agar limpo e geometrica: com a mesma contagem de bracos, bracos mais grossos ocupam mais do
perfil angular R(theta) e da area das baias.
