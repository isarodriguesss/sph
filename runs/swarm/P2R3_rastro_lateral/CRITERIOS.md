# P2R3 — rastro largo so dos lados (W = 2 dx, sem converter a frente)

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica sobre o P2R2: `RASTRO_COS_FRENTE`
1.0 -> 0.5 — o agar num cone de +-60 graus a frente da velocidade do lider NAO e convertido.
W = 2 dx, `RASTRO_R_MIN` = 0.7, t=50, `SEED` fixo, 10 threads. Medido com o renderizador final
(regiao sem agar, licao #92).

## Diagnostico que motiva (P2R2, t de 30 a 50)
Pontas serpenteiam: tortuosidade 1.22 (p75 1.44) contra 1.01-1.03 do P2/P2R; mudanca de rumo
25 graus por quadro contra 5-10; velocidade-raio cos 0.82 contra 0.98-0.99; vizinhas filler
dos lideres 18% contra 1-5%. 71% do agar a < 2 dx do lider esta A FRENTE dele: a ponta avanca
sobre filler que ela mesma criou, que nao entra no gradiente de `cs`.

## Predicao
- Tortuosidade p50 <= 1.06, mudanca de rumo <= 12 graus, cos velocidade-raio >= 0.95.
- Conversoes por chamada ~30-40% das do P2R2 (so o agar lateral e de tras).
- Bracos mais largos que o P2R e menos que o P2R2; relevo mais proximo do P2R.

## Criterios
1. **Objetivo (retos):** tortuosidade p50 <= 1.08 e cos velocidade-raio >= 0.93.
2. Nao perder a largura: bracos visivelmente mais grossos que no P2R (visual, §11).
3. Forma: sem fusao de bracos em alcas; relevo >= 0.8x o P2.
4. Numerica: iter/t <= 1.5x o P2, massa sem runaway.

## Resultado (2026-09-14, t=50, 1738 s) — renderizador final, janela t in [35, 50]

| | P2 | P2R (W=1) | P2R2 (W=2) | **P2R3** |
|---|---:|---:|---:|---:|
| tortuosidade p50 (t 32-50) | 1.01-1.03 | 1.01-1.03 | 1.22 | **1.03** |
| mudanca de rumo / cos velocidade-raio | 5-10° / 0.98 | 5-10° / 0.99 | 25° / 0.82 | **8° / 0.97** |
| largura p50 r/R99=0.65 / 0.8 (dx) | 1.6 / 2.1 | 2.2 / 2.7 | 4.3 / 3.9 | **3.2 / 3.3** |
| bracos em 0.65 / 0.8 | 10.0 / 11.3 | 24.8 / 18.5 | 18.0 / 18.8 | 23.4 / 18.0 |
| pedacos soltos >= 4 dx² / > 20 dx² | 35.0 / 8.7 | 25.7 / 15.0 | 0.5 / 0.0 | **14.4 / 10.6** |
| amplitude (R99=2.4) / dedos / agar limpo | 0.312 / 36 / 78% | 0.286 / 33 / 78% | 0.163 / 26 / 55% | **0.255** / 31 / 67% |
| R99 / Rmin (t=50) / dR/dt t[25,50] | 3.99 / 1.77 / 0.078 | 4.48 / 1.49 / 0.090 | 3.98 / 1.41 / 0.073 | 4.71 / 1.26 / 0.091 |
| vivas (t=50) | 1206 | 903 | 581 | 616 |
| iter/t / massa t=50 / `a_press_med` | 56.3 / 208.0 / — | 62.9 / 205.1 / 1.8e-3 | 69.0 / 202.3 / 2.1e-3 | 63.8 / 203.8 / 1.7e-3 |

1. Objetivo (retos): **APROVADO** — tortuosidade 1.03, cos 0.97, rumo 8°, nivel P2/P2R.
2. Largura: **APROVADO** — 3.2/3.3 dx contra 2.2/2.7 do P2R (1.5x/1.2x) e 1.6/2.1 do P2.
3. Forma: sem alcas; amplitude 0.82x o P2 — **aprovado no limite**; 31 dedos contra 36.
4. Numerica: **APROVADO** — iter/t 1.13x, massa abaixo do P2.

**Defeito novo, fora dos criterios:** os bracos saem CORTADOS no painel (a) — 14.4 pedacos
soltos contra 0.5 do P2R2. Causa medida em t≈43: 30 particulas de agar no corredor dos bracos,
0.73 dx a frente das posicoes de chamada (p50), 47% dentro da cunha frontal excluida: o agar
logo a frente do lider cai no cone na chamada k e ja esta a mais de 2 dx dele na k+1 (o lider
anda 2.9 dx por chamada, p90 3.5). Correcao em teste: P2R4 (conversao ao longo do segmento
percorrido). O P2R tem o mesmo defeito por outro motivo (vao entre discos de 1 dx com passo
2.6 dx): 25.7 pedacos — a "continuidade" da licao #93 era do renderizador com ponte de 1 dx.

Figura (classe das particulas nos cortes do painel (a), t≈43): `plots/zoom_quebras_P2R3.png`.
