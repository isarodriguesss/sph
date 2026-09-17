# P2R8 — P2R6 + lideres internos em movimento (pontas dos dendritos novos)

Pre-registrado em 2026-09-15, antes de lancar. Base escolhida pela usuaria: P2R6 (W = 3,
`RASTRO_SEG`, guarda de baia G = 2 dx na versao ORIGINAL, `RASTRO_BAIA_CHAMADA`=False). Alavanca
unica: `RASTRO_R_MOV` 0 -> 0.4 — alem dos lideres a r >= 0.7 R99, entram as vivas (`rho_b`>=0.1)
entre 0.4 e 0.7 R99 que andaram >= 0.5 dx desde a chamada anterior com componente radial >= 0.5
do passo. t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (P2R5/P2R6, regioes A, C, E marcadas pela usuaria)
Dendritos novos nascem por dentro (r/R99 ≈ 0.5-0.6), abaixo do limiar de 0.7 R99 da conversao: o
lider anda (|v| 0.02-0.05) e deixa so o rastro tracejado do wake, com agar entre os grumos. Base no
P2R6 (t 30-50): 12 de 90 passos de lideres internos atravessam agar (eixo do caminho > 50% agar).
Checagem da selecao nos quadros do P2R6 (it 2400 e 2800): 13 lideres internos novos por chamada,
em (-2.4, 0.8), (2.6, -0.4), (-0.4, -2.5) etc. — as regioes marcadas.

## Predicao
- Passos internos em agar: 12/90 -> <= 3; pedacos soltos >= 4 dx²: 1.3 -> <= 0.5.
- Os dendritos novos aparecem como bracos continuos: bracos em 0.65 R99 23.8 -> 24-28.
- Baias vizinhas a eles estreitam (guarda G=2 limita a ~2 fileiras): frestas < 2 dx na base
  8.0 -> 8-11; baias abertas 36.7 -> 32-37.
- Filler t≈41: 7299 -> +5-15%. Largura externa 5.5 inalterada; tortuosidade ~1.01.

## Criterios (janela t in [35, 50] / t 30-50 para os passos, contra P2R6)
1. **Objetivo:** passos de lideres internos com eixo > 50% agar <= 3 E pedacos soltos <= 0.5.
2. Nao piorar a base: frestas < 2 dx <= 1.4x o P2R6 (11.2) e baias abertas >= 0.85x (31.2).
3. Largura p50 em 0.8 R99 >= 4.5 dx; tortuosidade <= 1.08.
4. `Rmin`(t=50) <= 2.0; numerica: iter/t <= 1.5x o P2, massa sem runaway.
5. Visual (§11 + estilo viewer): pontas dos dendritos novos continuas; classe (b).

Nota: com `RASTRO_BAIA_CHAMADA`=False e `RASTRO_R_MOV`=0 o codigo atual segue a logica do P2R6
(a guarda foi reescrita sem mudar a regra; empates exatos na distancia g poderiam diferir).

## Resultado (2026-09-15, t=50, 1468 s) — renderizador final, janela t in [35, 50]

| | P2R6 | **P2R8** |
|---|---:|---:|
| passos de lideres internos com eixo > 50% agar (t 30-50) | 12 / 90 | **3 / 83** |
| pedacos soltos >= 4 dx² / > 20 dx² | 1.3 / 0.0 | **0.0 / 0.0** |
| frestas < 2 dx / baias abertas >= 2 dx (0.4+0.5 R99) | 8.0 / 36.7 | 8.6 / 33.8 |
| bolsas de agar fechadas | 1.8 | 1.4 |
| largura p50 0.65 / 0.8 R99 (dx) | 5.4 / 5.5 | 5.4 / 5.5 |
| bracos em 0.65 / 0.8 | 23.8 / 18.2 | 25.8 / 18.2 |
| tortuosidade / Rmin (t=50) | 1.01 / 1.28 | 1.00 / 1.28 |
| amplitude / agar limpo / vivas | 0.156 / 46% / 609 | 0.156 / 46% / 590 |
| filler t≈40 / iter/t / massa / `a_press_med` | 7299 / 60.4 / 204.2 / 1.5e-3 | 7224 / 59.7 / 203.2 / 1.4e-3 |

1. Objetivo: **APROVADO** — 3 passos internos em agar (limite 3) e 0.0 pedacos soltos (limite 0.5).
2. Base: **APROVADO** — frestas 8.6 (limite 11.2), baias abertas 33.8 (limite 31.2).
3-4. **APROVADOS** — largura 5.5, tortuosidade 1.00, Rmin 1.28, numerica igual ao P2R6.
5. Visual (§11): **APROVADO** — dendritos novos continuos (ex. (-2.8, 0.9), (0.3, 2.5),
   (2.5, -0.5)), sem traços tracejados, baias abertas. Figuras: `plots/cmp_colonia_P2R8_t50.png`,
   `plots/fig_P2R8_t50.png`.
Predicoes: bracos em 0.65 R99 25.8 (24-28, ok); frestas e baias dentro do previsto; filler −1%
(previ +5-15% — a conversao nova e pequena perto do total).
