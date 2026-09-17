# P2R9 — P2R8 + piso de `rho_b` no filler (0.4)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R8:
`FILLER_RHO_B_FLOOR` 0 -> 0.4 — todo filler novo (agar convertido pelo rastro, wake e insert)
nasce com `rho_b` = max(herdado, 0.4). t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (P2R8, t=50; regioes roxas no viewer marcadas pela usuaria)
Zonas so de filler na base dos bracos com `rho_b` baixo: nas regioes (-1.7, 0) e (1.8, 0.25) o
filler tem `rho_b` p50 0.136-0.143 -> `fade` da EOS 0.02-0.03 (2-3% da coesao/pressao). Origem:
o filler herda o `rho_b` de quem o cria, e a base foi construida cedo (t≈10-25) por lideres com
`rho_b` ≈ 0.13; os bracos, depois, por lideres a ~0.5. Na base (0.2-0.6 R99) 6% do filler tem
`fade` < 0.5 (p10 0.335); nos bracos (0.6-1.0 R99) 0%. Precedente: licao #67-N (piso 0.4 no
filler do wake do C4): AR 6.59 com 24 bracos, coesao da juncao 2.7x, `cs` intacto.

## Predicao
- Filler da base (0.2-0.6 R99): `rho_b` p10 0.335 -> >= 0.40; frac `fade` < 0.5: 0.06 -> 0.
- Viewer: manchas roxas da base viram turquesa (visual).
- Morfologia igual ao P2R8 dentro do ruido (largura 5.5, 18 bracos em 0.8 R99, 0 pedacos).
- `cs` inalterado (filler fora da quimica). Colonizacao: doador mais denso na base ->
  vivas 590 -> 590-800; `Rmin` 1.28 -> 1.28-1.5.

## Criterios (janela t in [35, 50], contra P2R8)
1. **Objetivo:** filler da base com `fade` < 0.5 = 0% e `rho_b` p10 >= 0.40.
2. Nao virar disco (licao #90/#94-E): vivas <= 1000; `Rmin`(t=50) <= 1.6; baias abertas na base
   >= 0.85x o P2R8 (28.7).
3. Manter P2R8: pedacos soltos <= 0.5; largura p50 em 0.8 R99 >= 4.5 dx; tortuosidade <= 1.08.
4. Numerica: iter/t <= 1.5x o P2, massa sem runaway, `a_press_med` <= 2x o P2R8.
5. Visual (§11 + viewer): base sem manchas de `rho_b` baixo; classe (b).

## Resultado (2026-09-15, t=50, 1550 s) — renderizador final, janela t in [35, 50]

| | P2R8 | **P2R9** |
|---|---:|---:|
| filler da base (0.2-0.6 R99): `rho_b` p10 / p50 / frac `fade` < 0.5 | 0.335 / 0.410 / 6% | **0.400 / 0.412 / 0%** |
| vivas (t=50) / Rmin (t=50) | 590 / 1.28 | 622 / 1.17 |
| baias abertas / frestas na base / bolsas fechadas | 33.8 / 8.6 / 1.4 | 34.2 / 9.8 / **0.3** |
| pedacos soltos / largura 0.8 R99 / tortuosidade | 0.0 / 5.5 / 1.00 | 0.0 / 5.4 / 1.01 |
| bracos 0.65 / 0.8 R99 / amplitude / agar limpo | 25.8 / 18.2 / 0.156 / 46% | 26.5 / 20.0 / 0.154 / 47% |
| iter/t / massa / `a_press_med` / `contrast_cs` t=50 | 59.7 / 203.2 / 1.4e-3 / 19.5 | 62.1 / 204.2 / 2.2e-3 / 18.4 |

1. Objetivo: **APROVADO** — nenhum filler da base sem coesao; p10 0.400.
2. Nao virar disco: **APROVADO** — vivas 622, Rmin 1.17, baias abertas 34.2 (limite 28.7).
3-4. **APROVADOS** (0 pedacos, largura 5.4, tortuosidade 1.01; `a_press_med` 1.6x o P2R8 < 2x).
5. Visual (estilo viewer): **APROVADO** — a mancha de `rho_b` baixo da base sumiu; mesma forma.

**Achado fora dos criterios — o nucleo pinado DESCOLA do corpo na escala de contato.** C5 com
ligacao de 1.05 dx: P2R8 C5a 0.29-0.48 / C5b 9-11%; P2R9 **0.07-0.14 / 0.7-1.2%**. O componente do
centro vira so o nucleo (76 particulas, r <= 0.35, 43 das 44 pinadas), separado do resto por uma
fresta de ~1.05-1.2 dx (7 vizinhos a < 1.2 dx contra 59 no P2R8; densidade de colonia em r 0.20-0.35
cai para 0.8 part/dx²). A 1.4 dx os dois runs se equivalem (C5b 56-97% contra 63-91%) e a 2 dx
sao identicos (99.7%). Espacamento de vizinho mais proximo igual nos dois (nao e efeito de
compressao). Leitura: o filler agora tem pressao da EOS e e empurrado para longe do nucleo, que
esta pinado — o descolamento do nucleo da licao #73 (E11 tardio), antecipado. Nao aparece no painel
(a) (nao ha agar na fresta).

Figuras: `plots/fig_P2R9_t50.png` (figura de tese).
