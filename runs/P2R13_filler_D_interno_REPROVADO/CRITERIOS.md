# P2R13 — a matriz conduz `c_s` melhor POR DENTRO, sem abrir a troca com o agar

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R10** (o melhor do ranking contra a
literatura, licao #95). Alavanca unica: a difusividade de `c_s` entre pares **filler-filler** passa
de `D_ext` = 0.08 para **0.25**; o par filler-agar e o par filler-viva continuam em `D_ext`
(`FILLER_CS_D` = 0.25 com `FILLER_CS_D_INTERNO` = 1, parametro novo; com `INTERNO` = 0 o codigo
reproduz o P2R11/P2R12, que aplicavam o valor a qualquer par com filler). t=50, `SEED` fixo,
10 threads.

## Por que esta alavanca, e nao mais conducao

A serie P2R10/11/12 e UM eixo — a difusividade do filler (0.08 / 0.02 / 0.01) — e ele ja esta no
limite: o `c_s` do corpo sobe (0.19 -> 0.22 -> 0.47 do teto) mas o nucleo cresce junto
(0.29 -> 0.31 -> 0.38 Rmax) porque a conducao tambem vaza para o agar da base, que cruza o gate da
colonizacao e realimenta o laco do disco (licao #90/#94-H). Como a licao #95 mediu que o nucleo do
P2R10 (0.38) JA e o alvo da literatura (0.36), subir o mesmo eixo passaria do alvo.

Separando os dois papeis: a matriz de EPS conduz bem internamente, e a troca com o agar segue
limitada pela difusividade do agar (o par misto fica no valor mais lento, no espirito da media
harmonica). Alvo da licao #95: Trinschek tem Gamma ≈ Gamma_max em toda a pegada da colonia
(corpo 0.99, baias 0.95 do teto) com o gradiente so na borda externa.

## Predicao
- Comprimento de conducao ao longo do braco: sqrt(D/0.5 lambda) = 1.83 contra 1.03 do P2R10 (1.8x).
- `c_s` do filler por faixa (0.2-0.4 / 0.4-0.6 / 0.6-0.8 R99): 0.41 / 0.22 / 0.12 (P2R10) ->
  **0.42 / 0.30-0.38 / 0.18-0.25**. Corpo (mediana sobre a colonia, `tools/lit_rank.py`):
  0.47 -> **0.60-0.75** do teto.
- Baias: 0.21 -> 0.25-0.35 do teto (Trinschek 0.95).
- Agar da base (0.2-0.4 R99): sobe pouco, 0.29 -> 0.31-0.36 (o par com o agar nao mudou, mas o
  interior mais cheio vaza mais); gate da base em t≈36: 67% -> 70-85%.
- Nucleo: 0.38 -> 0.40-0.45 Rmax. Forma no resto igual ao P2R10 dentro do ruido.

## Criterios (t=50 e janela t in [35, 50], contra P2R10 e a faixa da literatura da licao #95)
1. **Objetivo:** `c_s` do corpo >= 0.60 do teto (P2R10 0.47; alvo 0.99) e nas baias >= 0.25
   (P2R10 0.21; alvo 0.95), com o perfil do filler decrescente do nucleo para fora.
2. **Guarda do nucleo:** nucleo <= 0.45 Rmax (P2R10 0.38; alvo 0.36) — acima disso a colonia sai da
   faixa da literatura e a alavanca esta esgotada.
3. Forma dentro da faixa da licao #95: area/disco 0.33-0.42, dedos >= 15, AR <= 20, R99 >= 4.1.
4. Quimica: `contrast_cs` >= 12.9; frac(`c_s` > 0.45) nas vivas >= 0.6 (P2R10 0.66).
5. Numerica: iter/t <= 1.6x o P2, massa sem runaway.
6. Visual (§11): classe (b); painel (d) com o corpo mais uniforme que o do P2R10.
