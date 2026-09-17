# P2R12 — P2R11 com conducao de `c_s` pelo filler ainda mais lenta (D = 0.01)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R11: `FILLER_CS_D`
0.02 -> 0.01 (1/8 do agar, 7x o biofilme). Resto identico ao P2R11 (`src/` igual byte a byte).
t=50, `SEED` fixo, 10 threads.

## Motivo
O P2R11 manteve o painel (d) e quase fechou o laco do disco, mas o nucleo ficou intermediario
(Rmin 1.36-1.55 conforme a regua, contra 1.05-1.17 do P2R9), com 3.2 bolsas de agar fechadas na
base, e frac(`c_s` > 0.45) nas vivas em 0.77 (limite 0.8). Pergunta: ha um D em que o nucleo volta
ao do P2R9 sem perder o painel?

## Predicao (queda do nucleo `c_s` = 0.47·exp(−d/L), L = sqrt(D/0.075) = 0.365; distancias efetivas
d medidas no P2R10 e P2R11)
- Filler 0.4-0.6 R99 (d ≈ 0.8): **~0.05-0.07** — abaixo do limite de 0.08 do criterio 2, e perto do
  fundo de ~0.05 que o agar ja tem no P2R9. **Previsao: o criterio 2 reprova**; o painel deve ficar
  com o meio dos bracos proximo do P2R9.
- Filler 0.2-0.4 R99 (d ≈ 0.2-0.25): ~0.22-0.27.
- Agar da base 0.2-0.4 R99 (d ≈ 0.55): 0.15 -> **~0.09-0.11**.
- Gate da base (agar+limbo com `c_s` > 0.3): fecha entre t=19 (P2R9) e t=24-30 (P2R11); em t≈36 ~3-6%.
- Rmin (`rank_runs`) ~1.2-1.3; bolsas de agar fechadas 1-2.5; bracos em 0.8 R99 18-20; R99 4.6-4.8.
- frac(`c_s` > 0.45) nas vivas (interpolado em sqrt(D): 0.94 / 0.77 / 0.66 em D = 0 / 0.02 / 0.08):
  ~0.80-0.85.

## Criterios (t=50 e janela t in [35, 50], contra P2R9 e P2R11)
1. **Objetivo — nucleo como o P2R9:** Rmin (`rank_runs`) <= 1.3, reportando tambem as reguas de
   `largura_final` e `cmp_colonia`; bolsas de agar fechadas na base <= 2; bracos em 0.8 R99 >= 18;
   R99 >= 4.5; gate da base em t≈36 <= 10%.
2. **Objetivo — painel:** `c_s` do filler decrescente do nucleo para fora e >= 0.08 em 0.4-0.6 R99;
   painel (d) sem descontinuidade entre nucleo e bracos.
3. Forma: pedacos soltos <= 0.5; largura p50 0.8 R99 >= 4.5 dx; tortuosidade <= 1.08.
4. Quimica: frac(`c_s` > 0.45) nas vivas >= 0.8; `contrast_cs` t=50 >= 12.9.
5. Numerica: iter/t <= 1.5x o P2, massa sem runaway.
6. Visual (§11): classe (b); nucleo e baias da base comparaveis ao P2R9.

## Resultado (2026-09-15, t=50, 2090 s) — renderizador final, janela t in [35, 50]

| | P2R9 (D=0) | P2R11 (D=0.02) | **P2R12 (D=0.01)** |
|---|---:|---:|---:|
| `c_s` do filler 0-0.2 / 0.2-0.4 / **0.4-0.6** / 0.6-0.8 / 0.8-1.0 R99 | congelado (0.29 / 0.18 / 0.14 / 0.12 / 0.10) | 0.46 / 0.29 / **0.098** / 0.07 / 0.07 | 0.44 / 0.21 / **0.080** / 0.06 / 0.07 |
| `c_s` do agar em volta 0.2-0.4 / 0.4-0.6 R99 | 0.045 / 0.049 | 0.152 / 0.091 | 0.114 / 0.073 |
| gate da base (agar+limbo com `c_s` > 0.3) em t≈18 / 22 / 26 / 30 / 36 | 1.00 / 1.00 / 0.61 / 0.12 / 0.05 | 1.00 / 1.00 / 0.53 / 0.18 / 0.07 | 1.00 / 1.00 / 0.76 / 0.25 / 0.07 |
| **Rmin**: `rank_runs` / media da janela / `cmp_colonia` | 1.17 / 1.05 / 1.12 | 1.36 / 1.44 / 1.55 | **1.40 / 1.34 / 1.40** |
| bolsas de agar fechadas na base | 0.3 | 3.2 | **2.6** |
| frestas / baias abertas (0.4+0.5 R99) | 9.8 / 34.2 | 13.0 / 30.0 | 12.3 / 31.1 |
| R99 / Rmax / bracos 0.8 R99 / largura p50 | 4.74 / 4.99 / 20.0 / 5.4 | 4.68 / 4.96 / 18.0 / 5.7 | 4.81 / 4.99 / 18.7 / 5.6 |
| pedacos / tortuosidade / amplitude / vivas | 0 / 1.01 / 0.154 / 622 | 0 / 1.01 / 0.159 / 554 | 0 / — / 0.165 / 553 |
| frac(`c_s` > 0.45) nas vivas / `contrast_cs` | 0.94 / 18.4 | 0.77 / 18.2 | **0.80** / 19.2 |
| `a_mar_bio_p95` (mediana t>30) / dR/dt / iter/t / massa | 2.05 / 0.096 / 62 / 204.2 | 2.23 / 0.093 / 65 / 204.5 | 2.83 / 0.098 / 80 / 206.0 |

**Predicao vs medido:** acertou o agar da base (0.114, previsto 0.09-0.11), os bracos (18.7) e o
R99 (4.81); acertou o filler em 0.2-0.4 R99 (0.213, previsto 0.22-0.27, 3% abaixo). **Errou para
baixo** o filler em 0.4-0.6 R99 — previ 0.05-0.07 e deu **0.080**, exatamente no limite, entao o
criterio 2 que eu previa reprovar PASSOU por uma casa. **Errou o Rmin**: previ 1.2-1.3 e deu
1.34-1.40 — o nucleo NAO voltou ao do P2R9.

1. Objetivo (nucleo): **REPROVADO** — Rmin 1.40 (`rank_runs`, limite 1.3), 1.34 pela media da
   janela; bolsas 2.6 (limite 2). Passa em bracos (18.7), R99 (4.81) e gate em t≈36 (7%).
2. Objetivo (painel): **APROVADO no limite** — 0.080 em 0.4-0.6 R99; painel (d) contínuo
   (`plots/fig_P2R12_t50.png`), com o halo um pouco mais fraco no meio dos bracos que o do P2R11.
3. Forma: aprovado (0 pedacos, 5.6 dx).
4. Quimica: **APROVADO no limite** — 0.80 exatamente no limite; `contrast_cs` 19.2, o melhor dos tres.
5. Numerica: aprovado (iter/t 1.42x o P2, massa 206.0).
6. Visual: classe (b); o nucleo e visivelmente maior que o do P2R9 e parecido com o do P2R11
   (`plots/cmp_colonia_P2R12_t50.png`).

**O achado: o tamanho do nucleo NAO e proporcional a D.** Cortar D pela metade (0.02 -> 0.01)
mexeu o Rmin em ~0.1 e custou 18% do `c_s` do filler no meio dos bracos. Motivo medido: o gate da
base ja esta 100% aberto ate t≈22 nos TRES runs — inclusive no P2R9, que nao conduz nada. Enquanto
a colonia e pequena, o `c_s` produzido pelas vivas chega a base sozinho; a conducao pelo filler so
atrasa o fechamento de ~21-29 (P2R9) para ~22-29 (P2R12) e ~24-31 (P2R11). **A conducao nao e o que
abre a janela — ela apenas a mantem aberta um pouco mais, e por isso D nao e a alavanca para
devolver o nucleo do P2R9.** Reduzir D abaixo de 0.01 leva o painel ao nivel do P2R9 (o filler
encosta no fundo de 0.05 do agar) sem recuperar o nucleo: **a rota do D esta esgotada entre
0.01 e 0.02.**
