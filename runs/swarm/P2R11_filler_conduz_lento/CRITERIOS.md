# P2R11 — P2R10 com conducao de `c_s` pelo filler mais lenta (D = 0.02)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R10: `FILLER_CS_D`
0 -> 0.02 — pares que envolvem filler difundem com D = 0.02 (1/4 do agar, 13x o biofilme) em vez
de `D_ext` = 0.08. Decaimento do filler igual (0.5 lambda); filler segue sem produzir e sem
Marangoni. Resto = P2R9 + `FILLER_CS_CONDUZ`=1. t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (P2R10)
A conducao pelo filler como agar (L_D = 1.03) levou o `c_s` do nucleo ate o agar da base (0.29,
67% acima do gate da colonizacao 0.3 em t≈36) e reabriu o laco do disco (licao #90/#94-H).

## Predicao (queda exponencial do nucleo, calibrada no P2R10: `c_s` = 0.47·exp(−d/L), L = sqrt(D/0.075))
- L = 0.52. `c_s` do agar da base (0.2-0.4 R99): 0.29 -> ~0.18 (0.12-0.22); frac acima de 0.3 em
  t≈36: 67% -> <= 20%.
- `c_s` do filler: 0.2-0.4 R99 ~0.20-0.30; 0.4-0.6 R99 ~0.08-0.14 (P2R9 agar 0.05; P2R10 0.22).
- Sem disco: Rmin ~1.2-1.35, bracos em 0.8 R99 ~18-20, R99 ~4.5-4.7.
- Painel (d) com `--cs-filler`: continuo, intermediario entre P2R9 e P2R10.

## Criterios (t=50 e janela t in [35, 50], contra P2R9 e P2R10)
1. **Objetivo — nao reabrir o laco:** agar+limbo da base acima do gate em t≈36 <= 25%; Rmin(t=50)
   <= 1.4; bracos em 0.8 R99 >= 18; R99(t=50) >= 4.5.
2. **Objetivo — painel:** `c_s` do filler decrescente do nucleo para fora e >= 0.08 em 0.4-0.6 R99
   (acima dos 0.05 do agar no P2R9); painel (d) sem descontinuidade entre nucleo e bracos.
3. Forma: pedacos soltos <= 0.5; largura p50 0.8 R99 >= 4.5 dx; tortuosidade <= 1.08.
4. Quimica: frac(`c_s` > 0.45) nas vivas >= 0.8; `contrast_cs` t=50 >= 12.9.
5. Numerica: iter/t <= 1.5x o P2, massa sem runaway.
6. Visual (§11): classe (b); painel (a) sem o disco do P2R10.

## Resultado (2026-09-15, t=50, 1656 s) — renderizador final, janela t in [35, 50]

| | P2R9 | P2R10 (D=0.08) | **P2R11 (D=0.02)** |
|---|---:|---:|---:|
| `c_s` do filler 0-0.2 / 0.2-0.4 / 0.4-0.6 / 0.6-0.8 / 0.8-1.0 R99 | congelado (0.29 / 0.18 / 0.14 / 0.12 / 0.10) | 0.47 / 0.41 / 0.22 / 0.12 / 0.09 | **0.46 / 0.29 / 0.10 / 0.07 / 0.07** |
| `c_s` do agar em volta 0.2-0.4 / 0.4-0.6 / 0.6-0.8 / 0.8-1.0 R99 | 0.045 / 0.049 / 0.059 / 0.049 | 0.291 / 0.185 / 0.112 / 0.069 | **0.152 / 0.091 / 0.070 / 0.055** |
| agar+limbo da base com `c_s` > 0.3, t≈36 | 5% | 67% | **7%** |
| mesma fracao t≈18 / 24 / 30 | 100% / 61% / 12% | 100% / 100% / 100% | 100% / 100% / 18% |
| R99 / Rmax / Rmin (t=50, `rank_runs`) | 4.74 / 4.99 / 1.17 | 4.31 / 4.54 / 1.64 | **4.68 / 4.96 / 1.36** |
| Rmin por outras reguas: media da janela (`largura_final`) / `cmp_colonia` | 1.05 / 1.12 | 1.58 / 1.65 | 1.44 / 1.55 |
| bracos 0.65 / 0.8 R99 / largura p50 0.8 R99 | 26.5 / 20.0 / 5.4 | 24.8 / 16.0 / 5.7 | 25.5 / **18.0** / 5.7 |
| frestas / baias abertas / bolsas de agar fechadas (0.4+0.5 R99) | 9.8 / 34.2 / 0.3 | 7.3 / 18.7 / 1.9 | 13.0 / 30.0 / **3.2** |
| pedacos / tortuosidade / amplitude / vivas | 0 / 1.01 / 0.154 / 622 | 0 / 1.03 / 0.170 / 668 | 0 / 1.01 / 0.159 / 554 |
| frac(`c_s` > 0.45) nas vivas / `contrast_cs` t=50 | 0.94 / 18.4 | 0.66 / 15.1 | **0.77** / 18.2 |
| `a_mar_bio_p95` (mediana t>30) / dR/dt / iter/t / massa | 2.05 / 0.096 / 62 / 204.2 | 1.72 / 0.089 / 81 / 205.6 | 2.23 / 0.093 / 65 / 204.5 |

**Predicao:** acertou dentro da faixa em todos os itens — agar da base 0.15 (previsto 0.12-0.22),
filler 0.29 e 0.10 (previsto 0.20-0.30 e 0.08-0.14), bracos 18 (18-20), R99 4.68 (4.5-4.7),
Rmin 1.36 (1.2-1.35, no limite de cima).

1. Objetivo (laco): **APROVADO no limite** — gate da base 7% em t≈36 (limite 25%), R99 4.68, bracos
   18.0 (limite 18), Rmin 1.36 pelo `rank_runs` (limite 1.4), mas 1.44-1.55 pelas outras duas
   reguas: o nucleo e intermediario entre o P2R9 e o P2R10, visivelmente maior que o do P2R9.
2. Objetivo (painel): **APROVADO** — `c_s` do filler decrescente, 0.10 em 0.4-0.6 R99 (limite
   0.08); painel (d) continuo, com o halo acompanhando os bracos (`plots/fig_P2R11_t50.png`).
3. Forma: aprovado (0 pedacos, 5.7 dx, tortuosidade 1.01).
4. Quimica: **REPROVADO no limite** em frac(`c_s` > 0.45) nas vivas, 0.77 (limite 0.8);
   `contrast_cs` passa.
5. Numerica: aprovado (iter/t 1.15x o P2, massa estavel).
6. Visual: classe (b), sem o disco do P2R10; base mais cheia que a do P2R9
   (`plots/cmp_colonia_P2R11_t50.png`).

**Onde a base fechou:** nos tres runs o gate de colonizacao esta aberto na base ate t≈18 (colonia
pequena, `c_s` alto em toda parte). No P2R9 ele fecha entre t=19 e 29; no P2R11 fica aberto ate
t≈24 e fecha entre 24 e 30; no P2R10 so fecha depois de t=36. A janela extra de ~5 s explica o Rmin
intermediario e as 3.2 bolsas de agar fechadas na base (P2R9: 0.3). Depois de t=30 a conducao lenta
nao reabre o laco.
