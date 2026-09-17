# P2R7 — guarda de baia corrigida e G = 4 dx (sobre o P2R6)

Pre-registrado em 2026-09-15, antes de lancar. Contra o P2R6: (a) CORRECAO de implementacao —
a guarda passa a considerar tambem o agar convertido por outros lideres NA MESMA chamada (no
P2R6 a arvore era montada no inicio da chamada; 46% dos pares dos dois lados das frestas da base
foram convertidos no mesmo quadro); (b) alavanca: `RASTRO_BAIA` 2 -> 4 dx (com 2 dx a guarda
preserva ~2 fileiras de agar = fenda de 1.6-2.6 dx no painel (a), exatamente a "baia unida").
W = 3, `RASTRO_SEG`, t=50, `SEED` fixo, 10 threads.

Teste sintetico: dois lideres paralelos a 8 dx convertendo na mesma chamada deixam 3 fileiras de
agar entre si (sem guarda: 1); vizinho fixo com G=2 reproduz o resultado do P2R6.

## Predicao
- Baias da base (0.4 + 0.5 R99): frestas < 2 dx 8.0 -> <= 3; baias abertas >= 2 dx 36.7 -> >= 36;
  bolsas de agar fechadas 1.8 -> <= 0.5.
- Baia p50 na base 2.8 -> 4-5 dx.
- Largura p50 em 0.8 R99: 5.5 -> 5.0-5.6 (baias de 10-18 dx, guarda age pouco); em 0.65: 5.4 ->
  4.5-5.4.
- Amplitude 0.156 -> 0.17-0.21; agar limpo 46% -> 50-60%; filler t≈41 7299 -> 10-25% menos.

## Criterios (janela t in [35, 50], contra P2R5/P2R6)
1. **Objetivo:** frestas < 2 dx em 0.4+0.5 R99 <= 3.6 (0.5x o P2R5) E bolsas fechadas >= 4 dx²
   <= 0.5 E baias abertas >= 2 dx >= 30.3 (P2R5).
2. Largura p50 em 0.8 R99 >= 4.5 dx.
3. Tortuosidade p50 <= 1.08; pedacos soltos <= 2 por quadro.
4. `Rmin`(t=50) <= 2.0; numerica: iter/t <= 1.5x o P2, massa sem runaway.
5. Visual (§11 + estilo viewer): baias da base abertas, sem frestas; classe (b).

## Resultado (2026-09-15, t=50, 1506 s) — renderizador final, janela t in [35, 50]

| | P2R5 | P2R6 (G=2) | **P2R7 (G=4, corrigida)** |
|---|---:|---:|---:|
| frestas < 2 dx em 0.4+0.5 R99 | 7.2 | 8.0 | **2.5** |
| baias abertas >= 2 dx em 0.4+0.5 R99 | 30.3 | 36.7 | **6.7** |
| bolsas de agar fechadas | 1.0 | 1.8 | **0.0** |
| largura p50 0.65 / 0.8 R99 (dx) | 5.5 / 5.5 | 5.4 / 5.5 | 5.2 / 5.4 |
| bracos em 0.65 / 0.8 | 21.8 / 18.0 | 23.8 / 18.2 | 20.3 / **14.5** |
| Rmin (t=50) / amplitude / agar limpo | 1.17 / 0.168 / 42% | 1.28 / 0.156 / 46% | **1.78** / 0.202 / 53% |
| vivas (t=50) / C5a / C5b | 598 / 0.40 / 8.0 | 609 / 0.34 / 8.4 | **1703** / 0.63 / 22.6 |
| borda das baias da base: convertido / depositado / RECRUTADO | 89 / 8 / 2% | 88 / 10 / 2% | **49 / 13 / 38%** |
| filler t≈41 / iter/t / massa / `a_press_med` | 7564 / 60.3 / 204.3 / 1.4e-3 | 7299 / 60.4 / 204.2 / 1.5e-3 | 5882 / 57.9 / 206.1 / 2.3e-3 |

1. Objetivo: **REPROVADO** — frestas 2.5 e bolsas 0.0 passam, mas baias abertas 6.7 (alvo >= 30.3):
   as baias da base nao ficaram abertas, FECHARAM de vez.
2-4. Largura 5.4, tortuosidade 1.00, 1.7 pedacos, Rmin 1.78 <= 2.0, numerica — aprovados.
5. Visual: **REPROVADO** — nucleo virou disco de raio ~2.2 (o do P2 e 1.8); bracos saem da borda.

**Causa (medida):** a guarda cumpriu o papel — a borda das baias da base caiu de 89% para 49% de
agar convertido. O agar que ela preservou foi RECRUTADO pela colonizacao (38% da borda; vivas
598 -> 1703; 173 no limbo): e o laco filler-doador -> recruta da licao #90, que no P2R5 nao tinha
o que recrutar porque o rastro ja tinha convertido esse agar em filler (nao recrutavel). Na base,
o agar que o rastro nao converte vira colonia pela colonizacao — a guarda troca um preenchimento
pelo outro. Predicoes que acertaram: amplitude 0.202 (0.17-0.21), agar limpo 53% (50-60%), filler
−19% (−10 a −25%). A que errou: baias abertas (previ >= 36, deu 6.7).

Figura: `plots/cmp_colonia_P2R7_t50.png` (P2R5 / P2R6 / P2R7, painel (a), t=50).
