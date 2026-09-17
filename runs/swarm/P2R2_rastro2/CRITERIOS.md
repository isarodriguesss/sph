# P2R2 — rastro largo com W = 2 dx

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica sobre o P2R: `RASTRO_W` 1.0 -> 2.0
dx (`RASTRO_R_MIN` = 0.7). t=50, `SEED` fixo, 10 threads. Medir com o renderizador corrigido
(licao #92).

## Motivacao (licao #93)
Com W = 1 dx a conversao pegou o agar a frente do lider e as falhas do rastro: bracos continuos
(soltos 9.3 -> 0.3), mas pontas com a mesma largura (2.5 dx). O canal que a ponta abre ja tem
~2 dx; W = 2 dx alcanca o agar dos LADOS do canal.

## Predicao
- Largura p50 em r/R99 = 0.8: 2.5 -> 3.5-4.5 dx; a meia altura 2.2 -> 3-4 dx.
- Conversoes por chamada ~3x o P2R (area do disco 4x, menos o canal ja convertido).
- Custos maiores que no P2R: vivas abaixo de 903 (P2: 1206), amplitude 0.8-0.9x o P2, dedos
  28-33, frac(sigma_a<0.85) acima de 12%. Continuidade igual ao P2R (soltos <= 1).

## Criterios (contra P2 e P2R, janela t in [35, 50])
1. **Objetivo:** largura p50 >= 3.1 dx em r/R99 = 0.8.
2. Forma: amplitude em R99=2.4 >= 0.25 (0.8x o P2); bracos em r/R99=0.8 >= 12; soltos <= 2.6.
3. Disco nao antecipado: R_min(t=50) <= 1.8 (P2 1.23, P2R 1.53, com o renderizador corrigido).
4. Numerica: iter ate t=50 <= 1.5x o P2; massa sem runaway.
5. Visual (§11): classe (b) Fingering, bracos visivelmente mais grossos, baias abertas.
