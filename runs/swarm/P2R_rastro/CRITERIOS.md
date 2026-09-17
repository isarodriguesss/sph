# P2R — P2 com rastro largo (conversao lateral do agar em volta das vivas da frente)

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica: `RASTRO_W` 0 -> 1.0 dx
(`RASTRO_R_MIN` = 0.7) sobre o P2. t=50, `SEED` fixo, 10 threads.

## Mecanismo
A cada chamada do wake (100 iteracoes), o agar (`rho_b` = 0, nao filler) a menos de W = 1 dx
de uma viva (`rho_b` >= 0.1, nao filler) com r >= 0.7*R99 vira filler do wake, herdando o
`rho_b` da viva. Nao cria particula (massa inalterada), nao cria limbo, filler e transparente
ao `cs` (licao #39). O gate radial evita as bordas das baias perto do nucleo, onde a conversao
alimentaria o laco do disco (licao #90).

## Diagnostico que motiva (licao #91 e perfil transversal)
Em r/R99 = 0.6-0.88 o braco e 80% filler do wake e 1-2% viva, com 1-2 dx de largura: e o
rastro de 1-2 vivas. O agar dos lados tem `cs` 0.10-0.17 < `COL_CS_MIN` = 0.3, entao 0% e
recrutavel. Largura (janela t in [35,50]): 2.5 dx nas pontas, 2.2-2.5 a meia altura.

## Predicao
- Largura p50 em r/R99 = 0.8: 2.4-2.5 -> 3.5-4.0 dx; em 0.65 idem proporcional.
- ~3 conversoes por lider por chamada, ~100 por chamada em regime: o filler nos bracos quase
  dobra. Massa identica ao P2 (conversao nao cria massa).
- Risco (lei L, licao #67-L): amplitude e dedos caem. Filler e barreira para a difusao de
  `cs`: bracos mais largos sao barreiras mais largas.

## Criterios (contra P2, janela t in [35, 50])
1. **Objetivo:** largura p50 >= 3.1 dx em r/R99 = 0.8 (e >= P2 em 0.65).
2. Forma: amplitude em R99=2.4 >= 0.8x o P2 (>= 0.25); bracos em r/R99=0.8 >= 12; bracos soltos <= 2.6.
3. Disco nao antecipado: `R_min(t=50)` <= 2.1.
4. Numerica: iter ate t=50 <= 1.5x o P2; massa sem runaway.
5. Leitura visual (§11): classe (b) Fingering, bracos visivelmente mais grossos.

## Resultado (2026-09-14, t=50, 1612 s) — medido com o renderizador CORRIGIDO (licao #92)

| janela t in [35,50] | P2 | P2 (2a realiz.) | **P2R** |
|---|---:|---:|---:|
| largura p50 r/R99=0.8 (dx) | 2.5 | 2.4 | 2.5 |
| largura p50 / media r/R99=0.65 (dx) | 1.9 / 2.0 | 1.7 / 1.8 | 2.2 / 2.5 |
| bracos em r/R99 = 0.65 / 0.8 | 16.5 / 12.5 | 14.4 / 13.0 | **25.8 / 18.3** |
| bracos soltos / fora do corpo | 9.3 / 14.8% | 8.8 / 15.8% | **0.3 / 0.4%** |
| C5a / C5b (particulas, 1.05 dx) | 0.45 / 17.0 | 0.43 / 12.8 | 0.30 / 9.6 |
| amplitude (R99=2.4) / dedos | 0.312 / 36 | 0.317 / 37 | 0.286 / 33 |
| R99 / dR/dt / vivas (t=50) | 3.99 / 0.078 / 1206 | 3.77 / 0.082 / 1149 | 4.48 / 0.090 / 903 |
| R_min em t=50 | 1.23 | 1.36 | 1.53 |
| frac(sigma_a<0.85) / void_10 (log, t~48-50) | 0.060 / 0.012 | — | 0.121 / 0.005 |

1. Objetivo (largura >= 3.1 dx em r/R99=0.8): **REPROVADO** — 2.5, igual ao P2. A meia altura +15-25%.
2. Forma: amplitude 0.92x, 18 bracos, soltos 0.3 — aprovado.
3. R_min 1.53 <= 2.1 — aprovado (baias ~0.3 mais rasas que o P2).
4. Numerica: iter/t ~1.13x o P2, massa 205 — aprovado.
5. Visual: dendritico, bracos CONTINUOS e ligados ao corpo; nao mais grossos.

Leitura: W = 1 dx converte sobretudo o agar A FRENTE do lider e as falhas do rastro — o canal
que a ponta abre ja tem ~2 dx. O resultado e continuidade (o colar de grumos vira braco
continuo), nao largura. Custos: -25% de vivas (agar convertido em filler deixa de ser
recrutavel), frac(sigma_a<0.85) 2x.
