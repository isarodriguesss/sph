# P2S — P2 com deposicao do wake ao longo do rastro (`WAKE_SEG`)

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica: `WAKE_SEG` False -> True
(`WAKE_SEG_MAX`=6) sobre o P2. t=50 (janela morfologica do P2, licao #90), `SEED` fixo,
10 threads — mesma condicao do `runs/P2_fillerdonor`.

## Objetivo (pedido da usuaria)
Engrossar um pouco os bracos antes da rugosidade.

## Diagnostico que motiva
Largura dos bracos no campo da figura de tese (arcos contiguos de colonia num circulo de
raio r), P2 em t~50: **2.5-2.8 dx** em r/R99 = 0.65 e 0.8, com 15-23 bracos e ocupacao
azimutal 0.11-0.25. Ruido entre realizacoes (P2_fillerdonor x P2_t100): ~0.1 dx.
Os bracos sao mais finos que o suporte do kernel (2h = 3.6 dx).

O wake deposita o grumo so na posicao ANTIGA da mae (`x_dep`), e a ponta anda 2-10 dx
entre chamadas. `WAKE_SEG` deposita tambem ao longo do segmento `x_dep -> x`. No C4 (licao
#65) isso deu bracos nitidamente mais grossos com AR 10.48 -> 7.83 (largura ~+33%), massa
+0.2%, `a_mar` na frente +18%, `frac(sigma_a<0.85)` 7.35% -> 5.04%, dedos 19 -> 18.

O teto por chamada `WAKE_MAX`=150 ja trava em 25 de 72 chamadas no P2 (media 114): o
efeito sera mais REDISTRIBUIR onde o wake cai (ao longo do rastro, nao em anel) do que
somar material.

## Predicao
- Largura p50 em r/R99 = 0.65 e 0.8: 2.5-2.8 -> **3.2-3.6 dx** (+25-35%).
- Numero de bracos em r/R99 = 0.8: 15 -> 13-15. `R99` em t=50 dentro de -5%.
- `a_mar` na frente sobe; `frac(sigma_a<0.85)` cai (braco mais largo = suporte melhor).
- Wake adicionado ate t=50: +10-30% sobre o P2 (teto por chamada), longe do orcamento.
- Risco (licao #90): wake e doador no laco que fecha as fendas perto do nucleo; bracos
  mais largos estreitam as baias. `R_min` pode subir mais cedo.

## Criterios (contra P2, janela t in [35, 50])
1. **Objetivo:** largura p50 >= 3.1 dx em r/R99 = 0.65 e 0.8 em t~50.
2. Forma: bracos em r/R99=0.8 >= 12; amplitude em R99=2.4 >= 0.8x o P2; bracos soltos <= 2.6.
3. Disco nao antecipado: `R_min(t=50)` <= 2.1 (P2: 1.69-1.81).
4. Numerica: sem colapso de dt (iter/t ate t=50 <= 1.5x o P2), massa sem runaway,
   `a_press_med` <= 2x o P2.
5. Leitura visual (§11): segue na classe (b) Fingering, bracos mais grossos, baias abertas.

## Resultado (2026-09-14, t=50, 1465 s)

| janela t in [35,50] | P2 | P2 (2a realiz.) | **P2S** |
|---|---:|---:|---:|
| largura r/R99=0.65: p50 / media / p75 (dx) | 2.5 / 4.3 / 4.4 | 2.2 / 3.3 / 3.5 | **2.9 / 5.6 / 6.9** |
| largura r/R99=0.8: p50 / media / p75 (dx) | 2.5 / 2.5 / 2.8 | 2.5 / 2.4 / 2.7 | 2.4 / 2.6 / 2.8 |
| bracos em r/R99=0.8 | 13.5 | 13.2 | 13.0 |
| amplitude (R99=2.4) / dedos | 0.312 / 36 | 0.317 / 37 | 0.293 / 40 |
| bracos soltos / C5a / C5b / C5b+limbo | 1.3 / 0.45 / 17.0 / 37.9 | 2.6 / 0.43 / 12.8 / 16.1 | 1.8 / 0.45 / **25.7 / 43.9** |
| R_min em t=50 | 1.81 | 1.53 | **1.34** |
| wake depositado ate t=50 | — | 2706 | 2942 (+9%) |

1. Objetivo (>= 3.1 dx em 0.65 E 0.8): **REPROVADO.** Engrossou so a meia altura (+20-25% na
   mediana, cauda +60%); as pontas nao mudaram.
2-4. Aprovados: forma preservada (amplitude 0.93x), disco ATRASADO (R_min 1.34), continuidade
   melhor, dt e massa como no P2.
5. Visual: mesma classe do P2; bracos nao visivelmente mais grossos.

Causa: o caminho da ponta nao e vazio — e AGAR, que a colonia atravessa sem deslocar (licao
#58). O `livre()` so aceita ponto sem vizinho a 0.7 dx, entao o segmento so enche vazio real:
+9% de wake. A largura do braco nao e fixada pela deposicao, e sim por ate onde o agar lateral
e RECRUTADO — e recrutamento lateral e o mecanismo do disco tardio (licao #90).
