# P2R22 — `WAKE_FREQ` 100 -> 50: fechar o vao entre o lider de ponta e o braco

Pre-registrado em 2026-09-17, ANTES de rodar. Base: **P2R21**, que passou a ser o baseline nesta data.
Alavanca unica: `WAKE_FREQ` 100 -> 50. t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva
Medido no P2R21 (frames 2200-2800, 83 pontas): o lider de ponta esta a **2.94 dx** do filler mais
proximo e anda **2.7 dx por chamada** do rastro — os dois numeros sao o mesmo. O vao NAO e
descolamento fisico: e a cadencia da conversao, que roda a cada `WAKE_FREQ` iteracoes. E a espicula
que aparece no viewer. (O outro defeito da ponta — o lider nao ter companhia, viva mais proxima a
**21.1 dx** — e estrutural, licao #98, e NAO e alvo desta rodada.)

## Teste de efeito (licao #102), feito antes — iteracao 600, t≈10
| | P2R21 (100) | teste (50) |
|---|---:|---:|
| vao lider -> filler p50 / p90 | 3.00 / 3.19 dx | **1.14 / 1.69 dx** |
| filler | 710 | 1007 (+42%) |
| wake+insert acumulado | 308 | 398 (+29%) |
| massa | 198.45 | 198.71 (+0.13%) |
| vivas / pontas / R99 | 175 / 9 / 1.10 | 181 / 13 / 1.18 |

## Predicao para t=50
- **Objetivo:** vao lider -> filler **1.0-1.6 dx** (P2R21 2.94); espicula visivelmente menor no painel (a).
- Filler +20-40%; massa dentro de 1% do P2R21 (teto do wake 0.12 x 197.5 = 23.7, hoje usado ~12).
- Largura EDT em 0.6 R sobe de 5.4 para 5.4-6.5 dx; ciclo 0.72 -> 0.72-0.85.
- **Risco 1 (lei (L), licao #67):** mais material por comprimento de frente engrossa e reduz relevo.
- **Risco 2 (licao #104):** o teste deu 13 pontas contra 9 em t≈10; mais pontas ativas fecham as baias
  a meia altura e matam a colonizacao — foi assim que o P2R19 reprovou.
- `dt` nao muda (a cadencia nao entra no CFL): iter/t ~60.

## Criterios (janela t in [35, 50]; frames em fase, iter 1600-2800)
1. **Objetivo:** vao lider -> filler (mediana, iter 2200-2800) **<= 1.8 dx**.
2. **Colonizacao viva:** limbo > 0 em TODOS os frames com t > 32 **e** vivas crescendo entre t=32 e t=50
   (P2R21: 501 -> 584).
3. **Relevo:** amplitude >= 0.8x a do P2R21; dedos entre 16 e 26 (P2R21: 20).
4. **Guarda de forma:** R99 >= 4.3 (P2R21 4.45); largura EDT em 0.6 R <= 6.5 dx; ciclo <= 1.0;
   nucleo <= 0.38 Rmax (P2R21 0.32); pedacos soltos <= 0.5.
5. **Motor:** `a_mar_bio_p95` (mediana t>30) >= 0.37 (P2R21).
6. **Agar a meia altura:** agar em 0.5-0.7 R99 na iteracao 1600 >= 74 (P2R21).
7. **Numerica:** `max_cs` <= 0.5; iter/t <= 78; massa sem runaway; teto do wake nao esgotado antes de t=50.
8. **Visual (§11):** classe (b) fingering preservada; espicula menor que a do P2R21 no mesmo instante.

**Nota de processo:** `tools/archive_run.sh` faz `rm -rf` no destino — ele APAGOU este arquivo quando
arquivei a rodada. Escrever o CRITERIOS.md so depois de arquivar, ou guardar copia fora de `runs/`.

## Resultado (2026-09-17, t=50) — OBJETIVO ATINGIDO, mas REPROVADO pelo motor

| | P2R21 (baseline) | **P2R22 (FREQ=50)** | criterio |
|---|---:|---:|---|
| **vao lider -> filler (p50 / p90)** | 2.94 / 4.47 dx | **1.37 / 2.24 dx** | <= 1.8 ✅ (previ 1.0-1.6) |
| limbo em t>32 | 101-215 | 148-175 | > 0 ✅ |
| vivas (t=32 -> 50) | 543 -> 600 | 530 -> 559 | crescendo ✅ |
| amplitude (R99=2.4) / dedos | 0.132 / 22 | **0.155 / 20** | >= 0.8x e 16-26 ✅ |
| R99 / Rmax / Rmin | 4.45 / 4.90 / 1.42 | **4.74** / 4.99 / 1.26 | >= 4.3 ✅ |
| largura EDT 0.6R / ciclo / nucleo | 5.4 dx / 0.72 / 0.32 | 6.0 dx / 0.68 / **0.24** | <= 6.5 / <= 1.0 / <= 0.38 ✅ |
| pedacos soltos / agar limpo na baia | 0.0 / 24% | 0.0 / **31%** | <= 0.5 ✅ |
| **motor `a_mar_bio_p95` (mediana t>30)** | **0.37** | **0.24** | >= 0.37 ❌ |
| agar 0.5-0.7 R99 na iter 1600 | 74 (t=25.2) | 328 (t=29.3) | >= 74 ✅ (ver ressalva) |
| massa / `max_cs` / iter/t | 202.4 / 0.4943 / 60 | 203.8 / 0.4943 / 59 | ✅ |
| filler | 10921 | 12650 (+16%) | previ +20-40% |
| `a_press_med` / picos > 4 | 2.73 / 13% | 2.99 / **27%** | nao pre-registrado ⚠ |

**1. O objetivo foi atingido e a predicao acertou:** o vao caiu 2.94 -> 1.37 dx, dentro da faixa
prevista de 1.0-1.6. A espicula da ponta e metade da do baseline.

**2. Os dois riscos pre-registrados NAO se materializaram.** Risco 1 (lei (L), engrossar e perder
relevo): a largura subiu so 0.6 dx e a **amplitude SUBIU** (0.132 -> 0.155, 1.17x). Risco 2 (licao
#104, mais pontas matam a colonizacao): as pontas ficaram em 20-22, iguais, e o limbo nunca zerou.

**3. Reprova no motor:** `a_mar_bio_p95` 0.37 -> **0.24 (-35%)**, com 41 vivas a menos no fim. Junto
vem `a_press_med` 2.73 -> 2.99 e picos > 4 em 27% das amostras contra 13% — mais material depositado
por unidade de tempo perto da frente, mesma classe do que a serie ja viu ao adensar.

**Ressalva no criterio 6:** a iteracao 1600 cai em t=29.3 aqui e em t=25.2 no P2R21 (o `dt` adaptativo
diverge), entao os 328 contra 74 NAO sao comparaveis — e o mesmo erro de fase que o proprio
CRITERIOS do P2R21 registra. O que se pode afirmar e que a colonizacao seguiu viva, pelo criterio 2.

**Veredito: REPROVADO pelo criterio 5.** `main.py` voltou a `WAKE_FREQ` = 100, bit-a-bit identico ao
snapshot do P2R21. A alavanca fica documentada: **ela resolve o vao da ponta e melhora relevo, raio e
agar limpo, ao custo de 35% do motor** — se o motor for recuperado por outra via, vale reconsiderar.
Comparacao visual em `plots/cmp_P2R21_P2R22.png`.
