# Do C4 ao E5 — o que mudou no baseline

> **Estabelecido 2026-09-01.** O `runs/C4` foi o baseline do projeto de 2026-08-07 até
> aqui. A série E o substitui. Este documento lista **o que mudou, por quê, e o efeito
> medido** de cada mudança — inclusive uma que não é calibração, mas correção de bug.
>
> Configuração final: `k_src=0.3`, `HILL_K=0.25`, `sigma=11.1`, `k_col=0.03`, com a
> `BiomassEOS` em `Group` próprio. Rodadas: `runs/E5_hillK025` (t=50) e `runs/E11_t100`
> (mesma config, t=100).

**Índice.** [Resumo](#resumo-das-mudancas) · [O bug da EOS](#1-correcao-de-bug--a-eos-nunca-repelia) ·
[Nutriente](#2-k_src--fonte-de-nutriente) · [Hill K](#3-hill_k--o-joelho-do-quorum-sensing) ·
[Colonização](#4-k_col--recrutamento) · [Parametrizações](#5-parametrizacoes-sem-efeito-fisico) ·
[Comparação](#comparacao-medida) · [Em aberto](#o-que-permanece-em-aberto)

---

## Resumo das mudanças

| # | mudança | de → para | natureza | rodada |
|---|---|---|---|---|
| 1 | `BiomassEOS` em `Group` próprio | — | **correção de bug** | E1 |
| 2 | `k_src` (fonte de nutriente) | 0.0 → **0.3** | calibração | E1 |
| 3 | `HILL_K` (joelho do quorum sensing) | 0.1 → **0.25** | calibração | E2/E5 |
| 4 | `sigma` (recalibrado para o `HILL_K` novo) | 10.0 → **11.1** | consequência de #3 | E5 |
| 5 | `k_col` (colonização) | 0.0 → **0.03** | calibração | E3 |
| 6 | `SEED_MODE`, `LAMBDA_BIO_RATIO` | hardcoded → parâmetro | sem efeito físico | — |
| 7 | `use_floor`, `use_bridge`, `use_matrix`, `use_promo` | — → `False` | código preservado desligado | — |

---

## 1. Correção de bug — a EOS nunca repelia

**A mudança mais importante não é calibração.** A `BiomassEOS` estava no mesmo `Group`
que o `SummationDensity`, então lia `d_rho` **durante** a acumulação do somatório — que
começa em zero no `initialize` — e não o valor final. O ramo repulsivo `if ratio > 1`
**nunca disparava**.

Medido nos estados finais:

| | `p > 0` | `\|p\|` máximo | partículas com `rho > 1` |
|---|---:|---:|---:|
| **C4** | **0 de 70982** | 1.575e-3 | 65961 |
| E5 | 1887 | **1.645e-1** (104×) | 66606 |
| E11 (t=100) | 6903 | 1.002e-1 | 71380 |

No C4, **65961 partículas estavam comprimidas (`rho > 1`) e nenhuma empurrava de volta**.
O `|p|` máximo bate exatamente em `B_tension·0.3 = 1.575e-3`, o teto do ramo **atrativo** —
ou seja, a única pressão que existia era coesão.

**Atenção ao ler o histórico:** a coluna `a_pressure` do log é uma **estimativa**
(`B·excess²`), calculada à parte, e reportava ~3.0 no C4. Ela não vinha do campo `p` real.
Toda conclusão anterior sobre "pressão dentro do orçamento §8" no C4 foi tirada dessa
estimativa, não da EOS efetiva.

É a mesma correção que o `KernelSum` já tinha recebido, pelo mesmo motivo (ver §12 do
CLAUDE.md, Pass N v2.4).

---

## 2. `k_src` — fonte de nutriente

`0.0 → 0.3`. Adiciona `dc_n/dt += k_src·(1−c_n)` ([T2] Srinivasan: swarming é regime
*nutrient-rich*, ao contrário do biofilme).

**Calibração:** o equilíbrio local é `c_n = k_src/(k_src + k_n·rho_b)`. Com `k_src=0.1` o
valor nos braços seria 0.40 — exatamente o zero do smoothstep `[0.4, 0.8]` do gate de
crescimento. `0.3` dá 0.67 (gate 0.74).

**Efeito medido:** `min_c_n` nunca cai abaixo de **0.376** (no C4 vai a zero em t≈55);
`a_mar` na frente **+23%**; vazio >1.5dx 0.16% → **0.02%**; `frac(sigma_a<0.85)` 7.4% →
2.6%. E o principal: **elimina o congelamento**. No C4 em t=96, `mean_v` cai para 0.00002
(25× menos que em t=50); no E11 em t=100 é **0.00125**, *acima* do valor em t=50.

**Efeito colateral que desbloqueou o resto:** o gate de nutriente da `BiomassColonization`
estava fechado no C4 — só **36** partículas elegíveis. Com `k_src=0.3` são **712** (lição
#71-A). Toda a série J havia calibrado `k_col` contra um gate quase fechado.

---

## 3. `HILL_K` — o joelho do quorum sensing

`0.1 → 0.25`, em `qs = rho_b²/(rho_b² + K²)`. Era `0.1` desde o início do projeto, o que
punha o joelho em **10% da densidade de saturação**: `qs(0.1) = 0.5`, ou seja meia
produção de surfactante com um décimo da biomassa. É o que tornava a banda sub-quórum
quimicamente barulhenta e o que afogou as rotas de recrutamento (lições #48, #52, #53, #59).

**Inversão de sinal documentada (lição #71-B):** sobre o C4 puro, `K=0.3` **reprova** em
AR (5.92 → 4.97, lição #70). Sobre `k_src=0.3` o mesmo lever **sobe** o AR. Medido em duas
rodadas independentes, e nas duas a previsão foi de queda: E2 (`K`=0.3) previsto 5.1–5.4,
medido **6.31**; E5 (`K`=0.25) previsto 6.0–6.5, medido **7.61**.

**Interpolar entre rodadas isoladas é inválido quando os mecanismos compõem** — silenciar
a banda sub-quórum só custa gradiente na zona de transição quando essa zona está faminta;
com nutriente sustentado, não está.

---

## 4. `sigma` — recalibração

`10.0 → 11.1`. Não é uma alavanca independente: com `K=0.25` o `qs` cai, e `11.1` é o
valor que preserva a produção da biomassa **madura** (`rho_b > 0.5`) constante.

---

## 5. `k_col` — recrutamento

`0.0 → 0.03`. Reativa a `BiomassColonization`, o único termo **aditivo** em `rho_b`
(`BiomassGrowth` é multiplicativo, então `rho_b = 0` é estado absorvente e nenhum gate o
ressuscita — lição #49).

O valor `0.03` vem da varredura da lição #49: elimina 100% dos zeros absolutos custando
1.2% do motor, enquanto `0.1` custa 49% e `0.3` custa 81%.

**Efeito medido:** VIVAS 167 → **276** em t=50, e **1209** em t=100.

---

## 6. Parametrizações sem efeito físico

- `SEED_MODE = 8` — o modo azimutal do inóculo era `cos(8·theta)` hardcoded; virou
  parâmetro. O M5 testou `16` e não se sustentou (dedos finos mas ocos, ágar cercado 417
  contra 181).
- `LAMBDA_BIO_RATIO = 2.0` — a razão do decaimento de `cs` no biofilme, também hardcoded.

---

## 7. Código preservado desligado

Conforme §10 do CLAUDE.md, mecanismos reprovados ficam no código com a flag em `False`,
nunca removidos: `use_floor` (piso, séries E7–E10 — movido para a renderização, lição
#72-K), `use_bridge` (ponte, lição #68), `use_matrix` (série M), `use_promo` (série P),
`AGAR_FADE`, `AGAR_DRAG_RATIO`, `FILLER_RHO_B_FLOOR`.

---

## Comparação medida

| | C4 (t=50) | C4 (t=96) | **E5 (t=50)** | **E11 (t=100)** |
|---|---:|---:|---:|---:|
| **AR** | 5.92 | 5.08 | **7.61** | 4.16 |
| dedos | 20 | 14 | **24** | 22 |
| `R99` | 4.62 | 6.53 | 4.15 | 5.37 |
| VIVAS | 167 | 167 | 276 | **1209** |
| biomassa real | 0.362 | 0.358 | 0.472 | **0.983** |
| **`mean_v`** | 0.00049 | **0.00002** | 0.00088 | **0.00125** |
| `a_mar` | 11.35 | 8.09 | **12.31** | 10.68 |
| `a_press` mediana / picos>4 | 3.00 / 13% | 3.00 / 10% | **2.58 / 7%** | 2.85 / 20% |
| `contrast_cs` | 12.9 | 22.8 | 12.7 | 10.9 |
| `frac(cs>0.45)` | **9.7%** | **6.1%** | 52.4% | 52.4% |
| `min_c_n` | 0.000 | 0.000 | **0.376** | **0.376** |
| **`C5a`** @1.05dx | 0.12 | — | 0.20 | **0.36** |
| buraco cercado (dx²) | 3357 | — | **2751** | — |
| janela útil | t ≲ 55 | — | — | **t > 100** |

**Instante de referência morfológica: t=50.** Entre t=50 e t=100 o AR cai de 7.84 para
4.16 — a biomassa dobra e os dedos engordam (lei (L), lição #67-L). O t=100 serve para
medir C5 e provar que a janela existe, não como a figura da morfologia.

---

## O que permanece em aberto

**`frac(cs>0.45)` = 52.4%**, contra 9.7% do C4 — é o único eixo em que o C4 ganha, e é a
única reprovação do E5. Significa `cs` saturado sobre metade da biomassa, que é o alerta
da lição #43 sobre gradiente achatado. Ainda não cobra (`a_mar` é *maior* no E5), mas é
provavelmente o que bloqueia a rota A do
[PLANO_C5_CONTINUIDADE](PLANO_C5_CONTINUIDADE.md), e o problema aberto da lição #64.

**C5 (continuidade da expansão)** continua violado — `C5a` = 0.20 em t=50. É o objeto do
plano acima, e o E5 é o ponto de partida escolhido (custo da rota C: 43% do material,
contra 54% no C4).
