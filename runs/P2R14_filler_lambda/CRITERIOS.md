# P2R14 — a matriz PROTEGE o surfactante (decaimento do filler 0.5 lambda -> 0.25 lambda)

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R10**. Alavanca unica: o sumidouro de
`c_s` no filler passa de `0.5 lambda` para `0.25 lambda` (`FILLER_CS_LAMBDA`, parametro novo;
0.5 = comportamento anterior). Conducao inalterada (`FILLER_CS_D` = 0, ou seja `D_ext` = 0.08 nos
pares com filler, como no P2R10). t=50, `SEED` fixo, 10 threads.

## Por que esta alavanca, depois do P2R13
O P2R13 (conducao interna `D` = 0.25 so entre pares filler-filler) **divergiu numericamente**: o
`c_s` maximo saltou de 0.494 para 5.05 ja na iteracao 200. Causa medida: o limite explicito de
difusao `dt <= h^2/(4D)` vale 0.0094 s com `D` = 0.25, e o `dt` do baseline e 0.017-0.022 s. Com
`D_ext` = 0.08 o limite e 0.029 s, e o baseline passa com folga de 1.4x. **Subir `D` acima de ~0.10
exige limitar o `dt` (custo ~3x em iteracoes).**

O comprimento de conducao ao longo do braco e `L = sqrt(D / lambda_filler)`. O mesmo alongamento
pode vir do DENOMINADOR, que nao tem custo numerico nenhum (o sumidouro e local, `lambda*dt` ~ 2e-4).
Fundamentacao: e o mesmo argumento do Pass J e do §7 — a matriz de EPS protege o ramnolipideo da
degradacao que ele sofre no agar livre. Aqui o filler deixa de decair como agar e passa a decair
como material protegido.

## Predicao
- `L` = sqrt(0.08/0.0375) = **1.46** contra 1.03 (1.4x).
- `c_s` do filler por faixa (0.2-0.4 / 0.4-0.6 / 0.6-0.8 R99): 0.41 / 0.22 / 0.12 (P2R10) ->
  **0.42-0.46 / 0.28-0.35 / 0.16-0.22**.
- Corpo (mediana sobre a colonia, `tools/lit_rank.py`): 0.47 -> **0.55-0.65** do teto (alvo 0.99);
  baias 0.21 -> **0.25-0.30** (alvo 0.95).
- Agar da base (0.2-0.4 R99) 0.29 -> 0.31-0.36; nucleo 0.38 -> **0.40-0.45** Rmax.
- `max_cs` CRAVADO em 0.494 (o teto) o run inteiro — sem instabilidade difusiva.

## Criterios (t=50 e janela t in [35, 50], contra P2R10 e a faixa da literatura da licao #95)
0. **Aborto precoce:** se `max_cs` > 0.6 na iteracao 200, o run e numericamente invalido (repete o
   P2R13) e e interrompido.
1. **Objetivo:** `c_s` do corpo >= 0.55 do teto (P2R10 0.47) e nas baias >= 0.25 (P2R10 0.21).
2. **Guarda do nucleo:** nucleo <= 0.45 Rmax (P2R10 0.38; alvo da literatura 0.36).
3. Forma na faixa da licao #95: area/disco 0.33-0.42, dedos >= 15, AR <= 20, R99 >= 4.1.
4. Quimica: `contrast_cs` >= 12.9; frac(`c_s` > 0.45) nas vivas >= 0.6.
5. Numerica: iter/t <= 1.6x o P2, massa sem runaway.
6. Visual (§11): classe (b); painel (d) com o corpo mais uniforme que o do P2R10.

## Resultado (2026-09-16, t=50, 3398 s) — janela t in [35, 50], reguas da licao #95

| | P2R9 (sem conducao) | P2R10 (conduz) | **P2R14 (protegido)** | alvo da literatura |
|---|---:|---:|---:|---:|
| `c_s` do filler 0.2-0.4 / **0.4-0.6** / 0.6-0.8 R99 | congelado | 0.41 / **0.22** / 0.12 | 0.45 / **0.32** / 0.18 | — |
| **corpo / teto** (`lit_rank`) | 0.16 (artefato) | 0.47 | **0.63** | 0.99 |
| **baias / teto** | 0.11 | 0.21 | **0.28** | 0.95 |
| halo `L` / Rmax | 0.12 | 0.12 | 0.14 | 0.18 |
| **nucleo / Rmax** | 0.23 | 0.38 | **0.42** | 0.36 |
| dedos / area-disco / AR | 20 / 0.33 / 23.2 | 19 / 0.39 / 15.4 | 18 / **0.43** / **12.7** | 9-15 / 0.35-0.39 / 6.5-11 |
| R99 / dR/dt | 4.74 / 0.096 | 4.31 / 0.089 | **3.89 / 0.057** | — |
| gate da base acima de 0.3 (t≈36 / t=50) | 5% / 6% | 67% / 13% | **95% / 80%** | — |
| vivas / massa / iter/t | 622 / 204 / 62 | 668 / 206 / 81 | **783 / 210 / 137** | — |
| `a_mar_bio_p95` (mediana t>30) / `contrast_cs` | 2.05 / 18.4 | 1.72 / 15.1 | **1.40** / 14.4 | — |

**Predicao:** acertou o campo — filler em 0.4-0.6 R99 previsto 0.28-0.35, medido **0.32**; corpo
previsto 0.55-0.65, medido **0.63**; baias previstas 0.25-0.30, medidas 0.28; nucleo previsto
0.40-0.45, medido 0.42. Errou o custo: nao previ a queda de R99 (3.89) nem o salto de iter/t.

0. Aborto precoce: **passou** — `max_cs` cravado em 0.494 o run inteiro (o P2R13 tinha 5.05 na
   iteracao 200). O sumidouro e local, entao a alavanca nao tem custo de estabilidade.
1. Objetivo: **APROVADO** — corpo 0.63 (limite 0.55) e baias 0.28 (limite 0.25). E o melhor campo
   ja medido no projeto, e o painel (d) fica com a cara do painel (b) de Trinschek
   (`plots/fig_P2R14_t50.png`).
2. Guarda do nucleo: aprovado no limite (0.42 <= 0.45), mas ACIMA do alvo da literatura (0.36) e do
   P2R10 (0.38).
3. Forma: **REPROVADO** em R99 (3.89 < 4.1) e por pouco em area/disco (0.43 > 0.42). Dedos 18 e
   **AR 12.7 — o melhor da serie** (P2R10 15.4), mas em parte porque o raio encolheu.
4. Quimica: aprovado (0.67 e 14.4).
5. Numerica: **REPROVADO** — iter/t 137, 2.2x o P2 (limite 1.6x), 3398 s de parede.
6. Visual: classe (b), nucleo visivelmente maior, bracos mais curtos.

**A TROCA E MONOTONICA, e e o achado da serie.** Nos tres passos (sem conducao -> conduz ->
protegido) o campo dentro da colonia sobe (0.16 -> 0.47 -> 0.63) e TUDO que depende da frente cai
junto: `a_mar_bio_p95` 2.05 -> 1.72 -> 1.40, R99 4.74 -> 4.31 -> 3.89, dR/dt 0.096 -> 0.089 ->
0.057, nucleo 0.23 -> 0.38 -> 0.42, gate da base 5% -> 67% -> 95%. **O mecanismo e um so:** todo
`c_s` que se poe dentro da colonia vaza para o agar (o agar em volta dos bracos vai de 0.049 a 0.243
em 0.4-0.6 R99), e o agar mais rico (a) enfraquece o gradiente na frente e (b) abre o gate da
colonizacao, que engorda o nucleo (licao #90). **A colonizacao le exatamente o campo que estamos
tentando levantar** — esse acoplamento, e nao a quimica, e o que limita a rota.
