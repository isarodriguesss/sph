# P2R15 — decaimento de `c_s` no agar / 4: o comprimento de onda dos dedos

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R9** (a que tem folga: nucleo 0.23 Rmax
contra o alvo 0.36 da literatura, gate da base fechado em 5%, motor mais forte da serie e R99 4.74).
Alavanca unica: `AGAR_CS_LAMBDA` 0.5 -> 0.125 (parametro novo; 0.5 reproduz o anterior), ou seja o
decaimento de `c_s` no agar passa de `0.5 lambda` para `0.125 lambda`. Difusao (`D_ext` = 0.08),
producao e tudo o mais inalterados. t=50, `SEED` fixo, 10 threads.

## Hipotese (§2.4-C): o passo entre bracos escala com `L_D_ext`
A licao #95 mediu que o material ja esta certo — o ciclo de trabalho azimutal e 0.39-0.42 contra
0.32-0.46 das referencias — e que o erro esta no PASSO: 24 bracos onde Trinschek tem 9 e a PA14 12,
com passo 0.157 Rmax contra 0.31-0.42. Com o mesmo material repartido em metade dos bracos, a
largura dobra e o AR cai para a faixa da literatura sem adicionar nada.

`L_D_ext` = sqrt(`D_ext`/`lambda_agar`) = 1.03 hoje, e o passo medido em 0.6 R99 e 0.71 (mesma
ordem). Com `lambda_agar`/4, `L_D_ext` = **2.07**.

## Predicao
- Passo em 0.6 R99: 0.157 -> **~0.31 Rmax**; bracos **24 -> ~12**; largura 0.066 -> **~0.13 Rmax**;
  **AR 23 -> ~11** (topo da faixa 6.5-11.2).
- `c_s` do agar ~2x (o nivel perto de uma fonte escala com 1/sqrt(`lambda D`)): base 0.045 -> ~0.09,
  baias -> ~0.10, frente -> ~0.10. **Todos ainda abaixo do gate da colonizacao (0.3)**, entao o gate
  da base continua fechado (5% -> 10-15%) e o nucleo fica em **0.25-0.30 Rmax**.
- Halo `L`/Rmax 0.12 -> **0.20-0.24** (alvo 0.18).
- Gradiente na frente aproximadamente preservado (`c_s` da borda ×2 e `L` ×2 se cancelam):
  `a_mar_bio_p95` (mediana t>30) 2.05 -> 1.6-2.2; R99 4.3-4.7.
- `mean_cs` sobe ~2x, entao `contrast_cs` cai para ~9. **Isso NAO e criterio** (licao #75-E: sob o
  teto, `contrast_cs` = 0.494/`mean_cs` e so uma funcao da media).

## Criterios (t=50 e janela t in [35, 50])
1. **Objetivo (o teste da hipotese):** bracos em 0.6 R99 <= 16 (hoje 24) e largura >= 0.10 Rmax.
   **Se ficar acima de 16, a hipotese "passo ∝ `L_D`" esta REFUTADA** e o caminho passa a ser
   geometrico (filtro de separacao angular nos lideres).
2. Guarda do nucleo: nucleo <= 0.40 Rmax (P2R9 0.23; alvo 0.36); gate da base acima de 0.3 em
   t≈36 <= 25%.
3. Forma na faixa da licao #95: AR <= 14, area/disco 0.33-0.45, dedos (R(theta)) >= 12, R99 >= 4.2.
4. Motor: `a_mar_bio_p95` (mediana t>30) >= 1.6 (P2R9 2.05; P2R14 1.40).
5. Numerica: `max_cs` <= 0.5 (o agar nao produz, entao nao deve passar do teto); iter/t <= 1.6x o P2;
   massa sem runaway.
6. Visual (§11): classe (b), com dedos visivelmente mais largos e menos numerosos que os do P2R9.

## Resultado (2026-09-16, t=50, 1675 s) — a hipotese do comprimento de onda esta REFUTADA

| em 0.6 R99 / 0.75 R99 | P2R9 | **P2R15** | alvo |
|---|---:|---:|---:|
| **bracos** | 24 / 19 | **25 / 20** | 9-15 |
| **passo** | 0.157 / 0.248 R | **0.151 / 0.236 R** | 0.31-0.42 R |
| ciclo de trabalho | 0.36 / 0.23 | 0.41 / 0.26 | 0.32-0.46 |
| largura | 5.4 / 5.5 dx | 5.6 / 5.6 dx | 0.10-0.19 R |

**O passo NAO se moveu** (−4%, dentro do ruido), apesar de a alavanca ter engajado: `L_D_ext`
dobrou e o `c_s` do agar subiu ~3x (base 0.045 -> 0.132, baias 0.049 -> 0.133, frente 0.049 ->
0.091 — previ 2x). **Criterio 1 REPROVADO** (limite: <= 16 bracos).

**Por que, medido no mesmo dia:** cada dedo e o rastro de UMA celula viva isolada — 19-28 lideres
formando 18-21 grupos de 1.1-1.4, com o vizinho mais proximo a **14-21 dx** (nenhum tem companheiro
dentro de 3h = 5.4 dx, o alcance de toda forca do modelo). Nao ha onda para alongar: a contagem de
dedos e a contagem de celulas na frente, e nenhuma alavanca de campo a altera.

**O que a alavanca ENTREGOU (e e util):**

| | P2R9 | **P2R15** | P2R10 (conduz) | alvo |
|---|---:|---:|---:|---:|
| `c_s` do corpo / teto | 0.16 (artefato) | **0.36** | 0.47 | 0.99 |
| `c_s` das baias / teto | 0.11 | **0.25** | 0.21 | 0.95 |
| halo `L` / Rmax | 0.12 | **0.15** | 0.12 | 0.18 |
| nucleo / Rmax | 0.23 | **0.27** | 0.38 | 0.36 |
| area/disco | 0.33 | **0.38** | 0.39 | 0.35-0.39 |
| AR | 23.2 | 18.4 | 15.4 | 6.5-11.2 |
| R99 / `a_mar_bio_p95` | 4.74 / 2.05 | 4.57 / 1.67 | 4.31 / 1.72 | — |
| gate da base (t=50) / vivas / iter/t | 6% / 622 / 62 | 15% / 783 / 69 | 13% / 668 / 81 | — |

2. Guarda do nucleo: **aprovado** (0.27 <= 0.40; gate 15% <= 25%).
3. Forma: aprovado em area/disco (0.38, dentro da faixa), dedos e R99; **reprovado em AR** (18.4 > 14).
4. Motor: aprovado no limite (1.67 >= 1.60).
5. Numerica: aprovado (`max_cs` 0.494 o run inteiro, iter/t 69 = 1.2x o P2, massa 205).

**Leitura:** `AGAR_CS_LAMBDA` nao e uma alavanca de comprimento de onda — e uma alavanca de CAMPO, e
barata: sobe o `c_s` do corpo 2.2x e o das baias 2.3x **sem conduzir nada pelo filler**, mantendo o
nucleo em 0.27 (contra 0.38 do P2R10 e 0.42 do P2R14) e o R99 em 4.57 (o maior dos runs com campo).
No ranking da licao #95 ele fica em 3o (0.407), atras do P2R14 (0.262) e do P2R10 (0.308), mas e o
unico que melhora o campo mantendo o nucleo ABAIXO do alvo — ou seja, e a base com folga para
receber as outras alavancas.
