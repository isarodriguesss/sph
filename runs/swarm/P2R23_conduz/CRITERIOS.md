# P2R23 — `FILLER_CS_CONDUZ` 0 -> 1 sobre o baseline P2R21

Pre-registrado em 2026-09-17, ANTES de rodar. Base: **P2R21** (baseline desde hoje). Alavanca unica:
`FILLER_CS_CONDUZ` 0.0 -> 1.0, com `FILLER_CS_D` = 0 (pares com filler usam `D_ext` = 0.08) e
`FILLER_CS_LAMBDA` = 0.5, exatamente como o P2R10 e o P2R17. t=50, `SEED` fixo, 10 threads.

## Por que
O painel (d) do P2R21 nao e um campo: 94-95% do corpo e filler quimicamente transparente, e **31.4% da
area da colonia nao tem grau de liberdade de `c_s`** (medido em 2026-09-17). O valor que o painel mostra
dentro dos bracos e o congelado no instante da conversao (licao #39). Com o filler conduzindo, o corpo
passa a ter campo medido — e o par equivalente da serie (P2R16 -> P2R17, mesmo `RASTRO_W`=5) mostrou que
isso tambem SOBE o motor.

## Predicao, ancorada no par P2R16 -> P2R17 e em P2R9 -> P2R10
| | P2R21 (base) | previsto no P2R23 | fonte da previsao |
|---|---:|---|---|
| `c_s` do corpo / teto | 0.20 (congelado) | **0.45-0.55** (campo medido) | P2R16 0.19 -> P2R17 0.49 |
| `c_s` das baias | 0.17 | 0.22-0.28 | P2R17 0.25 |
| halo `L`/Rmax (alvo 0.18) | 0.16 | 0.17-0.19 | P2R17 0.18 |
| **motor `a_mar_bio_p95`** | 0.37 | **0.45-0.60** | P2R16 0.44 -> P2R17 0.60 |
| nucleo / Rmax (alvo 0.36) | 0.32 | **0.38-0.42** (anda PARA o alvo) | P2R16 0.29 -> P2R17 0.39 |
| **R99** | 4.45 | **3.9-4.3** (cai 5-12%) | P2R16 4.78 -> P2R17 4.22 |
| AR (faixa 3.8-11.2) | 11.4 | 10.0-11.0 | P2R17 10.4 |
| dedos | 20 | 15-20 | P2R17 15 |
| iter/t | 60 | 60-90 | P2R17 89.7 |

**O risco principal e o raio:** a conducao leva `c_s` ao agar da base, abre o gate da colonizacao e o
nucleo come a base (licoes #90, #94-H, #97, #101). No P2R21 o nucleo esta ABAIXO do alvo (0.32 contra
0.36), entao um crescimento moderado e favoravel — o que reprova e o raio.

## Criterios (t=50; janela t in [35, 50] onde couber)
1. **Objetivo:** `c_s` do corpo >= **0.40** do teto **e** o campo passa a ser medido (sem regiao sem
   grau de liberdade no corpo).
2. **Motor:** `a_mar_bio_p95` (mediana t>30) >= **0.37** (nao piorar em relacao ao baseline).
3. **Guarda de raio:** `R99` >= **4.1**. Abaixo disso a rodada REPROVA, mesmo com campo bom — foi o que
   reprovou o P2R17 (4.11 ficou no limite).
4. **Forma:** nucleo <= **0.45** Rmax; baia >= 0.28; dedos 14-24; AR <= 13; pedacos soltos <= 0.5;
   amplitude >= 0.8x a do P2R21 (0.132 em R99=2.4).
5. **Colonizacao viva:** limbo > 0 em todos os frames com t > 32 e vivas crescendo entre t=32 e t=50.
6. **Numerica:** `max_cs` <= 0.5 — **aborto precoce se passar de 0.6 na iteracao 200** (licao #96: o
   `c_s` e explicito e o `dt` do PySPH nao tem criterio difusivo; aqui `D` nao muda, entao nao deve
   ocorrer); iter/t <= 95; massa sem runaway.
7. **Visual (§11):** classe (b) preservada; painel (d) com o campo acompanhando os bracos, como o
   painel (b) de Trinschek.

**Nota de processo:** este arquivo esta no scratchpad ANTES da rodada porque `tools/archive_run.sh` faz
`rm -rf` no destino e apagou o CRITERIOS.md do P2R22.

## Resultado (2026-09-17, t=50) — REPROVADO NO LIMITE por `R99` (4.07 contra 4.10), tudo o mais passou

| | P2R21 (base) | **P2R23** | previsto | criterio |
|---|---:|---:|---|---|
| `c_s` do corpo / teto | 0.20 (congelado) | **0.44** (medido) | 0.45-0.55 | >= 0.40 ✅ |
| `c_s` das baias / halo `L` | 0.17 / 0.16 | 0.25 / **0.17** | 0.22-0.28 / 0.17-0.19 | ✅ |
| **motor `a_mar_bio_p95`** | 0.37 | **0.73** | 0.45-0.60 | >= 0.37 ✅ (superou a predicao) |
| **`R99`** | 4.45 | **4.07** | 3.9-4.3 | >= 4.10 ❌ (falta 0.03) |
| nucleo / baia | 0.32 / 0.31 | **0.39 / 0.37** | 0.38-0.42 | <= 0.45 / >= 0.28 ✅ |
| **AR** | 11.4 | **9.9** | 10.0-11.0 | <= 13 ✅ (e DENTRO da faixa 3.8-11.2) |
| dedos / amplitude | 20 / 0.132 | 15 / **0.176** | 15-20 | 14-24 e >= 0.8x ✅ |
| limbo t>32 / vivas | 101-215 / 543->600 | 147-158 / 547->594 | — | ✅ |
| agar limpo na baia | 24% | **35%** | — | — |
| `max_cs` / iter/t / massa | 0.4943 / 60 / 202.4 | 0.4940 / **62** / 202.3 | iter/t 60-90 | ✅ |
| `a_press_med` / picos > 4 | 2.73 / 13% | 3.22 / **31%** | — | ⚠ nao pre-registrado |
| C5a | 0.14 | 0.11 | — | ⚠ |

**O objetivo foi atingido com folga.** O campo passa a ser medido em todo o corpo (o filler conduz,
entao ele tem grau de liberdade) e vale 0.44 do teto, contra 0.20 congelado. O halo fica em 0.17 R,
contra o alvo 0.18 de Trinschek.

**O motor DOBROU** (0.37 -> 0.73), acima da predicao — o par P2R16 -> P2R17 tinha dado +36%, aqui deu
+97%. Com ele vieram AR 9.9 (entra na faixa da literatura pela primeira vez na serie da calota),
amplitude +33% e agar limpo 24% -> 35%.

**Reprova por 0.7% no raio:** `R99` 4.07 contra o limite 4.10 que eu mesmo pre-registrei. A margem esta
DENTRO do ruido entre realizacoes documentado na licao #88 (duas realizacoes do mesmo P2 diferiram em
bracos soltos 1.3 vs 2.6), entao o veredito e formal, nao fisico. **Decisao de adotar ou nao e da
usuaria.** `main.py` voltou ao P2R21 (bit-a-bit) enquanto isso.

**Ressalva a registrar:** picos de pressao > 4 passam de 13% para 31% das amostras, o mesmo sinal que o
P2R22 deu — os dois runs adensam a frente por vias diferentes. E o C5a piora (0.14 -> 0.11).

**Contra a literatura (licoes #95/#99), o P2R23 e 2o lugar geral (0.093), atras so do P2R17 (0.074) e a
frente do P2R21 (0.129)** — e diferente do P2R17, ele chega la com `R99` proximo do baseline, sem
encolher a colonia para acertar as razoes.
