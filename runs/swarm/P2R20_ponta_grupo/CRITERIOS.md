# P2R20 — calota so na PONTA de cada grupo de lideres (`RASTRO_PONTA_LINK` = 6 dx)

Pre-registrado em 2026-09-17, ANTES da rodada completa. Base: **P2R16**. Alavancas: `RASTRO_PONTA` 0 ->
1.0 e `RASTRO_PONTA_LINK` 0 -> 6.0 (a segunda so tem sentido com a primeira; e UMA mudanca de desenho:
"calota so no lider mais externo de cada grupo"). t=50, `SEED` fixo, 10 threads.

## Motivacao (licao #104-D)
No P2R19 a calota em todo lider poupou o agar lateral nos primeiros 5 dx atras da ponta; co-lideres
nessa zona nao ficaram cercadas de filler e nao pararam -> +30% de pontas -> baias a meia altura fecharam
-> colonizacao morreu. Aplicar a calota so no lider mais externo de cada grupo (grupo = lideres a
<= 6 dx, o mesmo raio do diagnostico) deve manter o corte plano em volta das co-lideres.

## Teste de efeito (t=18, regra #102) — passou no mecanismo, parcial nos sinais
| co-lideres t≈10->17 | 0-3 dx atras | 3-5 dx atras | >5 dx |
|---|---:|---:|---:|
| P2R16 | 82% | **6%** | 0% |
| **LINK=6** | 92% | **12%** | 0% |
| P2R19 | 94% | **83%** | 0% |

| t≈17 | P2R16 | LINK=6 | P2R19 |
|---|---:|---:|---:|
| lideres externos | 20 | 24 | 29 |
| agar 0.5-0.7 R99 | 28 | 6 | 0 |
| vivas | 313 | 310 | 266 |

Estado diverge do P2R19. Na zona controlada pela calota (3-5 dx) o comportamento volta ao do P2R16. Os
+4 lideres estao na faixa 0-3 dx, que sobrevive em qualquer config.

## Predicao
- **Colonizacao viva:** limbo > 0 e vivas crescendo apos t=32 (P2R19: limbo 0, vivas congeladas 525).
- Lideres externos em t>=20: **20-24**. Agar na faixa 0.5-0.7 R99 em iter 1600: **>= 100** (P2R16 212,
  P2R19 21).
- **Ponta:** apos t≈17 quase todo lider e unico no grupo e recebe a calota -> perfil parecido com o
  P2R19: u=2 **3.5-5.5 dx** (em fase, janela iter 1600-2800; P2R16 7.1, P2R19 3.3).
- **Risco aberto:** o afinamento do corpo do P2R19 (9.7 -> 8.0) nao tem mecanismo; se vier da calota em
  si e nao do excesso de pontas, repete aqui.

## Criterios (frames em fase, iter 1600-2800; mais t=50)
1. **Objetivo:** u=2 <= 6.0 dx **e** colonizacao viva (limbo > 0 em todos os frames com t > 32 e vivas
   crescendo entre t=32 e t=50).
2. Falseabilidade: corpo em u=5 >= 9.0 dx (P2R16 9.7). Se reprovar, o afinamento e da calota, nao das
   pontas extras.
3. Guarda: lideres externos em t>=20 <= 24; bracos em 0.5R dentro de ±2 do P2R16 (17); R99 >= 4.4;
   nucleo <= 0.40 Rmax; pedacos soltos <= 0.5; agar 0.5-0.7 R99 (iter 1600) >= 100.
4. Motor: `a_mar_bio_p95` (mediana t>30) >= 0.40.
5. Numerica: `max_cs` <= 0.5; iter/t <= 90; massa sem runaway.
6. Visual (§11), no mesmo renderizador e em frame em fase.

## Resultado (2026-09-17, t=50) — REPROVADO na falseabilidade (e por 0.01 em R99), mas resolve a colonizacao e isola a causa do afinamento

| | P2R16 | P2R19 | **P2R20** | criterio |
|---|---:|---:|---:|---|
| limbo (iter 1600/2000/2400/2800) | 91/143/145/142 | 2/0/0/0 | **75/36/69/100** | > 0 apos t=32 ✅ |
| vivas (iter 1600 -> 2800) | 512 -> 686 | 523 -> 525 | **508 -> 612** | crescendo ✅ |
| **ponta u=2** (em fase, iter 1600-2800) | 7.1 | 3.3 | **3.6** | <= 6.0 ✅ |
| u=3 | 9.2 | 4.8 | 6.4 | — |
| **corpo u=5** | 9.7 | 8.0 | **8.3** | >= 9.0 ❌ |
| largura EDT 0.75R | 5.0 | 4.5 | 4.1 | — |
| lideres externos (t >= 20) | 16-20 | 23-27 | **21-22** | <= 24 ✅ |
| bracos em 0.5R (iter 2400/2600/2800) | 16/17/17 | 15/19/23 | **16/16/19** | 17 ± 2 ✅ |
| agar 0.5-0.7 R99 (iter 1600) | 212 | 21 | **102** | >= 100 ✅ |
| R99 | 4.68 | 4.57 | **4.39** | >= 4.4 ❌ (por 0.01) |
| nucleo / pedacos soltos | 0.29 / 0 | 0.32 / 0 | **0.32 / 0** | ✅ |
| `a_mar_bio_p95` (t>30) | 0.44 | 0.35 | **0.79** | >= 0.40 ✅ |
| iter/t / `max_cs` / massa / picos p>4 | 57 / 0.494 / 203.2 / 27% | — | **58 / 0.494 / 203.0 / 13%** | ✅ |

1. **Objetivo: APROVADO.** Ponta afila (3.6) e a colonizacao SOBREVIVE (limbo > 0 em todos os
   frames, vivas crescendo). Separar a calota do corte das co-lideres funcionou.
2. **Falseabilidade: REPROVADO**, e com a leitura que o criterio previa: o corpo afina igual com e sem
   pontas extras (8.0 no P2R19, 8.3 aqui) -> **o afinamento e da CALOTA, nao das pontas**.
3. Guardas aprovadas, exceto R99 4.39 (limite 4.4).
4. Motor 0.79, **o maior da serie W=5** (P2R16 0.44).
6. Visual (`plots/cmp_P2R16_P2R20.png`, iter 2600): pontas pontiagudas, bracos em fuso 3.8 dx (P2R16
   4.3), mesma contagem (16 x 17), sem fragmentos.

## Causa do afinamento — confirmada por teste geometrico puro

Cada segmento converte largura plena `2w` so entre `p0` e `p1`; atras de `p0` cobre apenas um disco de
raio `w`. A zona que a calota poupou numa chamada fica logo atras do `p0` da chamada seguinte e **nunca
e reposta**. Minha conta de "afinamento transitorio" supunha o contrario, e estava errada.

Lider em linha reta, uniao das regioes convertidas:

| | passo 1.5 dx | passo 2.9 dx | passo 4.0 dx | ponta u=1/2/3/5 |
|---|---:|---:|---:|---|
| corte plano | 10.00 | 10.00 | 10.00 | 10/10/10/10 |
| **calota** | **8.59** | **9.25** (min 8.65) | 9.45 | 6.0/8.0/9.2/9.1 |
| **calota + recobre `a` atras de p0** | **10.00** | **10.00** | **10.00** | **6.0/8.0/9.2/10.0** |

O minimo de 8.65 = `sqrt(3)/2 * 2w` (ponto a `w/2` atras, onde o disco e a calota se cruzam). **Passo
menor afina mais** — explica o afinamento maior na simulacao (passos reais variam). Estender o
`clip` de `t` ate `-a/L` so para os lideres com calota repoe a largura plena e deixa a calota apenas na
ponta corrente, como era a intencao original.
