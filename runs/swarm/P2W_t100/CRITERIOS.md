# P2W — P2 com kernel Wendland C2 (h = 1.92 dx)

Pre-registrado em 2026-09-14, antes de lancar. Alavanca unica (numerica): `KERNEL`
"cubic" -> "wendland_c2" e `H_FACTOR` 1.8 -> 1.92, sobre o P2. t=100, threads padrao (10).

## Diagnostico que motiva (P2_t100, licao #88)
Aglomeracao em cascata apos t~75: pares a < 0.05 dx 779 (t=70) -> 5644 (t=87), dt ~1e-4.
Pressao POSITIVA nos pares (99%) — nao e tracao. Os pares que limitam o CFL em t=87 sao
filler-viva (32/50) e viva-viva (15/50), a r ~ 0; vivem onde ha 64-121 vizinhas em 2h
contra 37 no resto. E a instabilidade de PAREAMENTO do spline cubico sob compressao (Price
2012 JCP 231:759 §2.1: cubico instavel acima de h ~ 1.5 dx). Wendland C2 e estavel ao
pareamento para qualquer numero de vizinhas (Dehnen & Aly 2012, MNRAS 425:1068; Liu §6.4).

## Predicao
- h = 1.92 dx iguala o desvio-padrao do kernel (1.43 dx). Recalculado no frame do P2 em
  t=50: sigma_a mediana 0.963 -> 0.968, frac(<0.85) 4.4% -> 3.7%; |grad cs| na frente
  mediana 1.218 -> 1.185, p95 1.617 -> 1.585 (< 3%). Ate t=50 a morfologia deve cair
  dentro do ruido entre realizacoes do P2.
- ~42 vizinhas em 2h (contra 37): ~15% mais lento por iteracao.
- Os pares coincidentes param de crescer em cascata; o dt nao colapsa depois de t=75.
- NAO resolve o corpo virando disco depois de t~60 (fisico, nao numerico).

## Criterios
1. Pareamento: pares < 0.05 dx em t~87 abaixo dos 1246 do E11 (P2: 5644), e dt de saida
   >= 0.005 ate t=100; o run chega a t=100.
2. Morfologia ate t=50 dentro do ruido do P2 (realizacoes P2 e P2_t100): bracos soltos
   (janela 35-50) <= 2.6, C5a >= 0.40, amplitude em R99=2.4 entre 0.29 e 0.34, dedos 33-40.
3. a_pressure mediana <= 3, picos > 4 <= 10%; massa sem runaway.
4. Leitura visual dendritica (§11), incluindo t=75-100.

## Resultado (2026-09-14) — numerico APROVADO; morfologia tardia segue o defeito fisico do P2

Pares < 0.05 dx contados pelo MESMO metodo nos tres runs (os valores 1246/5644 acima vieram de
outra contagem): t≈87 — P2W **692**, E11 871, P2 cubico 4872 (e 3535 -> 4872 em 0.4 s antes de
parar). Em t=100 o P2W tem 1345 (E11: 1239). `rho/rho0` p99 em t=100: P2W 7.2, E11 4.8.

| criterio | medido | |
|---|---|---|
| 1 pares em t≈87 abaixo do E11 | 692 < 871 | passa |
| 1 chega a t=100 | sim, 9291 iteracoes, 6168 s | passa |
| 1 dt >= 0.005 ate t=100 | dt MEDIO 0.0053 (t 80-90) e 0.0052 (90-100); o dt INSTANTANEO de saida fica abaixo de 0.005 em 11 dos 30 frames apos t=60 | no limite |
| 2 bracos soltos <= 2.6 | 0.8 | passa |
| 2 C5a >= 0.40 | 0.42 | passa |
| 2 amplitude em R99=2.4 | 0.310 | passa |
| 2 dedos 33-40 | 40 | passa (no limite) |
| 3 a_pressure mediana <= 3 / picos > 4 <= 10% | 2.58 / 7% (P2: 0%) | passa |
| 3 massa sem runaway | 234.5 em t=100 (P2 cubico 232.0 em t=87) | passa |
| 4 leitura visual dendritica ate t=100 | t=50 dendritico; t=75 bracos + corpo largo; t=87-100 disco com bracos curtos (relevo/Rmax 0.63 -> 0.48) | parcial — e o defeito previsto |

Fora dos criterios: agar limpo nas baias 71% contra 76-78% das duas realizacoes do P2.
Figura tardia: `tardio_E11_P2_P2W.png`.
