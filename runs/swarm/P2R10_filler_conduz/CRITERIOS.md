# P2R10 — P2R9 + filler conduz surfactante (opcao A)

Pre-registrado em 2026-09-15, antes de lancar. Alavanca unica contra o P2R9:
`FILLER_CS_CONDUZ` 0 -> 1. No `SurfactantEquation`, o filler entra na difusao (fonte e destino)
tratado como agar (`rho_b` contado como 0 no D bi-escala -> `D_ext`) e decai a 0.5 lambda (como o
agar); NAO produz. Marangoni e flagelo continuam sem o filler (o gradiente das forcas usa so
nao-filler). Operador inalterado (Laplaciano de Brookshaw ja usado; so muda o conjunto de pares).
t=50, `SEED` fixo, 10 threads. O `equations.py` do run difere do P2R9 so pelo parametro novo.

## Diagnostico que motiva (P2R8/P2R9, painel (d) contra Trinschek (b))
Na referencia Γ ≈ Γmax em toda a colonia (dedos e baias) e o gradiente fica so na borda; no nosso
painel so o nucleo chega a 0.5. Causa medida: 95% da colonia e filler, que nao produz, nao recebe
e nao conduz `c_s` (licao #39). Fora do nucleo so ~40 lideres produzem; o agar em volta dos bracos
fica em `c_s` ≈ 0.05 (P2R9: 0.045-0.059 em 0.2-1.0 R99). Precedente: M3 (licao #67-D), matriz
conduzindo `c_s` eliminou o penhasco e restaurou o motor; M4 (filler sentindo Marangoni) espalhou
os bracos — por isso aqui o filler NAO sente Marangoni.

## Predicao (§3.3.6: producao inalterada; so muda transporte/decaimento no filler)
- `L_D` no filler = sqrt(D_ext / 0.5 lambda) = 1.03: o `c_s` do nucleo (0.48) difunde pela base.
- `c_s` do filler (agora fisico), mediana por faixa: 0.2-0.4 R99 ~0.15-0.30; 0.4-0.8 R99
  ~0.05-0.15; 0.8-1.0 R99 ~0.08-0.2 (pontas produzem). Perfil DECRESCENTE do nucleo para fora —
  pico dentro do corpo, sem patologia da §3.3.4.
- NAO satura como Trinschek: bracos bem abaixo de 0.45 (so o B, filler produzindo, faria isso).
- Agar em volta dos bracos sobe um pouco (fluxo do nucleo atravessa a base): 0.05 -> 0.05-0.10.
- Morfologia e motor iguais ao P2R9 dentro do ruido (forcas nao usam o filler).

## Criterios (t=50 e janela t in [35, 50], contra P2R9)
1. **Objetivo:** `c_s` do filler com mediana decrescente do nucleo para fora e >= 0.15 em
   0.2-0.4 R99; painel (d) com `--cs-filler` sem descontinuidade entre nucleo e bracos.
2. Forma: pedacos soltos <= 0.5; largura p50 em 0.8 R99 >= 4.5 dx; tortuosidade <= 1.08;
   bracos em 0.8 R99 >= 17; `R99`(t=50) >= 0.9x o P2R9 (4.27).
3. Quimica: `contrast_cs` t=50 >= 0.7x o P2R9 (12.9); frac(`c_s` > 0.45) nas vivas >= 0.8.
4. Numerica: iter/t <= 1.5x o P2, massa sem runaway.
5. Visual (§11): classe (b); painel (d) com o `c_s` acompanhando os bracos.

## Resultado (2026-09-15, t=50, 1970 s) — renderizador final, janela t in [35, 50]

| | P2R9 | **P2R10** |
|---|---:|---:|
| `c_s` do filler por faixa 0-0.2 / 0.2-0.4 / 0.4-0.6 / 0.6-0.8 / 0.8-1.0 R99 | congelado (0.29 / 0.18 / 0.14 / 0.12 / 0.10) | **0.47 / 0.41 / 0.22 / 0.12 / 0.09** |
| `c_s` do agar em volta 0.2-0.4 / 0.4-0.6 / 0.6-0.8 / 0.8-1.0 R99 | 0.045 / 0.049 / 0.059 / 0.049 | 0.291 / 0.185 / 0.112 / 0.069 |
| frac(`c_s` > 0.45) nas vivas / `contrast_cs` t=50 / `mean_cs` | 0.94 / 18.4 / 0.027 | **0.66** / 15.1 / 0.033 |
| R99 / Rmax / Rmin (t=50) | 4.74 / 4.99 / 1.17 | **4.31 / 4.54 / 1.64** |
| bracos 0.65 / 0.8 R99 / largura 0.8 R99 | 26.5 / 20.0 / 5.4 | 24.8 / **16.0** / 5.7 |
| pedacos / tortuosidade / amplitude / vivas | 0 / 1.01 / 0.154 / 622 | 0 / 1.03 / 0.170 / 668 |
| `a_mar_bio_p95` (mediana t>30) / dR/dt / iter/t / massa | 2.04 / 0.096 / 62.1 / 204.2 | 1.71 / 0.089 / 80.8 / 205.6 |

1. Objetivo (painel): **APROVADO** — `c_s` do filler decrescente do nucleo para fora, 0.41 em
   0.2-0.4 R99; painel (d) com `--cs-filler` continuo, parecido com Trinschek (b)
   (`plots/fig_P2R10_t50.png`). A predicao do perfil acertou a forma (decrescente, sem pico fora do
   corpo) e ERROU o nivel para cima: base 0.41 contra 0.15-0.30 previsto, bracos 0.12-0.22 contra
   0.05-0.15 — subestimei o fluxo que atravessa a base; e foi esse excesso que reabriu o laco abaixo.
2. Forma: **REPROVADO** — bracos em 0.8 R99 16.0 (limite 17); R99 4.31 passa no limite (4.27).
3. Quimica: **REPROVADO** em frac(`c_s` > 0.45) nas vivas 0.66 (limite 0.8); `contrast_cs` passa.
4. Numerica: aprovado (iter/t 1.43x o P2, limite 1.5x).
5. Visual: painel (d) aprovado; painel (a) **REPROVADO** — nucleo virou disco de raio ~1.7
   (Rmin 1.17 -> 1.64), bracos mais curtos (`plots/cmp_colonia_P2R10_t50.png`).

**Causa (medida): o laco da licao #90.** A conducao leva o `c_s` do nucleo ate o agar da base:
agar+limbo em 0.15-0.5 R99 com `c_s` acima do gate da colonizacao (`COL_CS_MIN`=0.3) em t≈36 —
P2R9 5%, P2R10 67%. A colonizacao (filler doador) recruta a base: agar+limbo restante na base em
t=50 1277 -> 475; vivas em 0.25-0.5 R99 107 -> 335. Mais o `c_s` mais alto a frente das pontas
(agar 0.07-0.11 contra 0.05) enfraquece o gradiente na frente: `a_mar_bio_p95` −16%, R99 −9%.

Figuras: `plots/fig_P2R10_t50.png` (tese, `--cs-filler`), `plots/cmp_colonia_P2R10_t50.png`.
