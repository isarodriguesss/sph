# P2R17 — W=5 com o filler CONDUZINDO: devolver o campo que a largura apagou

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R16** (W=5 sobre o P2R15). Alavanca unica:
`FILLER_CS_CONDUZ` 0 -> 1 (o mecanismo do P2R10; `FILLER_CS_D` = 0, entao os pares com filler usam
`D_ext` = 0.08, e o sumidouro do filler fica em 0.5 lambda, o padrao). t=50, `SEED` fixo, 10 threads.

## Diagnostico que motiva (P2R16)
Alargar o rastro (W 3 -> 5) melhorou a forma — AR 18.4 -> **12.0**, largura 5.6 -> 8.7 dx em 0.6 R99
e 9.6 dx em 0.75 R99, bracos em 0.75 R99 20 -> **12** por fusao — e **cobrou o motor**:
`a_mar_bio_p95` (mediana t>30) 1.67 -> **0.44** e `c_s` do corpo 0.36 -> **0.19**.

Causa: cada particula de agar que o rastro converte deixa de carregar `c_s` (congelada e
transparente, licao #39). Pintar mais largo **apaga o campo em volta do braco**, e e o gradiente ali
que move a frente. Ou seja `RASTRO_W` nao e alavanca puramente geometrica. Com o filler conduzindo,
o material convertido volta a participar do campo.

## Predicao
- `c_s` do corpo: 0.19 -> **0.45-0.60** do teto (no passo P2R9 -> P2R10 a conducao deu +0.31, e aqui
  ha muito mais filler para preencher); baias 0.15 -> 0.20-0.28.
- **Motor: 0.44 -> 0.9-1.5** — e a hipotese central: a conducao devolve o gradiente que a largura
  apagou. Se ficar abaixo de 0.9, a perda do P2R16 nao era do campo apagado.
- Forma quase intacta (a geometria do rastro nao muda): ciclo 0.49 -> 0.50-0.58, AR 12.0 -> 10-12,
  bracos em 0.75 R99 12-16.
- Nucleo 0.29 -> **0.40-0.48** (a conducao leva `c_s` ao agar da base e abre o gate, licao #94-H);
  R99 4.68 -> 4.5-5.0.
- `max_cs` cravado em 0.494: `D` nao muda, entao nao ha risco da instabilidade difusiva da licao #96.

## Criterios (t=50 e janela t in [35, 50])
1. **Objetivo:** `a_mar_bio_p95` (mediana t>30) >= 1.0 **e** `c_s` do corpo >= 0.40 do teto.
2. **Forma preservada:** ciclo em 0.6 R99 >= 0.45; AR <= 13; largura >= 0.09 Rmax.
3. **Guarda:** nucleo <= 0.45 Rmax; R99 >= 4.4; pedacos soltos <= 0.5.
4. Numerica: `max_cs` <= 0.5 (aborto precoce se passar de 0.6 na iteracao 200); iter/t <= 1.6x o P2;
   massa sem runaway.
5. Visual (§11): classe (b); painel (d) com o campo acompanhando os bracos largos.

## Resultado (2026-09-16, t=50, iter/t 89.7) — HIPOTESE REFUTADA, e a forma vai ao alvo por outro mecanismo

| janela t in [35,50] / t=50 | P2R15 (W=3) | P2R16 (W=5) | **P2R17 (W=5 + conduz)** | alvo |
|---|---:|---:|---:|---:|
| `c_s` do corpo / teto | 0.36 | 0.19 | **0.49** | limite 0.40 ✅ |
| `c_s` das baias / teto | 0.25 | 0.15 | **0.25** | — |
| **`a_mar_bio_p95` (mediana t>30)** | 1.60 | 0.44 | **0.60** | **limite 1.0 ❌** |
| ciclo em 0.6 Rmax (**EDT, corrigida**) | 0.46 | 0.60 | **0.48** | limite 0.45 ✅ |
| **largura 0.6R / 0.75R (EDT)** | 2.8 / 3.0 dx | 5.0 / 5.0 dx | **5.0 / 5.0 dx** | 0.064 R, limite 0.09 ❌ |
| largura do miolo (0.3R) | 19.0 dx | 18.8 dx | **22.5 dx** | — |
| AR | 18.4 | 12.0 | **10.4** | limite 13 ✅ (faixa 6.5-11.2 ✅) |
| nucleo / Rmax | 0.27 | 0.29 | **0.39** | limite 0.45 ✅ (alvo 0.36) |
| baia (Rmin/Rmax) | 0.26 | 0.29 | **0.38** | alvo 0.34 |
| dedos / area-disco | 21 / 0.38 | 18 / 0.41 | **15 / 0.38** | 9-15 / 0.35-0.39 ✅ |
| halo `L` / Rmax | 0.15 | 0.15 | **0.18** | alvo 0.18 ✅ |
| R99 / dR/dt | 4.57 / 0.089 | 4.68 / 0.095 | **4.11 / 0.075** | limite 4.4 ❌ |
| amplitude / pedacos soltos | 0.163 / 0.0 | 0.145 / 0.0 | **0.231 / 0.0** | limite 0.5 ✅ |
| `max_cs` / massa | 0.494 / 205 | 0.494 / 204 | **0.494 / 205** | ✅ |

**1. Objetivo: REPROVADO.** O campo voltou exatamente como previsto (corpo 0.19 -> 0.49, previsao
0.45-0.60) e o motor NAO acompanhou: 0.44 -> 0.60 contra o limite 1.0 e os 1.60 do P2R15. O proprio
criterio dizia "se ficar abaixo de 0.9, a perda do P2R16 nao era do campo apagado" — **e nao era.**
A perda de motor do P2R16 vem da LARGURA em si (mais agar convertido = menos vivas: 783 -> 696 -> 659,
e o gradiente distribuido sobre um braco 2x mais largo), nao do `c_s` que a conversao congelava.

**2. Forma: a LARGURA DO BRACO NAO MUDOU (5.0 dx nos dois) — reprova so contra o alvo absoluto.**
A primeira leitura ("a conducao desfez o engrossamento, 4.4 -> 2.6 dx") era **artefato da regua** —
ver a correcao no fim. Medido por EDT, P2R16 e P2R17 tem bracos identicos. O que a conducao fez foi
engrossar o MIOLO (18.8 -> 22.5 dx em 0.3R; nucleo solido 0.29 -> 0.39 Rmax) e **encolher o raio 12%
(4.78 -> 4.22)** — o laco do disco (licoes #90/#94-H): o `c_s` conduzido chega ao agar da base, abre
o gate da colonizacao e o nucleo come a base. O ciclo em 0.6R cai (0.60 -> 0.48) so porque `R`
encolheu; em 0.75R e igual (0.31 vs 0.32). Ciclo passa (0.48 >= 0.45); largura relativa reprova
(0.064 < 0.09) — **mas o P2R16 tambem reprova nesse eixo (0.056)**: nenhum dos dois atingiu.

**3-4. Guarda e numerica:** nucleo 0.39 <= 0.45 ✅, pedacos 0.0 ✅, `max_cs` cravado no teto ✅,
massa sem runaway ✅; **R99 4.11 < 4.4 ❌**; iter/t 89.7 contra o limite 90 — passou no fio.

**5. Visual:** classe (b) — nucleo redondo grande, ~15 bracos radiais, baias abertas ate a base,
halo continuo e amplo no painel (d). `plots/fig_P2R17_w5_conduz.png`.

**O ACHADO, apesar da reprovacao:** contra as reguas da literatura (licoes #95/#99) este e o melhor
conjunto morfologico ja medido no projeto — **nucleo 0.39 (alvo 0.36), baia 0.38 (0.34), dedos 15
(9-15), area/disco 0.38 (0.35-0.39), AR 10.4 (6.5-11.2) e halo L/Rmax 0.18 (0.18, exato)** — cinco
eixos dentro da faixa ou a menos de 12% dela. **Mas chegou la mexendo no DENOMINADOR:** o nucleo
cresceu comendo a base e o raio encolheu 12%, com a largura do braco parada em 5 dx. Acertar o numero
da literatura por um mecanismo que a licao #90 ja identificou como defeito nao e o mesmo que
reproduzir a morfologia.

**Regra:** `nucleo/Rmax`, `AR`, `area/disco` e `largura/Rmax` sao RAZOES — encolher `Rmax` produz os
mesmos numeros que acertar o numerador. Cruzar sempre com medidas absolutas (largura em dx, Rmax).

## Correcao de regua (2026-09-16, apontada pela usuaria ao comparar com o viewer)

A primeira leitura desta rodada dizia "ciclo 0.56 -> 0.39, largura 4.4 -> 2.6 dx". **Era artefato.**
A regua v1 amostrava particulas numa fatia `|r - R| < 0.6 dx` e contava arcos ocupados em 720 bins:
em 0.6 Rmax isso pega **172 particulas para 720 bins** (ocupacao crua 0.17) e acha **93 "bracos"**
onde a imagem tem ~14 — fragmentacao de amostragem. Varrendo a tolerancia (0.6 / 1.0 / 1.5 / 2.5 dx)
a contagem vai de 93 a 73 e a ocupacao de 0.17 a 0.43: **o numero era funcao da regua.**

Corrigido para EDT sobre a mascara rasterizada (regua da licao #75-A), em
[tools/ciclo_largura.py](../../tools/ciclo_largura.py). Dois erros de metodo a nao repetir:
(1) o teste barato que expoe isso e **contar o N da amostra** — sempre reportar junto da metrica;
(2) eu persisti a regua ERRADA em `tools/` como "correcao do erro de metodo", e estar num arquivo a
fez parecer confiavel. **Persistir nao valida — validar e conferir contra as particulas ou o viewer
(licao #92) ANTES de ranquear.**
