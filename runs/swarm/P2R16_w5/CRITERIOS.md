# P2R16 — `RASTRO_W` 3 -> 5: engrossar o braco ate o ciclo de trabalho da literatura

Pre-registrado em 2026-09-16, antes de lancar. Base: **P2R15**. Alavanca unica: `RASTRO_W`
3.0 -> 5.0 dx. Resto identico (inclusive `AGAR_CS_LAMBDA` = 0.125). t=50, `SEED` fixo, 10 threads.

## Por que esta alavanca agora (licao #99, fontes externas [T11]-[T14])
[T11] Deng et al. 2014 mede em PA14: braco **2-5 mm** de largura com espacamento **4-5.5 mm**, ou
seja razao largura/passo (ciclo de trabalho) **0.36-1.25, mediana ~0.64** — o nosso e 0.41, na
ponta baixa. [T12] da AR ≈ 3.8 (16.85 mm / 4.4 mm), contra os nossos ~12. E [T14] (MSI) mostra o
ramnolipideo delineando os tendrils e o agar ALEM deles, nao uniforme sobre a colonia — o que
despriorizou a rota de encher o campo interno (P2R14) e deixou a LARGURA como o desvio a atacar.

Como a contagem de dedos esta travada (licao #98: dedo = rastro de uma celula isolada), a unica
via para o AR e o ciclo de trabalho, pela identidade
`AR = (1 - nucleo) * n / (ciclo * 2 pi * 0.6)`.

## Predicao (tendencia medida no par limpo P2R4 W=2 -> P2R5 W=3: +2.4 dx de largura e +0.10 de
## ciclo por unidade de W; +0.07 de area/disco; R99 +8%)
- Ciclo em 0.6 R99: 0.41 -> **0.61**; largura 5.6 -> **10.4 dx (0.115 Rmax)**, dentro da faixa
  0.10-0.19 Rmax de [T11].
- **AR 18.4 -> 7.6-8.5** pela identidade acima (com n ~25 e nucleo ~0.30). Alvo: 3.8 ([T12]) a 11.2.
- area/disco 0.38 -> **0.45-0.52** — ACIMA da faixa 0.35-0.39 que extrai das imagens. Registro que
  essa faixa conflita com o ciclo de [T11]: uma colonia com ciclo 0.64 tem area/disco alta por
  construcao. Vou reportar as duas e julgar pelo par (ciclo, AR), conforme a licao #99.
- Nucleo 0.27 -> **0.30-0.38** (no par W=2->3 o nucleo quase nao mexeu: 0.25 -> 0.26, mas em W=5 os
  rastros vizinhos comecam a se fundir perto da base).
- R99 4.57 -> **5.0-5.3**; vivas 783 -> 600-700 (mais agar convertido deixa menos recrutavel).
- Campo: corpo 0.36 -> 0.30-0.36 e baias 0.25 -> 0.20-0.25 (mais filler congelado dentro).
- Baias abertas caem ~30%; `pedacos soltos` deve seguir 0.

## Criterios (t=50 e janela t in [35, 50])
1. **Objetivo:** ciclo em 0.6 R99 >= 0.55 **e** largura >= 0.10 Rmax **e** AR <= 9.
2. **Guarda:** nucleo <= 0.40 Rmax; R99 >= 4.4; pedacos soltos <= 0.5.
3. Forma: dedos (picos de R(theta)) >= 15 — se cairem abaixo disso por FUSAO, e ganho, mas tem de
   vir com ciclo >= 0.55, nao com bracos sumindo.
4. Campo: corpo >= 0.30 do teto; halo `L`/Rmax >= 0.13.
5. Numerica: `max_cs` <= 0.5; iter/t <= 1.6x o P2 (limite 90); massa sem runaway.
6. Visual (§11): classe (b), bracos visivelmente mais largos, baias ainda abertas ate a base.
