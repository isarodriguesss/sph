---
name: validar-morfologia
description: >-
  Executa o Protocolo Mandatorio de Analise Morfologica Comparativa (CLAUDE.md
  Sec.11) em 3 etapas: descricao quantitativa dos frames, validacao cruzada
  com log.csv, sintese. Use SEMPRE antes de declarar um Pass "bem-sucedido"
  ou "consistente com a referencia" — uma comparacao visual isolada e
  considerada uma FALHA ANALITICA neste projeto (ja aconteceu no Pass K.13 e
  foi revertido). Tambem use quando o usuario pedir para comparar frames
  com reference.jpg ou reference_result.png.
---

# Validar morfologia (Protocolo Sec.11 — 3 etapas obrigatorias, nesta ordem)

## O incidente que esta skill existe para prevenir

Pass K.13 (2026-04-20): o frame 110 foi declarado "morfologia indistinguivel
de reference.jpg" baseado em inspecao visual isolada. Analise temporal
subsequente revelou `contrast_cs` em colapso monotonico (131->32) e `mean_cs`
em crescimento sem limite (0.26->1.17) — afogamento global do motor que a
foto sozinha escondia. **O Pass foi revertido.** Nunca repita isso.

## Etapa 1 — Analise morfologica quantitativa/qualitativa

Leia os frames mais recentes de `main_output/movie/` (tool `Read` em
arquivos de imagem). Descreva, sem adjetivos vagos:
- Razao comprimento/largura dos dendritos.
- Morfologia da ponta: afiada / arredondada-bulbosa / bifurcando.
- Estrutura de ramificacao secundaria/terciaria (presente ou ausente).
- Densidade/espacamento dos bracos; baias limpas ou preenchidas por halo.

Compare com:
- **reference.jpg** (PA14 experimental): ~15-20 dendritos radiais, AR>=1:5,
  nucleo coeso, picos de surfactante nas pontas.
- **reference_result.png painel (b) Fingering** (Trinschek 2018 — meta
  morfologica primaria, Sec.2.2): 7-9 dedos finos, baias estaticas, halo de
  cs alem da biomassa.

## Etapa 2 — Validacao fisica cruzada com log.csv (NAO PULE)

Invoque a skill `analisar-log-csv` (ou leia `log.csv (raiz do repo, NAO main_output/)`
diretamente) e cruze CADA observacao visual com um numero:
- Crescimento/estagnacao visual <-> tendencia real de `mean_v`.
- "Motor saudavel" <-> `contrast_cs` alta e estavel (nunca so `mean_v`
  crescente — pode ser ativacao bulk espuria, licao #8).
- Pontas bulbosas/grossas <-> checar se `a_pressure` domina sobre
  `a_marangoni` fraco.
- `mean_v` baixo + frames com atividade concentrada nas pontas <-> pode ser
  regime tip-only saudavel (licao #21) se `a_flag`/`a_marangoni` > 5 —
  **nao declare "motor morto" sem essa checagem.**

## Etapa 3 — Sintese e diagnostico final

So depois de 1 e 2, emita o veredito. Formato exigido (Sec.11):

**PROIBIDO** (superficial):
> "Sim, parece consistente, ambos mostram multiplos bracos."

**EXIGIDO** (ancorado em dados):
> "Nao, a simulacao falha em replicar a referencia. Os bracos sao curtos e
> bulbosos ao inves de finos e ramificados. O log.csv confirma:
> `contrast_cs` colapsa de X para Y, enfraquecendo `a_marangoni`; isso reduz
> `mean_v` e permite que `a_pressure` domine, gerando a morfologia arredondada
> observada."

Se a conclusao for FALHA: identifique o mecanismo mais provavel usando o
historico de failure modes da Sec.9 (brittle neck, saturacao quimica Sec.2.4-B,
core pinning Sec.2.4-A, tip pumping, halo radial) — nao invente um mecanismo
novo sem descartar os ja catalogados primeiro.

Se a conclusao for SUCESSO: liste explicitamente quais dos criterios dos
Criterios Duplos (Sec.2.4, itens #1-8) estao satisfeitos, com os numeros do
log ao lado de cada um.

## Nunca faca
- Nunca pule a Etapa 2.
- Nunca declare "Pass L pode comecar" so por esta validacao — o criterio da
  Sec.2.2 exige sustentacao ao longo do tempo (AR>=1:5, baias estacionarias,
  tip-splitting visivel), nao um frame isolado bonito.
