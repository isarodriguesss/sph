---
name: formatar-pass-historico
description: >-
  Registra o resultado de um Pass de calibracao no CLAUDE.md Sec.9
  (Historico de Passes), seguindo exatamente o template usado em todo o
  historico existente (Resultado com data, Predicoes ex-ante vs. resultado
  real, Diagnostico/causa raiz, Licao numerada se aplicavel). Use DEPOIS que
  a skill/agent validar-morfologia (Protocolo Sec.11) e a analise de log.csv
  ja concluiram se o Pass foi sucesso/falha/parcial — nunca antes.
---

# Registrar um Pass no historico (Sec.9)

## Pre-requisito

Isto so roda DEPOIS de:
1. `analisar-log-csv` (ou `analista-log`) ter produzido as tendencias.
2. `validar-morfologia` (ou `analista-morfologia`) ter produzido o
   diagnostico das 3 etapas (Sec.11).
3. A predicao ex-ante (Sec.2.3) ja ter sido formulada ANTES da mudanca —
   se nao foi, reporte isso como um desvio de protocolo em vez de inventar
   uma predicao retroativa.

Se qualquer um desses faltar, pare e rode o que falta primeiro — nao
registre um Pass sem os dados que o sustentam.

## Template exato (copiar a estrutura, nao inventar formato novo)

```markdown
- **Pass <ID> — <resumo curto do lever> (<data ISO>, <VEREDITO EM CAIXA ALTA>):**
  <alavanca(s) mudada(s), com valores antes->depois e arquivo:linha>.

  **Predicoes ex-ante:** <o que foi previsto antes de rodar, com numeros>.

  **Resultado (t=0->Ts):** <tendencias de log.csv com numeros: mean_cs,
  contrast_cs, a_marangoni, a_pressure, mass_total, etc — inicio->fim/pico>.

  **Frames (Protocolo Sec.11):** <descricao dos frames-chave com t aproximado,
  comparado a reference.jpg / reference_result.png painel (b)>.

  **Diagnostico:** <causa raiz do resultado, citando licoes numeradas
  existentes se o mecanismo ja e conhecido (Sec.9 "Invariantes morfologicos
  descobertos"), ou propondo uma nova licao numerada se for descoberta
  inedita>.

  **Proximos passos / caminho de continuidade:** <o que fazer a seguir, uma
  alavanca por vez (Sec.10 proibicao de mudar >1 parametro de surfactante
  simultaneamente)>.
```

`VEREDITO EM CAIXA ALTA` segue o vocabulario ja usado no arquivo: `SUCESSO`,
`FALHA`, `SUCESSO PARCIAL`, `EM ANDAMENTO`, `DIAGNOSTICO ESTRUTURAL`,
`FALHA CATASTROFICA`, etc. — reutilize os termos existentes, nao invente
vocabulario novo.

## Onde inserir

- **Historico de Passes:** dentro de "### Em andamento" (Sec.9) ou mova para
  "### Resolvidos" se o Pass fechou definitivamente um bloqueio. Mantenha
  ordem cronologica.
- **Se descobriu um invariante/licao reutilizavel:** adicione a
  "Invariantes morfologicos descobertos" com o PROXIMO numero sequencial
  (a lista atual vai ate a licao mais alta ja registrada — confira o ultimo
  numero usado antes de numerar a nova).
- **Se mudou um parametro da tabela calibrada:** atualize tambem a Sec.7
  (tabela `Parametros Calibrados`) com o novo valor e uma nota curta
  referenciando o Pass.
- **Se o Pass desbloqueou/bloqueou um Bloqueio dos Criterios Duplos
  (Sec.2.4):** atualize o status (✅/❌/🔄) na secao correspondente.

## Nunca faca

- Nunca escreva o registro antes de ter os dados de log.csv E frames.
- Nunca reescreva ou apague entradas de Passes anteriores — o historico e
  append-only (e a fonte de verdade sobre o que ja foi tentado e falhou).
- Nunca invente uma licao numerada que duplique uma ja existente — confira
  a lista de licoes existentes primeiro (o `explorador` pode ajudar a
  verificar se um mecanismo parecido ja foi documentado).
