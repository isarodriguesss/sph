---
description: Fecha o ciclo de um Pass — roda o Protocolo Sec.11 (frames + log.csv) e compara predicao ex-ante vs. resultado real (Sec.2.3), emitindo veredito manter/reverter.
argument-hint: "(opcional) lembrete da predicao feita em /novo-pass, se nao estiver na mesma sessao"
---

# /validar-pass — fecha o ciclo do Pass (Protocolo Sec.11 + Sec.2.3)

Entrada do usuario (opcional): **$ARGUMENTS**

## Passos

1. **Recupere a predicao ex-ante.** Se foi feita nesta mesma sessao (via
   `/novo-pass`), use-a do contexto. Se `$ARGUMENTS` fornece a predicao
   (sessao nova), use isso. Se nao ha predicao disponivel de nenhuma forma,
   **pare e avise**: sem predicao registrada, a comparacao Sec.2.3 fica
   incompleta — registre isso explicitamente no resultado em vez de inventar
   uma predicao retroativa.

2. **Rode o protocolo completo de validacao morfologica.** Invoque o agent
   `analista-morfologia` (ou a skill `validar-morfologia`) para as 3 etapas
   obrigatorias da Sec.11: descricao dos frames, cruzamento com log.csv,
   sintese. **Nunca aceite uma comparacao so-visual como resposta final** —
   se o agent/skill nao mostrar numeros de log.csv, peca para completar.

3. **Compare predicao vs. resultado (Sec.2.3).** Se a discrepancia for
   > 2x, isso exige investigacao antes de prosseguir para o proximo Pass —
   nao ignore silenciosamente.

4. **Verifique os Criterios Duplos (Sec.2.4)** se a mudanca tocou o motor de
   surfactante: Bloqueio A (mecanico, `a_pressure`), Bloqueio B (quimico,
   `mean_cs`/`contrast_cs` plateau), Bloqueio C (geometrico Mullins-Sekerka).
   Reporte quais criterios #1-8 passam/falham com os numeros.

5. **Emita o veredito:** manter a mudanca / reverter / manter com ressalvas
   (e quais). Se reverter, diga exatamente o que desfazer (arquivo:linha).

6. **Aponte para `/registrar-pass`** para gravar o resultado no historico
   (CLAUDE.md Sec.9) — so depois deste passo ter concluido.

## Nunca faca
- Nunca declare "consistente com a referencia" sem ter cruzado com
  log.csv (violacao documentada do Pass K.13 — revertido por isso).
- Nunca pule direto para "sucesso" so porque um frame isolado parece bom.
