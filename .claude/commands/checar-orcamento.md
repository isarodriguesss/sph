---
description: Checagem rapida do log.csv mais recente contra o orcamento de aceleracoes da Sec.8 do CLAUDE.md, sem abrir um ciclo completo de Pass.
argument-hint: "(vazio — usa log.csv (raiz do repo, NAO main_output/))"
---

# /checar-orcamento — leitura rapida contra o orcamento Sec.8

Sem cerimonia de Pass — so uma leitura pontual. Use quando quiser saber "como
esta a simulacao agora" sem iniciar `/novo-pass`.

## Passos

1. Invoque o agent `analista-log` (ou a skill `analisar-log-csv`) sobre o
   log ativo.
2. Reporte, para as ultimas linhas disponiveis, cada termo contra o alvo da
   Sec.8:

   | Termo | Alvo | Valor atual |
   |---|---|---|
   | `a_marangoni` | ~6 (so interface/pontas) | |
   | `a_flag` | ~2 | |
   | `a_drag` | ~6 | |
   | `a_pressure` | <3 | |
   | `a_viscosa` | ~1.5 | |
   | Total \|a\| | 5-15 | |

3. Sinalize quais termos estao fora do alvo e, se possivel, aponte a causa
   provavel usando o historico de failure modes da Sec.9 (sem propor mudanca
   de parametro aqui — isso e escopo de `/novo-pass`).

## Nunca faca
- Nunca proponha ou implemente mudanca de parametro neste comando — e
  leitura pura. Para agir, use `/novo-pass`.
