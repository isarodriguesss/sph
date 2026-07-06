---
description: Inicia um novo Pass de calibracao seguindo o Protocolo Sec.10/Sec.2 do CLAUDE.md — le log.csv, exige predicao quantitativa, roda o protocolo de zonas de cs se aplicavel, e implementa UMA alavanca.
argument-hint: <hipotese ou parametro a mudar, ex: "reduzir D_ext para conter winner-steal">
---

# /novo-pass — inicia um Pass de calibracao (Protocolo Sec.10)

Entrada do usuario: **$ARGUMENTS**

Siga o protocolo do CLAUDE.md Sec.10 ("Ao receber uma tarefa de
calibracao") do inicio ao fim. Nao pule etapas.

## Passos

1. **Leia log.csv primeiro.** Use a skill `analisar-log-csv` (ou invoque o
   agent `analista-log`) para entender a fase atual do motor antes de
   propor qualquer coisa. Isso tambem satisfaz o hook `guard_param_change.py`
   — sem essa leitura fresca nesta sessao, a edicao do parametro sera
   bloqueada.

2. **Verifique se ja foi tentado.** Invoque o agent `explorador` com a
   hipotese de `$ARGUMENTS` para checar no historico de Passes (CLAUDE.md
   Sec.9) se uma alavanca igual ou muito parecida ja foi tentada e
   documentada como falha. Se sim, **nao repita** — proponha uma variacao
   explicitamente diferente ou pare e informe o usuario.

3. **Se a mudanca toca producao/sumidouro de cs** (`c_n_factor`,
   `growth_headroom`, `sigma`, `tip_boost`, `motile_boost`, `qs`,
   `k_consume`, `lambda`): rode OBRIGATORIAMENTE a skill `protocolo-cs-zonas`
   (Sec.3.3.6) antes de prosseguir. Sem a tabela de `cs_∞` em 4 zonas, a
   mudanca esta em violacao explicita do CLAUDE.md.

4. **Se a mudanca toca refinamento adaptativo, integrador customizado, nova
   equacao SPH, trigger de kernel, ou tratamento de fronteira:** invoque o
   agent `guardiao-literatura` para verificar a exigencia de citacao de
   Liu/Violeau (Sec.10) e conferir se conflita com uma proibicao explicita
   (hard pinning, D_ext/lambda_ext simultaneos, etc).

5. **Formule a predicao quantitativa (Sec.2.3)** ANTES de editar qualquer
   arquivo: "mudar X de A para B deve aumentar/reduzir <metrica> em ~N%".
   Escreva isso na conversa explicitamente — e o que sera comparado depois
   em `/validar-pass`.

6. **Implemente UMA alavanca por vez** (Sec.2.4-D, Sec.10 proibicao de
   combinacoes simultaneas de parametros de surfactante). Edite
   `main.py`/`src/equations.py`/`src/scheme.py` conforme necessario.

7. **Pare e pergunte antes de rodar a simulacao.** Runs SPH podem ser longos
   (o hook `warn_long_run.py` vai lembrar disso ao detectar `make run`/
   `python main.py`). Nao inicie em background sem combinar.

8. **Ao terminar o run, aponte para `/validar-pass`** para fechar o ciclo
   (Protocolo Sec.11 + comparacao predicao vs. resultado).

## Nunca faca
- Nunca mude mais de um parametro fisico nesta sessao de Pass.
- Nunca declare sucesso aqui — isso e escopo de `/validar-pass`.
- Nunca proponha Pass L (rugosidade) — permanece bloqueado (Sec.2.2) ate
  criterio morfologico sustentado.
