---
description: Grava o resultado de um Pass ja validado (via /validar-pass) no historico do CLAUDE.md Sec.9, no formato padrao do arquivo.
argument-hint: "(vazio — usa o resultado de /validar-pass desta sessao)"
---

# /registrar-pass — grava no historico (Sec.9)

Sem argumentos — usa o diagnostico produzido por `/validar-pass` nesta
mesma sessao.

## Passos

1. **Confirme que `/validar-pass` ja rodou nesta sessao** (ou que o usuario
   colou o resultado de uma validacao anterior). Se nao, **pare e peca**
   para rodar `/validar-pass` primeiro — nao invente o conteudo do registro.

2. **Use a skill `formatar-pass-historico`** para montar a entrada no template exato
   (Resultado com data, Predicoes ex-ante vs. real, Diagnostico, Proximos
   passos) e decidir onde inserir (Sec.9 "Em andamento" ou "Resolvidos").

3. **Edite `CLAUDE.md` diretamente** adicionando a entrada. Se um novo
   invariante/licao foi descoberto, adicione com o proximo numero
   sequencial na lista de "Invariantes morfologicos descobertos" (confira o
   ultimo numero usado antes).

4. **Se algum parametro da tabela Sec.7 mudou**, atualize a tabela tambem
   (valor novo + nota curta referenciando o Pass).

5. **Se algum Bloqueio dos Criterios Duplos (Sec.2.4) mudou de status**,
   atualize o status ali tambem.

6. **Mostre o diff** da edicao ao usuario antes/depois de aplicar (o Edit
   tool ja faz isso naturalmente) — este e um arquivo de governanca do
   projeto, precisao no registro importa mais que velocidade.

## Nunca faca
- Nunca reescreva ou apague entradas de Passes anteriores.
- Nunca registre um Pass que nao passou por `/validar-pass` (ou equivalente
  manual com log.csv + frames).
