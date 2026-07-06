---
name: analisar-log-csv
description: >-
  Calcula as tendencias quantitativas obrigatorias do Protocolo Sec.2.1 do
  CLAUDE.md a partir de log.csv (raiz do repo, NAO main_output/) (mean_v, n_fast, a_marangoni,
  a_flag, mean_cs, max_cs, contrast_cs, mass_total, a_pressure, mean_c_n).
  Use SOMENTE quando precisar ler/interpretar log.csv de um run — antes de
  calibrar qualquer parametro fisico, ou como insumo de dados do protocolo
  de validacao morfologica (skill validar-morfologia). NAO use para ler
  frames/imagens.
---

# Analisar log.csv (Protocolo Sec.2.1)

## Quando usar
- Antes de propor qualquer mudanca em parametro fisico (`beta`, `sigma`,
  `gamma`, `D`, `lambda_`, `r_growth`, `c0`, `f0`, `k_consume`, etc — CLAUDE.md
  Sec.10 proibe mudanca cega sem isso).
- Como parte do Protocolo Sec.11 (validacao morfologica), etapa 2.
- Quando o usuario pedir "como esta o log", "o motor esta vivo", "checar
  colapso".

## Quando NAO usar
- Analise de frames/imagens — isso e o agent `analista-morfologia` ou a
  skill `validar-morfologia`.
- Mudanca de fator de producao/sumidouro de cs especificamente — depois de
  rodar esta skill, se a mudanca proposta toca `c_n_factor`, `growth_headroom`,
  `sigma`, `tip_boost`, `motile_boost`, `qs`, `k_consume` ou `lambda`, va
  OBRIGATORIAMENTE para a skill `protocolo-cs-zonas` tambem (Sec.3.3.6).

## Procedimento

1. **Ache o log ativo.** O `main.py` grava em `LOG_FILE = "log.csv"`, um
   caminho relativo ao diretorio de trabalho de onde `python main.py` foi
   rodado — na pratica, a **raiz do repo** (`log.csv`), **NAO**
   `main_output/log.csv` (esse arquivo nao existe; `main_output/` so tem os
   `.hdf5`, `main.log`, `profile_info.csv` e a pasta `movie/`). Confirme com
   `stat -f "%Sm %N" log.csv main_output/main.log` se as duas datas batem —
   se nao baterem, o log.csv pode ser de um run diferente do output HDF5.

2. **Leia com Python, nunca awk cru.** Ha um bug de locale documentado neste
   Mac (memoria `reference_awk_locale_bug`): `awk` sem `LC_ALL=C` pode
   truncar decimais com ponto. Use:
   ```python
   import csv
   with open('log.csv (raiz do repo, NAO main_output/)') as f:
       rows = list(csv.DictReader(f))
   ```
   Isso tambem satisfaz o hook `mark_log_read.py` (marca a sessao como tendo
   lido log.csv fresco — necessario para o guard `guard_param_change.py`
   liberar a proxima edicao de parametro). **Prefira o tool `Read` do
   Claude Code no arquivo `log.csv (raiz do repo, NAO main_output/)` diretamente quando possivel**
   — e o que o hook detecta e marca; um `python -c "..."` via Bash tambem
   funciona mas nao aciona a marca (o hook so escuta o tool `Read`).

3. **Calcule tendencia (nao so valor final) para cada metrica:**
   `mean_v`, `n_fast`, `a_marangoni`, `a_flag`, `mean_cs`, `max_cs`,
   `contrast_cs` (=max_cs/mean_cs), `mass_total`, `a_pressure`, `mean_c_n`
   (se presente). Compare inicio vs. fim vs. pico, identifique se e
   monotonico, plateau, ou oscilante.

4. **Aplique os alarmes conhecidos (Sec.2.1, Sec.9 licoes):**
   - `contrast_cs` em queda monotonica = motor morrendo, mesmo se `mean_v`
     parece saudavel (licao #8).
   - `mean_cs` crescendo sem limite alem de ~0.5 = afogamento quimico
     iminente (Sec.2.4-B).
   - `n_fast=1` isolado ou spike unico de `max_v` = particula rogue, nao
     expansao real.
   - `mean_v`/`n_fast` baixos MAS `a_flag`/`a_marangoni` > 5 e frames
     saudaveis = regime tip-only saudavel, NAO motor morto (licao #21).
   - `mass_total` crescendo > ~2x no run = suspeita de tip pumping
     (Bloqueio E, historico Pass M-B).
   - `a_pressure` sustentado fora do orcamento Sec.8 (<3-10 conforme fase)
     = over-pack ou brittle neck (licao #15/#23/#27).

5. **Compare com o orcamento de aceleracoes (Sec.8):** `a_marangoni`~6,
   `a_flag`~2, `a_drag`~6, `a_pressao`<3, `a_viscosa`~1.5, total |a| 5-15.

6. **Reporte em tabela** metrica | inicio | fim | pico | tendencia | leitura,
   seguida da fase do motor (bootstrap/ativo/declinante/morto/tip-only) e
   quais termos estao fora do orcamento.

## Nunca faca
- Nunca conclua "motor morto" so por `mean_v`/`n_fast` baixos sem checar
  `a_flag`/`a_marangoni` e frames (licao #21).
- Nunca proponha mudanca de parametro sem este passo ter sido feito antes
  na mesma sessao — e o pre-requisito do hook `guard_param_change.py`.
