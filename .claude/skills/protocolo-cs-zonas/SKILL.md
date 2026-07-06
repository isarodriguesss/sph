---
name: protocolo-cs-zonas
description: >-
  Executa o protocolo OBRIGATORIO da Sec.3.3.6 do CLAUDE.md antes de mudar
  qualquer fator de producao ou sumidouro de cs (c_n_factor, growth_headroom,
  sigma, tip_boost, motile_boost, qs, k_consume, lambda). Calcula cs_infinito
  em 4 zonas representativas para achar onde o pico de cs cairia e se a
  mudanca produz push outward (correto) ou inward (patologico). Use SEMPRE
  antes de propor essas mudancas — e uma proibicao explicita do CLAUDE.md,
  nao uma sugestao.
---

# Protocolo cs em 4 zonas (Sec.3.3.6 — OBRIGATORIO, NAO PULE)

## Por que isso existe

CLAUDE.md Sec.3.3 documenta uma falacia recorrente: "cs concentrado nas
pontas puxa a colonia para fora". **Falso.** As forcas `MarangoniForce` e
`FlagellarForce` apontam na direcao `-∇cs` (de alto cs para baixo cs). Um
pico de cs no rim externo empurra para fora APENAS a biomassa que vem DEPOIS
do pico (radialmente) — se depois do pico e agar sem biomassa, a forca nao
faz nada util. TODA a biomassa ANTES do pico recebe push INWARD (compactacao).
Mudar um fator de producao/sumidouro sem calcular onde o pico cai e "mudar as
costas viradas para a fisica" — proibido explicitamente na Sec.10.

## Quando usar

Antes de propor mudanca em: `c_n_factor`, `growth_headroom`, `sigma`,
`tip_boost`, `motile_boost`, `qs`, `k_consume`, `lambda` (`lambda_ext`,
`lambda_ext_ratio`). Sem excecao — nem para mudancas "biologicamente
motivadas" (Sec.3.3.6: "a coerencia biologica nao implica push outward
correto").

## Procedimento

1. **Calcule `cs_∞` nas 4 zonas representativas** (formula Sec.3.3.6):
   ```
   cs_∞ = (sigma * qs * growth_headroom * c_n_factor) / (lambda_eff + k_consume * rho_b)
   ```
   - Interior profundo: `rho_b=1, c_n=0`
   - Mid-arm: `rho_b=0.5, c_n=0.4`
   - Outer rim: `rho_b=0.2, c_n=0.8`
   - Agar: `rho_b=0, c_n=1`

   Use os valores ATUAIS dos parametros (le `main.py`/`src/equations.py`
   para os valores calibrados correntes — Sec.7 tem a tabela mas o codigo e
   a fonte de verdade) e os valores PROPOSTOS, lado a lado.

2. **Identifique a zona com maior `cs_∞`** — e onde o pico de cs cairia no
   perfil radial.

3. **Verifique se ha biomassa contigua (`rho_b > 0.1`)** entre essa zona e a
   fronteira com o agar. Sem isso, o push outward nao tem em quem atuar.

4. **Classifique pelo mapa da Sec.3.3.4:**
   | Pico em | Push outward em | Morfologia |
   |---|---|---|
   | `rho_b≈1.0` (interior) | ~99% da colonia | Isotropica uniforme, sem seletividade |
   | `rho_b≈0.5` (mid-arm) | rim externo | Dendritica seletiva SE localizavel azimutalmente |
   | `rho_b≈0.2-0.3` (outer rim) | so `rho_b<0.2` ou agar | Colonia comprime |
   | `rho_b≈0` (agar virgem) | nenhuma | Patologica — colapso |

5. **Documente a predicao**: tabela `cs_∞` por zona ANTES e DEPOIS da
   mudanca, e para onde o pico se move. Isso e o artefato que justifica a
   mudanca — sem ele, a mudanca esta em violacao da Sec.10.

## Veredito

- Se depois da mudanca o pico cai no interior profundo ou no agar isolado:
  **rejeite** a mudanca (patologica).
- Se cai em mid-arm/swarmer zone com localizacao azimutal possivel (via
  motile_boost ou equivalente): **caminho valido** para seletividade
  dendritica — prossiga.
- Se cai no interior mas satura uniformemente: pode ajudar o motor mas nao
  a morfologia — sinalize essa limitacao explicitamente ao propor.

## Nunca faca
- Nunca proponha a mudanca so com justificativa biologica sem a tabela
  `cs_∞` de 4 zonas.
- Nunca assuma que "mais surfactante nas pontas" = "mais push outward" sem
  verificar a posicao do pico primeiro.
