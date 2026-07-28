# Plano de implementação — Preenchimento por rastro (wake-deposit)

**Proposta (confirmada 2026-07-28):** quando qualquer partícula da colônia se
desloca o suficiente para caber outra no espaço que deixou, inserir uma partícula
**inerte** no seu rastro. Alvo: preencher o vácuo cinemático dos dendritos (frontier
ativa, `rho_b<0.5`) sem retornar `NaN` na consulta do cliente, sem matar o fingering e
sem reacender o runaway.

Defaults confirmados pelo usuário:
1. **Gate de déficit geométrico** — só deposita se o ponto vago estiver genuinamente
   vazio ("cabe outra partícula": vizinho mais próximo > `prox`). Evita super-inserção
   no ágar (medição: gatilho puro daria +122% de N).
2. **Filler inerte-mas-livre** — não cresce, não produz cs, **não sente o motor**
   (Marangoni/flagelar), MAS participa da hidrodinâmica normal (pressão/viscosidade/
   arrasto). **Não pinado.** Evita o arrasto-de-obstáculo-congelado na frente.
3. **Escalares = valor local** — nasce herdando `cs, rho_b, c_n` da partícula que o
   depositou (= valor local no ponto), para não diluir/afogar o campo cs que o cliente
   consulta. Sem fonte/sumidouro de cs próprios (só decaimento passivo + difusão como
   carreador).

---

## 1. Escopo e método

- **Uma alavanca** (§2.3): flag nova `use_wake`. Baseline mantido idêntico (C3.4
  estrutural + KGC no estado atual). Rodar `use_wake=True` isola o mecanismo.
- **Custo de dt zero**: fillers nascem com `h=h0` (escapa do colapso de dt #29/#38 do
  Pass N).
- **Coexiste com C3.4**: a inserção estrutural (núcleo, `rho_b>0.5`, pinada) continua;
  o wake-deposit cobre a frontier (`rho_b<0.6`, livre). São dois mecanismos distintos.

## 2. Ancoragem teórica (§10 — obrigatório)

- **Liu §3.3** (consistência de partícula): restaurar a contagem de vizinhos nos braços
  recupera a consistência do operador SPH (o `∇cs` da Marangoni deixa de ser espúrio).
- **Violeau §3.6**: `σ_a → 1` é a medida direta do erro; é o critério de sucesso.
- **Liu §6.5** (superfície livre / partículas de suporte "dummy"): o filler inerte é
  dessa classe — suporte de kernel, **não** um regularizador de interface (shifting).
- **Distinção das lições que falharam:**
  - #31 (shifting matou fingering): aquilo **movia a frente**; aqui a frente se move
    pela própria física, intocada — só adicionamos suporte no rastro já vago.
  - #34/#36 (C4 ativo reacendeu runaway): aquele filler **produzia cs e crescia**; o
    nosso é inerte (sem cs-source, sem crescimento).
  - #29/#38 (dt colapsa com h pequeno): `h=h0`, sem colapso.
  - #35 (frontier vs estrutural): gatilho por deslocamento **não dispara no núcleo
    pinado** (`v≈0`), resolvendo o over-pack do C3 no centro (#27).

## 3. Design final do mecanismo

**Rastreio Lagrangiano por partícula.** Cada partícula guarda a posição do último
depósito `(x_dep, y_dep)`. A cada `WAKE_FREQ` passos, para partículas da colônia
(`rho_b > WAKE_RHO_B_MIN`):

```
disp = |(x,y) − (x_dep,y_dep)|
se disp ≥ WAKE_DISP (≈ dx):                        # moveu o bastante p/ caber outra
    spot = (x_dep, y_dep)                            # o rastro (espaço vago)
    se vizinho_mais_próximo(spot) > WAKE_PROX:       # "cabe outra" (déficit real)
        inserir filler inerte-livre em spot:
            escalares = valores locais da partícula (cs, rho_b, c_n, noise)
            u = v = 0 (nasce em repouso; livre depois)
            is_filler = 1, is_wake = 1
    (x_dep, y_dep) ← (x, y)                           # reset SEMPRE que moveu ≥ WAKE_DISP
```

Auto-limitante por construção: núcleo pinado (`v≈0`) → `disp` nunca cresce → 0
depósitos; ágar em deriva → o vago é refluído pela rede → `vizinho<WAKE_PROX` → 0
depósitos; **braços esticando → vago real → depósito.** Cap `WAKE_MAX` por call.

### Semântica das flags (duas classes de filler)

| | cresce? | produz cs? | sente motor? | hidro (P/visc/drag)? | pinado? |
|---|:---:|:---:|:---:|:---:|:---:|
| partícula normal | sim | sim | sim | sim | se `rho_b≥0.8` ou `c_n<0.6` |
| filler **estrutural** (C3.4) | não | não | (pinado) | (pinado) | **sim** |
| filler **rastro** (novo) | não | não | **não** | **sim** | **não** |

## 4. Mudanças por arquivo

### `src/particles.py` (+ `main.py` add_property)
- Nova propriedade `is_wake` (default 0.0).
- Novas propriedades `x_dep`, `y_dep` inicializadas com `x`, `y` (posição inicial =
  primeiro "último depósito").
- Registrar as três em `add_output_arrays` se quiser inspecioná-las nos frames.

### `src/equations.py` — tornar o filler inerte ao MOTOR (novo)
- `MarangoniForce.loop` (ou `post_loop`): `if d_is_filler[d_idx] > 0.5: return`
  (não recebe força de Marangoni). Ancorar: filler não é bactéria — não responde ao
  gradiente de surfactante.
- `FlagellarForce.post_loop`: idem `if d_is_filler[d_idx] > 0.5: return` antes de
  aplicar `-f0·gate·n̂(∇cs)`.
- **Sem mudança** em `BiomassGrowth` (já gateado por `is_filler`, [:30]),
  `SurfactantEquation` (qs=0 por `is_filler`, [:275]), `ViscousForce`, `LinearDrag`,
  `BiomassEOS` (filler PARTICIPA — é fluido de suporte).

### `src/scheme.py` — `CustomEulerStep.stage1` — não pinar o filler de rastro
- Adicionar `d_is_wake` à assinatura.
- Trocar a condição de pin:
  ```python
  # antes: if rho_b>=0.8 or c_n<0.6 or is_filler>0.5:
  if (d_rho_b_grown[d_idx] >= 0.8 or d_c_n[d_idx] < 0.6
          or d_is_filler[d_idx] > 0.5) and d_is_wake[d_idx] < 0.5:
      d_u = d_v = 0.0        # normais + filler estrutural
  else:
      integra normalmente    # inclui filler de rastro (is_wake=1)
  ```
  Filler de rastro (`is_wake=1`) nunca é pinado — nem pela regra de `c_n<0.6` (arm
  depletado), que reintroduziria o arrasto-congelado.

### `main.py` — constantes + bloco wake-deposit no `post_step`
- Constantes (junto às `INSERT_*`):
  ```python
  use_wake      = True
  WAKE_FREQ     = 100      # passos entre varreduras (rastreio é contínuo; check periódico)
  WAKE_DISP     = 1.0*dx   # deslocamento p/ disparar ("cabe outra")
  WAKE_PROX     = 0.7*dx   # vago só conta se vizinho+próximo além disso
  WAKE_RHO_B_MIN= 0.05     # gate barato p/ colônia (não rastreia ágar profundo)
  WAKE_MAX      = 150      # cap de depósitos por call
  ```
- Bloco novo em `post_step` (após o C3, análogo em estrutura):
  1. `disp = hypot(x−x_dep, y−y_dep)`.
  2. `cand = where((rho_b>WAKE_RHO_B_MIN) & (disp>=WAKE_DISP))`.
  3. `tree = cKDTree(positions)`; para cada `cand` (até `WAKE_MAX`): consultar
     `tree.query([x_dep,y_dep])`; se `> WAKE_PROX` e não colide com já-adicionados →
     acumular filler em `spot` com escalares locais da própria `cand`.
  4. `fluid.x_dep[cand] = fluid.x[cand]; fluid.y_dep[cand] = fluid.y[cand]` (reset).
  5. `add_particles(..., is_filler=1, is_wake=1, x_dep=spot_x, y_dep=spot_y, u=0, v=0)`;
     `append_parray`; `solver.nnps.update()`.
  6. Log: reusar/adicionar contador `wake_spawned` na linha do CSV.

## 5. Fases de implementação

- **Fase 0 — propriedades e flags** (particles.py, main.py add_property/output). Rodar
  1 passo p/ garantir que os arrays existem e o run não quebra (`use_wake=False`).
- **Fase 1 — inércia ao motor + integrador** (equations.py gates, scheme.py is_wake).
  Com `use_wake=False`, resultado deve ser **idêntico ao baseline** (fillers estruturais
  seguem pinados; nada muda). Regressão-guard.
- **Fase 2 — bloco wake-deposit** (main.py post_step) com `use_wake=True`.
- **Fase 3 — validação t=50s** (Protocolo §11: frames + log). Se limpo → **t=100s**.

## 6. Critérios de aceitação

Rodar `use_wake=True`, resto no baseline. Comparar vs **T2g / C3.4**:

| Critério | Alvo |
|---|---|
| σ_a dos braços (`mean_sig_arms`, frontier `rho_b∈[0.1,0.5]`) | **sobe p/ >0.85** (hoje ~0.95 agregado, mas com miolo em ~0.03) |
| Vácuo do miolo (scatter cru / `frac_lowsig_arms`) | **cai** — braços deixam de ter canal branco |
| Fingering vivo (`mean_v`, `contrast_cs`, morfologia) | **≈ T2g** (`contrast_cs` platô ~12–14; dendritos preservados) |
| Massa (`mass_total`) | **< ~110** em t=100s (bounded; sem runaway estilo #36) |
| dt / estabilidade | `dt ~0.02`; `a_pressure ≤ 4`; `max_v < 1` |

**Falha morfológica** (colônia congela / vira disco / dendritos borram) ⇒ o filler
livre ainda acopla demais à frente ⇒ **abortar e ir para A3** (resolução global fina).

## 7. Riscos e fallback

- **Arrasto na frente** (fingering degrada mesmo com filler livre + sem motor): acoplamento
  de pressão/viscosidade do filler à frente. Fallback: **A3** (dx 0.054→0.036), robusto.
- **Não densifica** (vácuo cinemático reabre entre depósitos): reduzir `WAKE_FREQ` /
  `WAKE_DISP`. Se estrutural → A3.
- **Massa cresce demais**: apertar `WAKE_PROX`, baixar `WAKE_MAX`, ou adicionar gate
  `σ_a<0.85` no ponto vago além do geométrico.
- **Bug de propriedade não-inicializada** nas filhas (cf. `mean_c_n>max_c_n` da v2.3 /
  `is_filler` lixo): **todo** campo novo (`is_wake, x_dep, y_dep`) DEVE entrar
  explicitamente no dict `add_particles`. Checar no log que `is_wake ∈ {0,1}` e que
  `mass` não salta.

## 8. Invariantes a preservar (não violar)

- Não remover o **hard pin `rho_b≥0.8`** (K.17) nem o pin estrutural do C3.
- Não deixar o filler de rastro **produzir cs / crescer / sentir motor** (senão = C4 #36).
- Não **mover partícula existente** (senão = shifting #31).
- `h=h0` nas filhas (senão = colapso de dt #29/#38).
- Uma alavanca por vez: só `use_wake` muda vs baseline (§2.3).
