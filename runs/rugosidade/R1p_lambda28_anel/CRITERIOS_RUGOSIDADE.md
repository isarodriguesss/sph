# CRITÉRIOS PRÉ-REGISTRADOS — série de rugosidade R1-R4 (Pass L)

> Escrito em **2026-09-18, ANTES de qualquer rodada da série**. §2.3 do CLAUDE.md.
>
> Fica em `docs/` e não em `runs/<R>/CRITERIOS.md` porque `tools/archive_run.sh` faz `rm -rf`
> no destino e já apagou o CRITERIOS.md pré-registrado do P2R22. Copiar para cada run
> **depois** de arquivar.

## 1. A pergunta

**A rugosidade atrapalha o swarm?** É uma comparação, não um mapa de fases. Por isso o padrão é
fixo (pilares), a alavanca é a **dose** (`φ`, via `Λ`, com `A = 10 dx` constante), e a evidência é
a **monotonicidade ao longo de R1→R3**, não a diferença liso × rugoso isolada.

Com uma realização por dose, **uma resposta não-monotônica não autoriza conclusão** — exigiria
réplicas. Isso está pré-registrado para não ser decidido depois de ver o resultado.

## 2. Referência — P2R23, medido em 2026-09-18

| grandeza | valor |
|---|---:|
| `R99`(t=50) | **4.067** |
| **`dR/dt` em t ∈ [25,50]** | **0.0756** |
| largura EDT em 0.6R / ciclo | 5.4 dx / 0.61 |
| tortuosidade | 1.00-1.02 |
| dedos / AR | 15 / 9.9 |
| `a_mar_bio_p95` (mediana t>30) | 0.73 |
| `void_15` | 0.02% |
| `a_press_med` / picos > 4 | 3.22 / 31% |
| **`frac(c_s > 0.45)` nas vivas** | **0.007** |
| `frac(c_s > 0.45)` no corpo | 0.164 |
| **`max c_s` no ágar a < 2h do corpo** | **0.300** |
| frac. desse ágar acima de `COL_CS_MIN` = 0.3 | **0.0000** |

## 3. Predição (eixo 5 da análise preditiva)

Dois canais de estorvo, com **limiar de dose diferente** — é isso que torna a predição
falsificável, e não só "a rugosidade atrapalha".

> **CORRIGIDA em 2026-09-18 pelo pré-passo (seção 8).** Eu havia previsto R1 "transparente"
> raciocinando pela fresta estática (18 dx = 3.3 larguras de braço). Medido: **100% dos braços
> interceptam pelo menos um pilar nas TRÊS doses**, com 1.50 / 2.38 / 4.02 contatos por braço.
> A fresta diz que o braço *pode* passar entre dois pilares, não que passe — ao percorrer 53 dx
> radiais numa rede 2D ele encontra pilares necessariamente. O zero da escada de dose é o
> **run liso**, não R1.

| | `Λ` | **contatos/braço** | `φ` | **τ** | **`dR/dt`** | dedos |
|---|---:|---:|---:|---|---|---|
| liso | — | **0** | 0 | 1.00-1.02 | 0.0756 | 15 |
| **R1** | 28 dx | **1.50** | 11.6% | 1.05-1.15 | 0.064-0.072 (−5 a −15%) | 12-15 |
| **R2** | 20 dx | **2.38** | 22.7% | 1.15-1.30 | 0.053-0.064 (−15 a −30%) | 9-13 |
| **R3** | 16 dx | **4.02** | 35.4% | > 1.30 | 0.038-0.053 (−30 a −50%) | 6-10 |

As faixas de `dR/dt` e de dedos são estimativas ancoradas na contagem de contatos, não medições;
o que está de fato pré-registrado é a **ordenação monotônica** e o desacoplamento dos dois canais
abaixo. **O resultado científico é a distribuição de desfechos dos contatos** (desvio / divisão /
parada), não um alvo numérico — qualquer das três respostas é informativa.

1. **Mecânico** — o braço é desviado antes de ser barrado: **τ responde em R2 com `dR/dt` ainda
   dentro do ruído; `dR/dt` só cai em R3.** Se os dois caírem juntos já em R1, a predição erra e
   a causa provável é atrito (ver R4), não bloqueio.
2. **Químico** — os pilares são fluxo zero para `c_s`; ele acumula a montante e o gradiente à
   frente do líder enfraquece. Apareceria como queda de `a_mar_bio_p95` **já em R1-R2**, antes de
   qualquer efeito mecânico. Se aparecer, é resultado próprio: nenhum trabalho levantado na
   revisão acopla obstáculo ao campo de surfactante.
3. **Perda de dedo** tem mecanismo específico (lição #104): líder que **para** deixa de ser líder,
   e cada dedo é o rastro de um líder (#98). Um pilar que intercepta um líder elimina um dedo
   inteiro — não é afinamento, é extinção.

## 4. Métricas — decididas ANTES, com população e estatística explícitas (§2.5, lição #46)

**Todas em média ± desvio na janela t ∈ [35, 50]**, contra o ruído entre as duas realizações do
P2 como piso. Nenhuma contagem ou massa **absoluta**: a rodada rugosa nasce com −11.6% a −35.4%
de fluido, então só grandezas **por área** ou mascaradas.

- **Primário (responde à pergunta):** `R99`(t=50) e `dR/dt` em t ∈ [25,50].
- **Primário 2 — ANISOTROPIA AZIMUTAL, acrescentada em 2026-10-06 (seção 12):** amplitude e
  **fase** do modo m=6 de `R(theta)`. Janela e estatística próprias, definidas em 12.5 — a janela
  t ∈ [35,50] acima **não** vale para ela. `dR/dt` é média azimutal e é cega a isto: em R1 a
  modulação chega a 22% enquanto a média difere 3% do liso.
- **Secundário (como atrapalha):** área da colônia / disco, largura EDT e ciclo em 0.6R e 0.8R,
  contagem de braços, **tortuosidade**.
- **Mecanístico (por que):** fração de líderes que passa a menos de uma largura de braço (5.4 dx)
  de um pilar; quantos deixam de ser líder após o contato (rastreio de identidade da lição #104);
  desfecho classificado em **desvio / divisão / parada**.

**INVALIDADAS para a série (constatado em R1, 2026-09-21) — núcleo sólido, `R_min`/baía e a
largura EDT em 0.3R.** Elas medem `PILAR_R_EXCL`, não a colônia. A zona de exclusão é um disco
sem pilar de raio **1.238** no centro; em R1 o "núcleo sólido" saiu em 0.29 × Rmax = **1.151** —
ou seja o contorno que a métrica encontra é a borda interna da rede de pilares, e ele se move
quando `PILAR_R_EXCL` ou `Rmax` mudam, não quando o núcleo biológico muda. É a armadilha da
lição #101 (entrar na faixa da literatura mexendo no **denominador**) com um denominador novo.
Não usar nenhuma das três para comparar doses enquanto houver zona de exclusão.

## 5. Vigilância do confundidor químico — o risco principal

O `c_s` máximo no ágar vizinho ao corpo do P2R23 é **exatamente 0.300**, encostado no gate
`COL_CS_MIN = 0.3`, com fração acima do gate **0.0000**. Não há folga nenhuma: **qualquer
acúmulo local de `c_s` abre o gate da colonização**, a base é recrutada, o núcleo cresce e a
colônia vira disco — pelas lições #90, #94-H, #97 e #101, por razão **química**, não por bloqueio.

Em R3, `φ = 35%` derruba `D_eff` a ~55% por tortuosidade de meio efetivo. Como o `c_s` do mid-arm
já está a 0.98 do teto (#73), a subida se dá no ágar das baias — exatamente onde o gate está no
limite.

**Reportar obrigatoriamente em toda dose, junto do resultado primário:**

| sentinela | referência P2R23 | leitura |
|---|---:|---|
| `frac(c_s > 0.45)` nas vivas | 0.007 | subiu? o campo está saturando |
| frac. do ágar a < 2h do corpo acima de `COL_CS_MIN` | 0.0000 | **> 0.02 = gate aberto** |
| `biomass_arms` / `biomass_total` | 0.317 | **caiu? é disco** (material foi para o núcleo) |

> A terceira sentinela era `R_min` / núcleo sólido (0.39 R), **trocada em 2026-09-21** pela razão
> `biomass_arms`/`biomass_total` — o núcleo sólido foi invalidado pela seção 4. A razão é imune ao
> orçamento de partículas (a rodada rugosa nasce com menos fluido) e não precisa de máscara:
> se a colônia vira disco, o material se concentra em `rho_b` alto e a razão cai.

**Regra de interpretação, pré-registrada:** se `dR/dt` cair **e** o gate abrir **e** a razão
`biomass_arms`/`biomass_total` cair, a conclusão **não pode** ser "a rugosidade atrapalha o swarm" — é "a rugosidade abriu o
gate da colonização". Nesse caso a série tem de ser repetida com `COL_CS_MIN` reescalado, e a
regra da lição #78-b se aplica: o invariante ao reescalar um limiar absoluto é a **fração da
população que ele seleciona**, nunca a razão ao pico.

## 6. Critérios de ABORTO (matar a rodada, não esperar t=50)

1. `max_cs` > 0.6 na iteração 200 — instabilidade difusiva explícita (lição #96). `D` não muda
   nesta série, então não deve ocorrer; custa zero verificar.
2. `n_pen` > 0 sustentado por mais de um quadro — a força de contato não segura, e a geometria
   perdeu sentido.
3. `dt` médio < 0.5× o do liso no mesmo `t` — colapso de `dt`.
4. `max_v` > 1.0 ou partícula ejetada do domínio.
5. Massa com runaway (razão 2ª/1ª metade > 1.3, guardrail do §2.5).
6. **Pares a menos de 0.05 dx > 4 000, OU `nn_median` < 0.15 — ACRESCENTADO em 2026-10-06.**
   Guarda **absoluta**, de propósito: o critério 3 é relativo ao liso no mesmo `t` e **não pega**
   a instabilidade de pareamento, porque na validação de t=100 foi o **liso** que colapsou (seção
   13). Sem esta guarda, a mesma cascata numa rodada rugosa seria lida como "a rugosidade colapsou
   o `dt`". Medir com `tools/pares_dt.py`; o joelho do liso foi em 3 429 pares / `nn_median` 0.17.

> **`n_pen` é estrito — vigiar `n_pen_sup` junto (medido em 2026-09-18, ANTES de R1).** A
> definição pré-registrada de `n_pen` é `d(centro) < R_p − 0.5 dx`, e ela lê **0 em todos os
> seis quadros do V1** — a rodada que correu **sem** a guarda de deposição e cujo vazamento
> motivou a guarda. As 6 partículas que o plano §5 registrou como "dentro dos pilares" (5 filler)
> estão a `d/R_p` entre **0.971 e 0.998**: nos vãos do anel de superfície, não no corpo sólido.
> A distância mínima fluido–partícula-de-pilar no V1 é **0.799 dx**, acima do `INSERT_PROX` de
> 0.7 dx.
>
> A definição pré-registrada **fica como está** — mexer nela depois de ver o dado seria escolher
> a régua pelo resultado. Foi acrescentada a coluna `n_pen_sup` (`d(centro) < R_p`, = 6 no V1)
> como **indicador precoce**: é ela que sobe primeiro se material começar a atravessar. Leitura:
> `n_pen` é o aborto; `n_pen_sup` subindo entre quadros é o aviso de que o aborto vem.

## 7. Como o resultado será lido

- **Atrapalha:** `dR/dt` cai monotonicamente com `φ`, τ sobe, **e** as sentinelas químicas da
  seção 5 ficam na referência. A magnitude do bloqueio geométrico é **R3 − R4**, não R3 − liso.
- **Não atrapalha:** `dR/dt` dentro do ruído em **todas** as doses, **com** a fração de
  interceptação de líderes demonstrando que houve contato. Sem essa fração medida, um negativo
  não é publicável — seria a lição #102 de novo.
- **Inconclusivo:** resposta não-monotônica, ou sentinela química fora da referência.

## 8. Pré-passo — EXECUTADO em 2026-09-18

**(A) Empírico, V1 (Λ = 28 dx).** Distância dos líderes (`~filler & rho_b>=0.1 & r >= 0.7·R99`,
a mesma seleção do `_alarga_rastro`) à *superfície* do pilar:

| t | R99 | líderes | d p50 (dx) | d mín (dx) | **< 1 braço** |
|---:|---:|---:|---:|---:|---:|
| 8.2 | 0.958 | 63 | 9.9 | 4.2 | 3% |
| 12.4 | 1.176 | 33 | 7.3 | 1.3 | **42%** |
| 16.6 | 1.403 | 24 | 4.5 | 1.0 | **62%** |
| 20.0 | 1.682 | 23 | 5.9 | 1.3 | **48%** |

Com a colônia mal tendo passado da zona de exclusão (1.2), metade dos líderes já está em alcance
de contato.

**(B) Geométrico, Monte Carlo sobre a rede real** (4000 raios radiais de r=1.2 a 4.067; braço de
5.4 dx toca pilar quando o centro fica a < R_p + w/2). A rede é *ordenada*, então a estimativa de
Poisson por livre caminho médio não se aplica — daí o Monte Carlo:

| Λ | φ | pilares | fresta/braço | ≥ 1 contato | **contatos/braço** |
|---:|---:|---:|---:|---:|---:|
| 28 dx | 11.6% | 92 | 3.3 | **100%** | **1.50** |
| 20 dx | 22.7% | 180 | 1.9 | **100%** | **2.38** |
| 16 dx | 35.4% | 256 | 1.1 | **100%** | **4.02** |

> ### CORREÇÃO OBRIGATÓRIA (2026-10-06) — a tabela acima usa o alcance do LISO
>
> Os 4.067 são o `R99` do **P2R23**. Cada dose alcança **menos** que isso, e os contatos têm de ser
> contados no alcance **de cada uma**. Varredura determinística em ângulo (20 000 raios, sem ruído
> de Monte Carlo), `R_p` = 5 dx, braço 5.4 dx, r de 1.2 a `r_max`:
>
> | `r_max` | R1 (Λ=28) | R2 (Λ=20) | R3 (Λ=16) | escada R3/R1 |
> |---:|---:|---:|---:|---:|
> | 2.500 | 0.53 — 53.2% | 0.80 — 79.8% | 1.70 — 100% | 3.2× |
> | 3.000 | 0.91 — 83.7% | 1.36 — 100% | 2.27 — 100% | 2.5× |
> | 3.500 | 1.10 — 83.7% | 1.60 — 100% | 3.02 — 100% | 2.7× |
> | **3.852** (`R99` real do R1) | **1.10 — 83.7%** | 1.81 — 100% | 3.44 — 100% | 3.1× |
> | 4.067 (`R99` do liso — a tabela acima) | 1.50 — 100% | 2.22 — 100% | 3.84 — 100% | 2.6× |
> | 4.451 (liso em t=60) | 1.50 — 100% | 2.41 — 100% | 4.02 — 100% | 2.7× |
>
> **O R1, no alcance dele, tem 1.10 contatos/braço e 83.7% de interceptação — não 1.50 e 100%.**
> Causa: as coroas de Λ=28 estão em 1.508, 2.611, 3.015, **3.989**, 4.523, e o R1 parou **3.4%
> antes** da coroa √7Λ = 3.989, que é a que leva a interceptação a 100%. **Logo 16% dos braços do
> R1 nunca tocaram um pilar** — e é exatamente essa fração que a seção 7 exige para um
> quase-nulo ser publicável.
>
> **O que NÃO muda:** a escada de dose é robusta, 2.5-3.2× em todo o alcance de 2.5 a 5.0. R2 e R3
> têm 5 e 8 coroas no alcance e 100% de interceptação mesmo encolhendo muito — a quantização só
> morde a dose mais esparsa.
>
> **O problema estrutural do R1 é quantização, não tempo.** Com 3-4 coroas o contato é uma
> ESCADA com degrau de 0.4 contatos na borda do alcance, então a dose do R1 é sensível a poucos
> por cento de raio — a ordem do ruído entre realizações. Corrigir isso é mexer na **geometria**
> (Λ de 28 para ~26, ou `A` maior, pondo a coroa √7Λ dentro do alcance), **nunca** na janela
> temporal: estender para t=60 cruzaria a coroa, mas compararia contra um liso já degradado
> (seção 13).

**Conclusão (revista em 2026-10-06):** o risco de falso negativo por transparência está afastado
para **R2 e R3** (100% de interceptação em qualquer alcance plausível) e **parcialmente para o R1**
(83.7%). A escada de dose vale 2.5-3.2×. O zero é o run liso.

Script: `scratchpad/interceptacao.py` (análise de geometria sobre run existente; não toca o solver).

## 9. Estado da infraestrutura (2026-09-18)

Bloqueantes fechados e verificados: máscara de pilar em `void_fraction`, em `fig_tese.campo` (e
por herança em `rank_runs`, `lit_rank`, `diag_disco`, `cmp_colonia`) e em `ciclo_largura`; guardas
geométricas em `_insere_vacuo`, no `livre()` do wake e no `ParticleShift`. Com `use_pilares=False`
o código é **bit-a-bit idêntico** ao anterior (`runs/swarm/_verif_mask_novo` × `_verif_mask_ctrl`,
10 threads, iterações 0/200/400/600).

**Instrumentação — FEITA em 2026-09-18:** colunas `n_pen`, `n_pen_sup`, `n_contato`, `a_rep_max`,
`a_rep_med` no FIM do `LOG_HEADER` (47 → 52 colunas), mais os arrays `ax_rep`/`ay_rep` acumulados por `ForcaContornoPilar` e **descontados** do resíduo `a_press_max`/`a_press_med` (plano §4: senão a coluna de pressão viraria a coluna de contato). `n_pen` também é impresso a cada quadro no console, com alerta, para o critério de aborto 2 ser verificável em tempo de execução e não só no CSV.

---

## 10. R1 — Λ = 28 dx, φ = 11.6% — EXECUTADO em 2026-09-21

`runs/rugosidade/R1_lambda28`, 2863 iterações, 1429 s, 16 frames, t = 50, SEED = 20260806,
10 threads. Fluido 61 127 partículas (68 121 − 6 994 de pilar). 92 pilares de `A` = 10 dx.

### 10.1 Aborto — os cinco critérios passam

| critério | limite | medido |
|---|---|---|
| 1. `n_pen` (fluido dentro do corpo do pilar) | qualquer > 0 | **0 em todos os 16 frames** |
| 2. `max_cs` > teto | > 0.5 já na iteração 200 | 0.4939 |
| 3. `max_v` | > 1.0 | 0.429 |
| 4. massa | runaway | +3.2% |
| 5. colapso de `dt` | ≫ 62 iter/t | 57.3 iter/t |

A força de contato de Monaghan & Kajtar segura: `a_rep_max` chega a 2.9, mesma ordem de `f0` = 3.0,
e nenhuma partícula cruza meio `dx` para dentro do sólido. `n_pen_sup` (centro além da superfície
nominal, ainda fora do corpo) sobe monotonicamente 0 → 43 — é a compressão elástica esperada do
kernel de contato, não penetração.

### 10.2 Primário — a expansão quase não muda

| | liso (P2R23) | R1 | Δ |
|---|---:|---:|---:|
| `R99`(t=50) | 4.067 | 3.852 | **−5.3%** |
| `dR/dt` em t ∈ [25,50] | 0.0756 | 0.0733 | **−3.1%** |

A predição pré-registrada era −5% a −15% em `dR/dt`. O medido (−3.1%) fica **abaixo da banda** e
plausivelmente dentro do ruído entre realizações.

> **QUALIFICAÇÃO (2026-10-06, seção 8):** recontados no alcance real do R1 (`R99` = 3.852, e não
> no 4.067 do liso), os contatos são **1.10 por braço com 83.7% de interceptação** — **16% dos
> braços nunca tocaram um pilar**. O −3.1% é portanto um quase-nulo medido sobre uma dose
> **parcialmente transparente**, não sobre contato universal. Pela seção 7 isso não invalida o
> ponto, mas obriga a reportar a fração junto do número. **Na dose mais baixa, a rugosidade não atrapalha
a expansão de forma destacável do ruído** — o que é exatamente por que a evidência da série é a
**monotonicidade R1→R3**, e não este par isolado (seção 7).

### 10.3 Sentinelas químicas — não dispararam, e o campo MELHOROU

| sentinela | liso | R1 | leitura |
|---|---:|---:|---|
| frac. do ágar a < 2h acima de `COL_CS_MIN` | 0.0000 | **0.0000** | gate fechado |
| máx. `c_s` no ágar vizinho | 0.300 | **0.300** | idêntico, ainda exatamente no gate |
| `frac(c_s > 0.45)` nas vivas | 0.0065 | 0.0174 | subiu, longe de saturar |
| `biomass_arms`/`biomass_total` | 0.317 | **0.519** | **subiu** — não é disco |

E o campo foi na direção da literatura: `c_s` no corpo 0.44 → **0.62**, nas baias 0.25 → **0.34**,
halo `L` 0.17 → **0.18 R — o alvo exato de [T11]/[T14]**. O confundidor que motivou a seção 5 não
se materializou em φ = 11.6%; ele continua sendo o risco a vigiar em R3 (φ = 35%).

### 10.4 O resultado estrutural, que não estava previsto

| | liso | R1 |
|---|---:|---:|
| portadores na banda de swarmer (`rho_b` ∈ [0.1,0.5)) | 469 | **1053** |
| `biomass_arms` | 0.209 | **0.503** |
| `biomass_total` | 0.658 | **0.970** |
| largura EDT em 0.6R | 5.4 dx | **7.3 dx** |
| ciclo em 0.6R | 0.61 | 0.78 |
| dedos | 15 | 14 |
| AR | 9.9 | 7.7 |
| `a_mar_bio_p95` | 0.78 | 0.56 |
| `void_15` (mascarado) / `mean_sig_all` | 0.0002 / 0.998 | 0.0000 / **1.001** |

R1 tem **2.25× mais portadores na banda ativa e 2.4× mais biomassa nos braços**, contra um
orçamento de partículas **10.3% menor** — a direção é o oposto do que o orçamento explicaria, e
`arms_mask` é uma banda de **densidade** (`rho_b` ∈ [0.1,0.5)), não uma região radial, então não há
o confundidor de `R99` menor. O suporte de kernel fica perfeito (`mean_sig_all` = 1.001).

**Leitura:** em φ = 11.6% os pilares não bloqueiam — eles **retêm**. Braços 35% mais largos, mais
densos, com motor mais fraco por partícula (`a_mar_bio_p95` −28%) e a mesma velocidade de frente.
Material que no liso avançaria e deixaria rastro fino, aqui fica. Isso é consistente com o
pré-passo (100% dos braços interceptam ≥ 1 pilar, 1.50 contatos/braço): o contato existe em todos
os braços e o efeito dele, nesta dose, é de espessamento, não de parada.

**Não generalizar para as doses altas.** Espessar o braço em φ baixo e bloquear em φ alto são
compatíveis; é justamente a curva que R2 e R3 têm de resolver. E o sinal a vigiar mudou: se o
espessamento continuar em R2/R3, a colônia caminha para disco por via **mecânica** (retenção), que
as sentinelas da seção 5 — todas químicas, exceto a razão nova — não detectam.

### 10.5 Defeito de instrumentação corrigido

`a_rep_med` saiu em 5.2e-42 porque era a mediana sobre `d_sup < dx`, conjunto de ~2285 partículas
dominado por ágar estático que sente força **exatamente zero** (o kernel de contato tem alcance
1 dx). Corrigido para a mediana sobre `a_rep > 0`. `n_contato` foi **mantido geométrico** de
propósito, para R1 seguir comparável com R2-R4 — `ax_rep` não vai para o HDF5, então o valor de R1
não pode ser recomputado offline. **`a_rep_med` só é válido a partir de R2**; `n_pen`, `n_pen_sup`,
`n_contato` e `a_rep_max` valem desde R1.

---

## 11. Ancoragem de literatura do R4 (§10) — acrescentada em 2026-09-21

O R4 (controle de atrito casado) estava sem forma funcional justificada: "subir `γ` até casar" é
tentativa e erro, que o §10 rejeita em revisão. **[T15] Verma et al. 2025** fornece a âncora.

No modelo deles o atrito **não é força** — é amortecimento nodal numa Langevin superamortecida,
com o coeficiente escalando pela **área de contato deformada**:

```
ζ(t) = ζ(0) · ( l(t) / l(0) )²
```

Dois motivos para isso importar aqui:

1. **Forma funcional.** O nosso `LinearDrag` usa `−(γ + 1.5·γ·ρ_b²)·v`. A dependência em `ρ_b²`
   nunca teve fundamentação publicada — é calibração herdada do Pass I.4. A forma de [T15] é
   área de contato, não densidade. Ao fixar o R4, declarar qual das duas está sendo usada e por quê.
2. **Aviso sobre o emaranhamento.** Os autores afirmam que atrito e adesão "são provavelmente
   correlacionados positivamente, mas a relação exata é desconhecida". É exatamente o nosso caso:
   o pilar é obstrução geométrica **e** fonte de `ViscousForce` ao mesmo tempo. Reforça a regra já
   pré-registrada na seção 7 — a magnitude geométrica é **R3 − R4**, nunca R3 − liso.

**O que [T15] NÃO autoriza:** o regime deles é um filme elástico aderido que acumula tensão até
flambar; o nosso é expansão advectiva superamortecida sem adesão, com a pressão da EOS 219× abaixo
da Marangoni (lição #85). A mesma resistência no substrato vira **tensão elástica** lá e **retenção
de material** aqui — que é a leitura correta do espessamento de R1 (seção 10.4). Não importar
conclusões sobre tensão crítica, raio crítico ou contagem de rugas.

**Prior qualitativo, não predição:** [T15] mede que em **adesão baixa** o biofilme é praticamente
insensível à heterogeneidade do substrato. O nosso modelo não tem termo de adesão nenhum — estamos
no extremo desse eixo, o que é coerente com o −3.1% de R1. Mas a heterogeneidade deles é em energia
de adesão e a nossa é obstrução geométrica: eixos diferentes, transposição sugestiva apenas.

> **Nota de arquitetura (registrada para não se perder):** no nosso regime superamortecido
> (`u = F/γ_eff`, lição #86-C), adesão e atrito seriam **degenerados** — os dois viram o mesmo
> coeficiente multiplicando `v`. A adesão só ganha assinatura própria quando há movimento fora do
> plano para ela resistir, que é o caso de [T15] e não o nosso. É a mesma degenerescência que
> aposentou a Opção 1 da análise arquitetural (campo de atrito contínuo ≡ `MOTOR_SCALE`, reprovado
> em #86-G). Acrescentar adesão a este modelo, como ele está, não acrescentaria física —
> acrescentaria um segundo botão para o mesmo efeito.

---

## 12. ANISOTROPIA AZIMUTAL — quarta métrica, acrescentada em 2026-10-06 (antes do R2)

### 12.1 A lição: a métrica primária escondeu o maior efeito do R1

`dR/dt` e `R99` são **médias azimutais**. Medido no R1: a modulação radial do modo m=6 chega a
**0.223 (22%)**, com `R(theta)` indo de ~1.17 a ~1.83 em t=25.6 — **57% de espalhamento entre a
direção mais rápida e a mais lenta** — enquanto a média difere 3% do liso. O R1 **não é uma
colônia igual um pouco mais lenta; é uma colônia hexagonalmente organizada cujo raio médio por
acaso ficou parecido.**

É a lição #101 por outro caminho: a régua escolhida não enxergava o fenômeno. Diferente de #101,
aqui a régua não estava errada — estava **incompleta**, e a correção é acrescentar uma métrica,
não trocar a existente.

Origem do enquadramento: sugestão do Prof. Cesar (UFPE, 2026-10-06) de descrever o transporte como
**caminhada enviesada sobre grafo com regiões proibidas**. Não é analogia: o líder É um caminhante
isolado (licao #98 — ~20 vivas além de 0.85 R, uma por dedo, a viva mais próxima a 14-21 dx, fora
do alcance de toda força do modelo, 3h = 5.4 dx), o viés é `-beta*grad(c_s)` radial, e o pilar é
região proibida verificada (`n_pen` = 0 nos 16 quadros). **Regime superamortecido** (`u = F/gamma`,
licao #86-C): transfere tortuosidade, livre caminho médio, percolação e travamento em direções da
rede; **não** transfere bilhar de Sinai nem difusão anômala balística.

### 12.2 A medição — R1 contra o liso (`tools/anisotropia.py`)

| t | R_med | a₆ liso | fase liso | a₆ R1 | **fase R1** |
|---:|---:|---:|---:|---:|---:|
| 16.6 | ~0.96 | 0.0365 | 3.7° | 0.0936 | **27.5°** |
| 20.0 | ~1.14 | 0.0281 | 2.2° | 0.1882 | **27.2°** |
| 25.6 | ~1.50 | 0.0142 | 2.8° | **0.2231** | **29.0°** |
| 29.1 | ~1.70 | 0.0119 | 7.1° | 0.2002 | **29.4°** |
| 32.7 | ~1.92 | 0.0142 | 15.3° | 0.1631 | **30.0°** |
| 39.2 | ~2.28 | 0.0084 | 38.0° | 0.0685 | **28.8°** |
| 46.2 | ~2.56 | 0.0272 | 40.6° | 0.0134 | **29.0°** |

**Fase travada em 27-30° por nove quadros consecutivos**, contra fase que passeia de 2° a 42° no
liso. Em t=25, a₆ do R1 é **16x** o do liso, e m=6 é o **modo dominante** do espectro de t=20 a
t=33. Depois passa para **m=12 = 2x6**, que cresce até 0.200.

**E a rede SUBSTITUIU o modo do modelo, não somou um.** O liso tem modo dominante próprio em
m=10-11 (0.075 -> 0.179, crescente); no R1 esse modo cai para 0.079-0.102.

### 12.3 O mecanismo, e por que a fase é a prova

A fase não é ajustável — é fixada pela geometria, incluindo `PILAR_R_EXCL`. Com rotação 0 e o
inóculo num sítio **vago**, a colônia enfrenta pilar de frente nos azimutes da primeira coroa
**presente** e garganta a 30° dela:

| | Λ | 1ª coroa presente | azimutes (mod 60°) | **máximos previstos** | medido |
|---|---:|---|---:|---:|---:|
| **R1** | 1.5077 | **coroa 1**, r = 1.0 Λ = 1.5077 | 0° | **30°** | **29°** ✓ |
| **R2** | 1.0769 | **coroa 2**, r = 1.732 Λ = 1.8653 | 30° | **0°** | — |
| **R3** | 0.8615 | **coroa 2**, r = 1.732 Λ = 1.4922 | 30° | **0°** | — |

Em R2 e R3 a primeira coroa (`r` = Λ = 1.077 e 0.862) cai **dentro** da zona de exclusão de 1.2 e
é removida; sobra a segunda coroa, que na rede triangular está a 30° da primeira.

> ### PREDIÇÃO PRÉ-REGISTRADA, falsificável e não-ajustável
>
> **A fase do m=6 deve VIRAR 30° de R1 para R2/R3: de ~29° para ~0° (mod 60°).**
>
> Se a fase do R2 sair perto de 0°, o travamento na rede está confirmado além de dúvida — nenhum
> ajuste produz uma inversão de fase prevista a partir de qual coroa sobrevive à exclusão. Se
> permanecer em ~30°, **a fase é fixada por outra coisa e o mecanismo está errado** — e aí a
> primeira suspeita é a semente azimutal `cos(8*theta)`, que não tem simetria de 6 dobras e
> portanto não explicaria nem o R1, mas tem de ser descartada por medição.
>
> Corolário: a amplitude deve **picar quando `R_med` cruza a 1ª coroa presente** — R1 picou em
> `R_med` = 1.499 contra coroa em 1.5077 (razão 0.99). Logo o pico do R2 é esperado em
> `R_med` ~ 1.87 e o do R3 em ~1.49 — **mais tarde** no R2 que no R1, apesar da dose maior.

### 12.4 A métrica — população, estatística e JANELA

- **População:** colônia = `rho_b >= 0.1` ou filler (§2.2), `R(theta)` das **partículas**, não do
  raster (licao #92), 720 bins.
- **Reportar sempre o par (amplitude, fase).** Amplitude sozinha não distingue assinatura da rede
  de flutuação do modelo — é a fase travada que distingue.
- **JANELA PRÓPRIA — a de t ∈ [35,50] do §4 NÃO vale aqui.** Nela o sinal do R1 já decaiu
  (0.040 ± 0.037 contra 0.022 ± 0.011 do liso: marginal), e usá-la **esconderia o efeito**. O
  imprint é um evento geométrico, não um regime estacionário. Janela pré-registrada:
  **`R_med` ∈ [0.7, 1.3] × raio da 1ª coroa presente**, com o pico e a trajetória completa
  reportados junto.
- **Piso de fase:** a fase só é lida onde `a₆ >= 0.03` (ou metade do pico). No R1 ela é estável em
  27-30° em todos os quadros acima disso e salta para 44-48° exatamente quando `a₆` cai a 0.007 —
  fase de amplitude nula é ruído, não resultado.
- **Régua calibrada antes de ranquear** (licoes #67-F, #101): `python tools/anisotropia.py
  --autoteste` recupera amplitude e fase exatas de sinais sintéticos m=4/6/12. Rodar antes de citar
  qualquer valor.

### 12.5 Consequência para a licao #107 — e é o resultado, não a correção

A licao #107 mediu que **o modelo não seleciona comprimento de onda**: `sigma(m)` sem máximo, 22
dos 40 modos com `sigma` negativo, `m*` ~ 13.8 fixo enquanto `R` cresce 4.55x, semente m=8 apagada
(0.293 -> 0.036). A conclusão foi que a contagem de dedos é fixada pelo número de portadores na
frente, que é número de **discretização**, não escala física.

**Em R1 a geometria externa IMPÕE o modo, com a fase fixada pela rede.** Isso não conserta a
limitação — **localiza-a**: a ausência de seleção é interna ao motor (teto de `c_s`, lições #64,
#77, #106), e um substrato estruturado a contorna por fora. Para a tese, a licao #107 deixa de ser
limitação confessada e passa a ser **resultado com controle**: medimos a ausência de seleção
interna e medimos a imposição externa de escala no mesmo modelo, pareado.

### 12.6 O que falta medir (pós-processamento, sem rodada nova)

1. **Classificação de desfecho de contato** (desvio / divisão / parada) nas trajetórias do R1 — é o
   que a §7 exige para um resultado quase-nulo ser publicável, e **não foi rodado**.
2. **Distribuição de `Delta_theta` COM SINAL.** Média diferente de zero = deriva transversal, a
   "quiralidade" do enquadramento. Mesma medição que o Chang usa e que o Postek motiva (velocidade
   azimutal): três motivos na mesma conta.
3. **Expoente de MSD dos líderes**, `<r²> ~ t^alpha` — separa deriva pura de deriva com
   espalhamento.
4. **Previsão para o padrão sorteado (padrão 3):** sem direções privilegiadas, a anisotropia deve
   **desaparecer** na mesma `phi`, com tortuosidade parecida. Dá ao padrão 3 uma previsão
   falsificável em vez de "vamos ver".

---

## 13. VALIDAÇÃO DO LISO ATÉ t=100 — INTERROMPIDA em t=66.7 (2026-10-06)

`runs/swarm/P2R23_t100_pareamento` (42 quadros, `main.py` e `src/` arquivados). Baseline P2R23,
`use_pilares = False`, `KERNEL = "cubic"`, `H_FACTOR = 1.8`, `total_sim_time = 100`. Morta em
`iter` = 8200, `t` = 66.72, após 1h09 de parede.

### 13.1 Não foi nenhum dos modos de falha que a seção 6 vigia

| sinal | numa falha típica | medido |
|---|---|---|
| `max_v` | sobe | **cai**: 0.097 → 0.043 |
| `n_fast` | sobe | **0** |
| `a_press_med` | sobe | 0.005 |
| `max_cs` | fura o teto | 0.4938, cravado |
| massa | runaway | +7%, linear |

**Os cinco critérios de aborto da seção 6 passariam todos.** Daí o critério 6 novo.

### 13.2 O que foi: instabilidade de pareamento (Liu §6.4)

| t | `dt`/passo | pares < 0.05 dx | `frac_clump` | `nn_median` |
|---:|---:|---:|---:|---:|
| 38.7 | 1.6e-2 | 463 | 0.53 | 0.47 |
| 53.6 | 5.4e-3 | 1 468 | 0.63 | 0.36 |
| 60.6 | 5.5e-3 | 3 429 | 0.77 | 0.17 |
| **62.4** | **1.3e-3** | — | 0.80 | 0.13 |
| 65.6 | 1.4e-3 | **12 243** | 0.88 | 0.020 |
| 66.7 | — | — | 0.90 | **0.006** |

Pares dobrando a cada ~2 unidades de `t`; `iter/t` acumulado **112.8** contra 62 do P2R23 a t=50.
Joelho em **t ≈ 61-62**. Não é falta de repulsão: o filler está em `rho_b` = 0.477, e na
`BiomassEOS` isso dá `s = t²(3−2t)` com `t` = 0.9425, ou seja **`fade_rep` = 0.99** — repulsão
praticamente plena, com pressão presente (`rho/rho0` p99 = 1.87). O que falha é o **kernel**: no
spline cúbico `DWIJ` → 0 quando `r` → 0, então a força entre o par desaparece justo quando eles se
encostam, e `dt_cfl ~ h·|v·r|/r²` estoura.

### 13.3 A população é FILLER — corrige a atribuição da lição #88

| em t=65.6, partículas em pares < 0.05 dx | n |
|---|---:|
| **filler** | **10 465 (98.2%)** |
| vivas ≥ 0.1 | 170 |
| das quais recrutas (0.1 ≤ `rho_b` < 0.5) | 135 |

A lição #88 atribuiu a cascata do P2 às **vivas recrutadas**. Aqui 98% são **filler**
(`rho_b` mediano 0.458). E a fonte é o **rastro**: filler 2 382 (t=21.7) → **17 444** (t=65.6),
enquanto o ágar cai 65 400 → 57 005 (**8 395 convertidos**). O mecanismo que dá continuidade aos
braços é o mesmo que mata a rodada longa.

### 13.4 A janela numérica do P2R23 é t ≲ 60, não t ≲ 75

O `t ≲ 75` do CLAUDE.md foi medido no **P2**, que tinha 779 pares em t=70. O P2R23 tem **12 243 em
t=65.6** — ~16× pior e mais cedo. **O rastro encurtou a janela numérica em ~15 unidades de `t`.**

### 13.5 E a colônia para de expandir em t ≈ 54, antes do colapso

| janela | `dR/dt` |
|---|---:|
| t ∈ [25, 50] | 0.0756 |
| **t ∈ [54, 66.7]** | **0.0149** |

`R99` = 4.455 (t=53.6) → 4.650 (t=66.7): **+4.4% em 13 unidades de `t`**.

### 13.6 Consequência para a série: a janela da série é t=50, e não se estende

Três razões, independentes:

1. **Não há expansão a comprar.** De t=50 a t=60 o `R99` do liso sobe 11%, quase tudo antes de
   t=54; depois disso a colônia está parada (13.5).
2. **O trecho t ∈ [54,60] está contaminado nas duas pontas** — motor morto (`n_fast` = 0) e
   pareamento em curso (1 468 pares em t=53.6, 3 429 em t=60.6).
3. **Toda métrica da série é rugoso × liso no mesmo `t`.** Comparar em t=60 usaria um baseline
   degradado: seria a mesma informação com referência pior.

> **O R1 SUPORTARIA t=60 — e isso não é argumento para estender.** Medido em t≈49: R1 tem
> `frac_clump` 0.482 e `nn_median` 0.512, contra 0.546 e 0.451 do liso em t=48.3 — o R1 está
> ~7 unidades de `t` **atrás** do liso no pareamento, porque tem menos fluido (pilares removidos) e
> menos filler. O limitante é o **baseline**, não a rodada rugosa.

**Se o comportamento numérico de longo prazo for necessário**, a saída é `KERNEL = "wendland_c2"`
com `H_FACTOR = 1.92`: já está no código e foi validado (`runs/swarm/P2W_t100` chegou a **t = 99.6
em 9 200 iterações**, com 7× menos pares — lição #89). A decisão de 2026-09-14 de ficar no cúbico
foi tomada **para o P2**, aceitando t ≲ 75; essa troca foi precificada numa janela que o P2R23 não
tem. **Se o que se quer é morfologia de longo prazo**, o t=100 do P2R23 não responde: pelas lições
#90/#92 a colônia já é disco a partir de t≈60-70.

---

## 14. DECISÕES FECHADAS ANTES DO R2 (2026-10-06, decisão da usuária)

### 14.1 O confundidor de LARGADA — `PILAR_R_EXCL` não é neutro na dose

`PILAR_R_EXCL` = 1.2 **não remove nada em R1** (Λ = 1.5077 > 1.2) e **remove a primeira coroa em
R2 e R3** (Λ = 1.0769 e 0.8615, ambos < 1.2). Como na rede triangular a 2ª coroa está a √3Λ, o raio
do primeiro obstáculo fica:

| | Λ | coroa 1 removida | 1ª coroa presente | **1º contato possível (superfície)** | colônia chega lá |
|---|---:|---|---:|---:|---:|
| **R1** | 1.5077 | não (0 pilares) | 1.5077 = 1.000 Λ | **1.2385** | t ≈ 21 |
| **R2** | 1.0769 | **sim (6 pilares)** | 1.8653 = 1.732 Λ | **1.5961** | **t ≈ 27** |
| **R3** | 0.8615 | **sim (6 pilares)** | 1.4922 = 1.732 Λ | **1.2230** | t ≈ 21 |

**Ordem do primeiro contato: R3 (1.223) ≈ R1 (1.239) ≪ R2 (1.596).** Não-monotônico na dose —
**R2, a dose intermediária, tem a largada mais livre das três**, com ~6 unidades de `t` a mais sem
obstáculo que as vizinhas.

Isto ataca a lógica da seção 1 (a evidência é a monotonicidade R1→R3, e uma resposta
não-monotônica não autoriza conclusão). É a mesma doença da seção 4 — métrica que mede
`PILAR_R_EXCL` em vez da colônia — num lugar novo: aqui ela não contamina uma régua, contamina a
**condição inicial** de cada dose.

**DECISÃO: manter `PILAR_R_EXCL` = 1.2** (opção (a) de quatro avaliadas). As alternativas — baixar
para < 0.8615, escalar com Λ, ou comprimir a escada para Λ todos > 1.2 — **invalidariam o R1**, que
já rodou com 1.2, e a última mataria o vão de dose (φ iria só até 18.8% em vez de 35.4%).

**Severidade medida, que é o que justifica a decisão:**

- **Leve no primário.** Em t ∈ [25,50] o R2 fica livre apenas nas 2 primeiras das 25 unidades. E a
  dose **total** preserva a ordem: 1.10 < 1.81 < 3.44 contatos/braço (seção 8).
- **Séria na anisotropia.** O pico ocorre em `R_med` ≈ raio da coroa — 1.508 / **1.865** / 1.492 —
  então a impressão do R2 acontece numa colônia **maior**, com mais dedos já formados competindo
  com o modo 6. Um a₆ menor em R2 poderia ser geometria da largada, não dose.

> ### PRÉ-REGISTRO — a anisotropia entre doses é comparada em `R_med` / coroa, nunca em `t` casado
>
> A impressão da rede é um **evento geométrico** disparado quando a frente atravessa a 1ª coroa
> presente, não um regime em `t`. Portanto:
>
> 1. **Abscissa da comparação entre doses: `R_med / r_coroa`**, com `r_coroa` = 1.5077 (R1),
>    1.8653 (R2), 1.4922 (R3). O pico do R1 ocorreu em `R_med` = 1.499, ou seja
>    `R_med/r_coroa` = **0.99** — e é esse 0.99 que as outras doses têm de reproduzir.
> 2. A janela da seção 12.4 (`R_med` ∈ [0.7, 1.3] × coroa) **já está nesta unidade** e não muda.
> 3. **Comparar a₆ em `t` casado entre doses é PROIBIDO** nesta série — em t=25 o R1 está no pico
>    e o R2 nem alcançou a primeira coroa.
> 4. **Se o R2 sair não-monotônico em qualquer métrica**, a primeira verificação obrigatória é o
>    artefato de largada: refazer a leitura em `R_med / r_coroa` e checar se a anomalia sobrevive.
>    Só então a seção 1 se aplica (réplicas).

### 14.2 `ax_rep` / `ay_rep` no HDF5 — a partir do R2

**DECISÃO: instrumentar.** Acrescentados ao `add_output_arrays` do fluido em
[main.py](../main.py). Sem eles não se pode recomputar offline **quem** sentiu contato, o que é o
que separa "desviou por causa do pilar" de "desviou por outro motivo" na classificação de desfecho
da seção 4.

**O R1 não os tem** — assimetria aceita e registrada, como já ocorre com `n_contato` (seção 10.5).
A alternativa seria ficar cego em R2, R3 e R4 para preservar simetria com uma cegueira.

### 14.3 Ordem de trabalho: desfecho de contato no R1 ANTES do R2

**DECISÃO: rodar primeiro.** Razão de pré-registro, não de conveniência: definir a métrica
mecanística **depois** de ver o R2 é escolher a régua com o resultado na mão — o erro das lições
#68 (métrica circular) e #101 (régua persistida sem validação). Exercitá-la no R1 fixa a régua e
produz a distribuição de referência contra a qual o R2 será lido.

### 14.4 O R1 fica como está

**DECISÃO: manter.** Os 83.7% de interceptação (seção 8) **não são erro, são propriedade da dose**
— é o degrau de *início de contato*, e isso é informativo. Refazer com Λ ≈ 26 compraria
interceptação cheia ao custo de encurtar a escada de dose e de uma rodada.

**Consequência para o relato:** a fração de interceptação (83.7%) tem de aparecer **junto** do
−3.1%, sempre. Pela seção 7, um quase-nulo sem a fração medida não é publicável; com ela, é um
quase-nulo sobre dose parcialmente transparente, que é um resultado diferente e mais fraco.

---

## 15. DESFECHO DO CONTATO NO R1 — a métrica mecanística, EXECUTADA em 2026-10-06

`tools/desfecho_contato.py runs/rugosidade/R1_lambda28`. Pós-processamento dos 16 quadros; régua
fixada **antes** de olhar o resultado (decisão 14.3). Identidade = índice do array, verificada: `n`
cresce monotonicamente 61 127 → 62 929 e o deslocamento máximo de índice comum entre quadros é
5-7 dx, compatível com movimento real.

**461 observações líder-quadro.** O controle é interno e é o ponto: líderes param também **sem**
pilar (lição #104), então a medida é P(parada | contato) contra P(parada | livre) na mesma rodada.

### 15.1 PARADA — o desfecho dominante, e atribuível

| limiar de contato | n | partículas distintas | **P(parada \| contato)** | P(parada \| livre) | razão | `dr` cont/livre |
|---:|---:|---:|---:|---:|---:|---:|
| **1.0 dx** (alcance do kernel) | 11 | 7 | **90.9% ±8.7** | 8.4% | **10.8×** | **0.12** |
| 1.5 dx | 19 | 12 | 63.2% ±11.1 | 8.1% | 7.8× | 0.13 |
| 2.0 dx | 30 | 17 | 43.3% ±9.0 | 8.1% | 5.3× | 0.29 |
| **3.0 dx** | **61** | **25** | **24.6% ±5.5** | 8.2% | **3.0×** | 0.71 |
| 4.0 dx | 96 | 32 | 20.8% ±4.1 | 7.7% | 2.7× | 0.90 |
| 5.0 dx | 126 | 34 | 15.9% ±3.3 | 8.4% | 1.9× | 1.02 |
| 7.0 dx | 190 | 35 | 12.1% ±2.4 | 9.2% | 1.3× | 1.17 |
| 10.0 dx | 311 | 55 | 10.6% ±1.7 | 10.0% | 1.1× | 1.30 |

**A forma da curva é a prova, não o valor isolado** (lição #68 — a varredura do parâmetro da métrica
é obrigatória): a razão **cresce monotonicamente quando o limiar aperta** e converge a 1 quando ele
abre, enquanto **o controle fica achatado em 7.7-10.0% em todos os limiares**. É assinatura de
efeito de curto alcance real; se fosse composição (líderes em região densa param mais), a razão não
escalaria com a distância ao pilar.

O ponto de operação estatístico é **3 dx (n=61, 25 partículas distintas, 3.0× ±5.5)**; a 1 dx o
efeito é 10.8× mas com n=11 de 7 partículas.

**E o líder que para deixa de ser líder: 45.5% ±15.0 contra 15.8% ±1.7.** É o mecanismo das lições
#104 e #98 medido diretamente — o limiar `r ≥ 0.7·r99` sobe com a colônia, então parar é sair do
conjunto, e cada dedo é o rastro de um líder.

### 15.2 DESVIO — praticamente ausente

| | em contato | livre |
|---|---:|---:|
| `\|arco\|` azimutal (p50) | 0.48 dx | 0.40 dx |
| avanço radial `dr` (p50) | **0.39 dx** | **3.24 dx** |

O arco azimutal muda 20% enquanto o avanço radial cai **8×**. **O líder não contorna o pilar — ele
para na frente dele.** Isso contraria a predição 1 da seção 3, que supunha "o braço é desviado antes
de ser barrado"; em φ = 11.6% o desfecho é barramento, não deflexão.

**Quiralidade: não detectada.** Arco com sinal +0.32 ±0.37 dx em contato contra −0.08 dx livre — o
erro é maior que o valor. Negativo honesto, e responde a sugestão do Prof. Cesar com a estatística
disponível.

### 15.3 DIVISÃO e FUSÃO — nível de grupo (ligação 6 dx = `RASTRO_PONTA_LINK`)

| | total | com pilar a < 3 dx |
|---|---:|---:|
| divisões | 13 | 6 |
| **fusões** | **0** | — |
| grupos extintos | 10 | 3 |
| grupos novos | 2 | — |

**As 13 divisões ocorrem todas entre t=8 e t=29** — é a fragmentação inicial da colônia, de 1 grupo
para ~23, não cisão por pilar. **Depois de t=29: zero divisões.** E **zero fusões em todo o run**,
coerente com as lições #106/#107 (não há competição lateral entre dedos).

### 15.4 A síntese — e é ela que fecha a §7 para o R1

O efeito do pilar sobre um líder, em φ = 11.6%, é **PARAR**, não desviar nem dividir. Decisivo
quando acontece (10.8× a 1 dx, avanço a 12%), **e raro**: 11 de 461 observações líder-quadro a 1 dx,
e apenas 3 das 10 extinções de grupo tiveram pilar a menos de 3 dx.

**Isso explica mecanicamente o −3.1%**, que é o que a seção 7 exigia para o quase-nulo ser
publicável: a dose mais baixa **não é transparente** (83.7% dos braços interceptam, seção 8) e o
contato **é severo** quando ocorre — mas a frequência é baixa o suficiente para a população de dedos
quase não mudar (23 grupos em t=29.1 → 17 em t=50). Não é "a rugosidade não afeta o líder"; é "a
rugosidade mata o líder que ela encontra, e em φ = 11.6% ela encontra poucos".

> **Predição para R2 e R3, decorrente:** se o mecanismo for este, a razão P(parada|contato) deve
> ficar **aproximadamente constante** (é propriedade do contato, não da dose) e o que deve crescer é
> a **frequência de contato** e, com ela, a taxa de extinção de grupo. Uma razão que CAIA com a dose
> indicaria mecanismo novo — provavelmente retenção (seção 10.4), em que o material fica preso sem o
> líder parar.
>
> **Ressalva de amostragem, pré-registrada:** com quadros a cada 200 iterações (3-4 unidades de `t`)
> e o líder andando 3-4 dx por quadro, a maioria dos contatos a 1 dx **não é amostrada**. Daí os
> n=11. Para R2/R3 o `n` sobe pela frequência de contato, não pela resolução — a leitura continua
> sendo a **forma da varredura**, não o valor num limiar.

### 15.5 CORREÇÃO (2026-10-28) — os números de 15.1–15.4 vêm de uma parede com dois artefatos

O R1 rodou com a parede em **grade** (contorno em degraus, atrito numérico de até 390% na
tangencial, §19.1). O R1′ corrigiu a parede mas trouxe um artefato oposto (o shifting empurrava o
líder contra ela, §20.1). **O R1″ (anel + projeção do shifting, §20.6) é a medida válida.** O que
muda em cada subseção:

| | R1 (grade) | R1′ | **R1″ (válido)** |
|---|---:|---:|---:|
| razão de parada a 1 / 1.5 / 2 dx | 10.8× / 7.8× / 5.3× | 1.6× / 1.2× / 1.1× | **4.6× / 2.7× / 2.2×** |
| razão a 3 dx (ponto de operação) | 3.0× | 2.2× | **1.8×** |
| P(parada \| contato) a 1 dx | 90.9% ±8.7 (n=11) | 20.0% ±10.3 | **53.3% ±12.9 (n=15)** |
| avanço em contato / livre, 1 → 4 dx | 0.12 → 0.90 | 0.24–0.29 | **0.18 → 0.63** |
| deixa de ser líder: contato / livre | 45.5% / 15.8% | 80.0% / 23.7% | **73.3% / 23.2%** |
| \|arco\| azimutal: contato / livre | 0.48 / 0.40 dx | 0.12 / 0.32 | **0.22 / 0.37** |
| arco com sinal (quiralidade) | +0.32 ±0.37 | +0.51 ±0.32 | **+0.35 ±0.24** |
| divisões / fusões / extinções de grupo | 13 / 0 / 10 | 17 / 3 / 19 | **17 / 1 / 19** |

- **15.1 (parada):** continua o desfecho dominante e continua de curto alcance (a razão cresce
  quando o limiar aperta), mas é **cerca de metade** do que o R1 mediu: 4.6× a 1 dx, não 10.8×. A
  outra metade era o atrito da parede em degraus.
- **15.1, "deixa de ser líder":** fica **mais** forte, não menos (73% contra 23%, 3.2×). O líder em
  contato quase sempre perde a liderança, mesmo quando não para por completo: fica para trás porque
  avança a 18–37% do livre nos últimos 2 dx.
- **15.2 (desvio):** continua ausente — o arco em contato fica **abaixo** do livre. O líder não
  contorna o pilar.
- **15.2 (quiralidade):** o arco com sinal ficou +0.32 / +0.51 / +0.35 nas três versões, sempre
  positivo e agora com erro menor que o valor (+0.35 ±0.24, ~1.5σ). **Ainda não é detecção**, mas
  deixou de ser negativo limpo. Vigiar em R2/R3.
- **15.3:** a parede corrigida dá mais extinções de grupo (19 contra 10) e algumas fusões.
- **15.4 (síntese):** "a rugosidade mata o líder que ela encontra" continua — pela perda de
  liderança, não pela parada. O "−3.1%" que 15.4 explicava não é mais o número: ver §20.7.
- **A predição para R2/R3 (razão aproximadamente constante com a dose, frequência de contato
  crescendo) se mantém**, ancorada agora no R1″: razão 1.8× a 3 dx e 4.6× a 1 dx.

---

## 16. MEIO COMO GRAFO — gerador estocástico implementado e caracterizado (2026-10-06)

Sugestão do Prof. Cesar: descrever o obstáculo como **grafo com regiões proibidas**, com a condição
de contorno **descrita estocasticamente**. O grafo é a triangulação de **Delaunay** dos centros
(nó = pilar, aresta = garganta, peso = distância de superfície a superfície) — a abstração padrão de
meio poroso (*pore network*), que vale igual para rede e para campo sorteado e por isso resolve o
que fazia a dose não ser comparável entre os dois.

### 16.1 O que foi implementado — aditivo, e verificado

`PILAR_JITTER` e `PILAR_SEED` em [main.py](../main.py). O jitter desloca cada centro por uma
gaussiana de `JITTER·Λ`, **depois** do clip de domínio e da exclusão, com rejeição de sobreposição /
saída de domínio / entrada na zona de exclusão, re-sorteando só o centro ofensor. Logo **a contagem
de pilares e `φ` ficam IDÊNTICOS ao caso regular** e só a distribuição de gargantas muda.

`pilares.json` passou a gravar `jitter`, `semente` e `gargantas_dx` (percentis, média, desvio, grau
médio) — porque no caso desordenado é a **distribuição de pesos de aresta** que é a dose, e
`lambda` perde sentido.

**Verificado sem rodar solver:** com `PILAR_JITTER = 0` o conjunto de centros é **bit-a-bit idêntico**
ao `pilares.json` arquivado do R1 (92 centros, max\|diferença\| = 0.000e+00), e o `default_rng` do
jitter **não toca o estado global do `np.random`**, então o ruído do inóculo fica intacto.

**O solver não mudou, e não precisa mudar:** `PILAR_LAMBDA` é usado em **uma linha** (só no gerador).
A `ForcaContornoPilar` é par-a-par em `XIJ`/`RIJ`, e as máscaras de métrica, as guardas de deposição,
o `n_pen` e o `void_fraction` leem `_pilar_c` / `_pilar_rp` / `_pilar_pts`. Nada disso sabe que a
rede é regular.

### 16.2 Caracterização do meio — `tools/compara_meio.py`, sem simulação

Λ = 28 dx (dose do R1), `r_max` = 3.852, braço 5.4 dx, semente 20261006:

| jitter | N | φ | garganta p10/p50/p90 (dx) | < braço | percola | contatos/braço | **a₆ do meio** | **fase** |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **0.00** | 92 | 10.7% | 18.00 / 18.00 / 18.00 | 0.0% | 100% | 1.10 (84%) | **0.465** | **30.0°** |
| 0.05 | 92 | 10.7% | 15.53 / 18.32 / 22.59 | 0.0% | 100% | 1.13 (84%) | 0.444 | 29.7° |
| 0.10 | 92 | 10.7% | 13.20 / 18.63 / 26.91 | 0.0% | 100% | 1.17 (83%) | 0.391 | 29.0° |
| 0.20 | 92 | 10.7% | 8.94 / 19.08 / 33.13 | 4.3% | 100% | 1.18 (83%) | 0.208 | 25.8° |
| **0.30** | 92 | 10.7% | 5.92 / 19.28 / 39.49 | 8.5% | 100% | 1.01 (79%) | **0.060** | 8.8° |

`a₆ do meio` = amplitude do modo 6 de `d(θ)`, a distância radial até o **primeiro** pilar. É a
**causa geométrica** do travamento de fase que o R1 mostrou (§12), e tem de desaparecer com desordem
para o travamento ser atribuível à rede.

> ### `PILAR_JITTER = 0.30` é o ponto de operação do controle de anisotropia
>
> Ele derruba a assinatura de 6 dobras do meio **7.8×** (0.465 → 0.060) mantendo **φ, N e contatos
> por braço praticamente fixos** (1.10 → 1.01, interceptação 84% → 79%). Muda a **ORDEM**, não a
> **quantidade** de obstáculo — que é exatamente o que um controle precisa fazer.

### 16.3 A previsão de virada de fase da §14.1 está CONFIRMADA pela geometria

Rodando nas outras duas doses, com jitter 0:

| | Λ | **a₆ do meio** | **fase** | gargantas < braço | percola |
|---|---:|---:|---:|---:|---:|
| R1 | 28 dx | **0.465** | **30.0°** | 0.0% | 100% |
| R2 | 20 dx | 0.118 | **0.0°** | 0.0% | 100% |
| R3 | 16 dx | 0.124 | **0.0°** | 0.0% | 100% |

**Fase 30.0° em R1 e 0.0° em R2/R3** — exatamente a virada que a §14.1 pré-registrou a partir de
qual coroa sobrevive à zona de exclusão, agora obtida por um caminho independente (azimute do
primeiro pilar, não espectro da colônia).

### 16.4 CONFUNDIDOR NOVO — o a₆ do meio é 4× menor em R2 e R3

0.465 (R1) contra 0.118 e 0.124. Causa: `d(θ)` é a distância ao primeiro pilar, e numa rede mais
densa o primeiro pilar está perto em toda direção, então o contraste angular cai.

**Consequência quantitativa.** No R1 a colônia atingiu a₆ = 0.223 com meio a 0.465, ou seja uma
**resposta de 0.48**. Se a resposta for proporcional, R2 e R3 dariam a₆ de colônia ≈ **0.057 e
0.060** — contra o piso do liso, que é 0.022 ± 0.011 com máximo 0.036. **O sinal ficaria quase
indistinguível do ruído, por razão puramente geométrica.**

> **PRÉ-REGISTRO, acrescentado à §14.1:** a anisotropia entre doses é comparada como **razão de
> resposta `a₆(colônia) / a₆(meio)`**, não como `a₆` da colônia. O valor do R1 é **0.48** e é esse
> que as outras doses têm de reproduzir. Reportar `a₆` de colônia cru entre doses está **proibido**,
> pelo mesmo motivo que comparar em `t` casado: as duas versões medem geometria do meio, não resposta
> da colônia.

### 16.5 Percolação NÃO morde na nossa faixa — negativo honesto

Em nenhuma combinação de Λ e jitter o subgrafo de gargantas largas deixa de atravessar: sempre
100%, caindo a 99% só em Λ=16 com jitter 0.30, e isso com **38.8% das gargantas mais estreitas que
o braço**. O critério de percolação que eu havia proposto como "limiar que cai entre R2 e R3"
**não se materializa** — a rede permanece conexa em gargantas largas mesmo muito desordenada.

E a correção que vem com isso: eu havia dito que a garganta do R3 era menor que o braço. **Errado** —
a garganta é `Λ − A` = 6.0 dx > 5.4 dx (1.11×). O que é menor que o braço em R3 é o **corredor
reto**, `0.866Λ − A` = 3.9 dx. São perguntas diferentes: "cabe entre dois pilares?" (sim, apertado)
contra "cabe em linha reta?" (não).

> **E essa distinção rende a melhor previsão para o R3.** O corredor reto (3.9 dx) é menor que o
> braço (5.4 dx), então o braço **não pode seguir reto** — tem de desviar. E a §15 mediu que **o
> líder não desvia, ele para** (arco azimutal muda 20% enquanto o avanço radial cai 8×). **Logo R3
> deve mostrar extinção de líder pesada, não braços sinuosos.** Deriva de duas medições
> independentes e é mais forte que a predição 1 da §3, que o R1 já contrariou.


---

## 17. PISO DO CAMPO DE `c_s` — e porque ele e o HALO, nao o nivel no corpo (2026-10-07)

**Decisao da usuaria.** Qualquer alavanca de morfologia tem de preservar o campo de `c_s` na
zona de expansao. O piso precisava de ancora de literatura, e a escolha da grandeza importa.

### 17.1 O nivel no corpo NAO e ancoravel

O Trinschek (b) da `corpo` = 0.99 do teto, e nos estamos em 0.44 (P2R23). Mas pela licao #99
esse 0.99 e **propriedade do MODELO dele** — Gamma e advectado com o filme, entao satura sobre
toda a pegada por construcao. **Nao e dado experimental.** A medida experimental direta,
[T14] (MSI de metabolitos), diz o oposto: os ramnolipideos ficam *"delineando os tendrils e o
espaco alem deles"*, nao uniformes sobre a colonia.

Portanto a distancia 0.44 x 0.99 **nao mede um defeito** — mede uma diferenca de formulacao.

### 17.2 O que a literatura ancora e o ALCANCE do halo

| | corpo | baia | 0.9R | ponta | **`L`/R** |
|---|---:|---:|---:|---:|---:|
| **Trinschek (b)** | 0.99 | 0.95 | 0.71 | 0.54 | **0.18** |
| P2R21 (sem conduzir) | 0.20 | 0.17 | 0.12 | 0.07 | 0.16 |
| **P2R23 (baseline)** | 0.44 | 0.25 | 0.14 | 0.08 | **0.17** |
| **R1 (Pass L)** | 0.62 | 0.34 | 0.19 | 0.11 | **0.18** |
| P2R26 (beta 6.5) | 0.30 | 0.17 | **0.07** | **0.04** | **0.11** |

Nessa grandeza **ja estamos no alvo**: 0.17-0.18 contra 0.18 do Trinschek.

> ### ~~PISO PRE-REGISTRADO: `L`/R >= 0.17~~ — SUBSTITUÍDO em 2026-10-28 (ver 17.5)
>
> ~~Nenhuma alavanca de morfologia pode levar o halo abaixo disso. Reportar junto `0.9R` e
> `ponta`, que sao o campo onde a Marangoni de fato trabalha.~~

### 17.3 Consequencia imediata: `beta` = 6.5 REPROVADO

O P2R26 entregou dedo +47%, AR 9.7 -> 14.1 e nucleo 0.385 -> 0.280 **sem abrir buraco**
(`void_07` ate caiu, `mean_sig_all` 0.9981). Mas derrubou o halo para **0.11** (-35%), o
`0.9R` de 0.14 para 0.07 e a ponta de 0.08 para 0.04 — **o campo some justamente onde a
Marangoni trabalha**. A colonia ficou 26% maior com a mesma producao e o campo diluiu.

**Revertido para `beta` = 5.0.** Arquivado em `runs/swarm/P2R26_beta65` como medido-e-reprovado.

### 17.4 E isso protege a serie de rugosidade

O R1 **ja rodou sobre o P2R23** e tem o melhor campo da serie (corpo 0.62, halo 0.18). O canal
2 da secao 3 e justamente o quimico, com a sentinela pre-registrada `max c_s` no agar vizinho
em **exatamente 0.300**, encostada no gate — margem zero.

Mudar o motor agora quebraria duas coisas: **R1 deixaria de ser comparavel** (a dose se mede
contra o liso, e o liso mudaria), e **o canal quimico ficaria cego** — medir a perturbacao do
pilar sobre um campo 35% mais curto e medir sinal sobre ruido.

**Enquanto a serie R1-R4 estiver aberta, o motor do liso nao se mexe.**

### 17.5 SUBSTITUIÇÃO DO PISO (2026-10-28): motor + réguas de forma; o halo vira grandeza reportada

**Por que o piso de `L`/R caiu.** Ele tinha uma única âncora numérica, o 0.18 do painel (b) do
Trinschek — e esse painel é resultado de **modelo**, desqualificado como alvo de campo pela mesma
razão do 17.1 (lição #99: Γ advectado com o filme). As âncoras experimentais disponíveis para o
alcance relativo do surfactante vão de **0.051 a 1.519** — uma faixa de 30×, que não sustenta piso
nenhum. Um critério sem âncora é escolha livre (lição #101: régua persistida não é régua validada).

**O que entra no lugar — dois critérios, os dois ancorados:**

| | critério | âncora |
|---|---|---|
| **(i) motor** | `a_mar_bio_p95` (mediana t > 30) **≥ 0.70** | o baseline P2R23 (0.73); é a força que de fato move os líderes, medida nas células vivas (§2.5.1 do CLAUDE.md), e não depende de formulação do campo |
| **(ii) forma** | as cinco réguas das lições #95/#99 **dentro da faixa da literatura**: dedos 9–15, área/disco 0.35–0.39, AR 3.8–11.2, núcleo ~0.36 R, baía ~0.34 R | fotos de PA14 e [T12] (morfometria experimental), medidas com o mesmo algoritmo dos runs |

**O halo (`L`/R, `0.9R`, `ponta`) continua sendo reportado** em toda rodada, sem piso: ele
descreve o campo, mas não decide aprovação.

**Escopo:** (i) e (ii) são **guardas para alavancas do modelo liso** (motor, campo, forma). Na série
de rugosidade o pilar é a perturbação que se está medindo; lá `a_mar_bio_p95` é **resposta** (o
canal químico da §3, item 2), não guarda. Ex.: o R1″ tem `a_mar_bio_p95` = 0.63 (−13% contra o
liso), e isso é resultado da série, não reprovação.

**O P2R26 (`beta` 6.5) continua REPROVADO pelos critérios novos** — por razões diferentes das
originais, e mais fortes, porque cada uma tem âncora:

| | P2R26 | faixa da literatura | |
|---|---:|---:|---|
| (i) `a_mar_bio_p95` | **0.587** | ≥ 0.70 | ❌ |
| dedos | **17** | 9–15 | ❌ |
| área / disco | **0.32** | 0.35–0.39 | ❌ |
| AR | **13.7** | 3.8–11.2 | ❌ |
| núcleo / R | **0.28** | ~0.36 | ❌ |
| baía / R | **0.25** | ~0.34 | ❌ |

Seis de seis fora. A reprovação original (halo −35%) segue registrada no 17.3 como descrição, não
como critério. O mecanismo por trás está na lição #109 do CLAUDE.md (*`beta` é auto-limitante*).

---

## 18. SOMBRA QUÍMICA EM VOLTA DOS PILARES — medida no R1 (2026-10-26)

A verificação que o plano §1 previa ("`c_s` e `c_n` num anel de 2 dx em volta dos pilares contra o
mesmo raio longe deles") e que nunca tinha sido feita. Responde ao canal químico da seção 3
(predição 2). `tools/sombra_quimica.py runs/rugosidade/R1_lambda28 runs/swarm/P2R23_conduz`,
pós-processamento puro.

### 18.1 Régua — e o que foi decidido antes e depois de olhar

Fixado **antes**: anel = fluido a (0, 2 dx] da superfície; montante/jusante = cos ≤ −0.5 / ≥ +0.5
em relação à direção radial do pilar; assimetria `S = ln(média montante / média jusante)`;
**pareamento** — o mesmo `S` no liso, nas mesmas posições, no quadro de `t` mais próximo, e o efeito
é `dS = S_R1 − S_liso` (o gradiente radial natural cancela); classes `dentro` (centro em r < 0.9 R99)
e `frente` (R99 ≤ r < R99 + 1.03 — o halo, onde o bloqueio ANTES do toque apareceria).

Acrescentado **depois** do primeiro resultado, e por isso declarado: um controle de **pilares
falsos** — as mesmas posições giradas 30° em torno da origem, mantidas só as cujo anel fica a mais
de 2 dx da superfície de qualquer pilar real (o primeiro filtro, 4 dx, não deixava nenhum falso
nas faixas de raio). Motivo: o R1 tem mais `c_s` que o liso em toda parte (seção 10.3), e um campo
apenas mais íngreme daria `dS > 0` sem pilar.

### 18.2 Na FRENTE (halo, antes do toque): o pilar deixa sombra, com a magnitude da teoria

| t | `dS` real | `dS` falso (controle) |
|---:|---:|---:|
| 4.2 | +0.861 | +0.035 |
| 8.2 | +0.588 | −0.004 |
| 12.4 | +0.484 | +0.022 |
| 16.6 | +0.594 | +0.152 |
| **janela [35,50]** | **+0.417 ± 0.028** (montante/jusante ×1.52) | — (nenhum falso nessa faixa de raio) |

Real 0.48–0.86 contra controle −0.004–0.15 nos quatro quadros em que os dois existem. O pilar
**enriquece o lado de dentro e sombreia o de fora** no halo, antes de a colônia chegar a ele.

**A magnitude é a da solução clássica**, não um número a explicar. Para um cilindro isolante num
gradiente uniforme, em regime permanente, a perturbação do campo é dipolar e multiplica a
assimetria local por `1 + R²/r²`. Com `R` = 5 dx e o anel em ~6 dx: **1.69** previsto. Medido na
janela: `(S_liso + dS) / S_liso = (0.631 + 0.417) / 0.631` = **1.66**. É ao mesmo tempo resultado
e **verificação da implementação** do fluxo zero (contorno adiabático de Cleary & Monaghan 1999):
omitir os pilares das somas de Brookshaw reproduz a física do contorno isolante.

**O nível não muda:** `N` (média do anel R1/liso) é +0.12 nos reais e +0.10 a +0.24 nos falsos — a
diferença de nível é global entre os runs (seção 10.3), não local ao pilar. Coerente com a
perturbação dipolar, que tem média zero em volta do cilindro.

**`c_n` na frente: nenhum efeito** (`dS` = −0.001 ± 0.003). O nutriente à frente é ~1 e uniforme; sem
gradiente não há o que sombrear.

### 18.3 DENTRO da colônia: inconclusivo

| janela [35,50] | real | falso |
|---|---:|---:|
| `c_s` `dS` | +0.108 ± 0.031 | **−0.171 ± 0.040** |
| `c_n` `dS` | −0.106 ± 0.070 | +0.086 ± 0.032 |

O controle **não é plano** dentro da colônia: o `dS` dos falsos deriva de +0.07 (t=20) a −0.22 (t=50),
porque a colônia do R1 é outra realização, com braços 35% mais largos (seção 10.4). O sinal real
(0.11) fica dentro da variação do controle, e os falsos são só 6 numa faixa de raio. **Não atribuível.**

### 18.4 O que isto corrige e o que deixa em aberto

**Corrige a estimativa ex-ante do plano §1** ("diferença montante/jusante de `c_s` < 30%"). Ela
raciocinava por `R_p/L`, mas o fator relevante é `1 + R²/r²`, que **não depende de `L`**: qualquer
obstáculo isolante amplifica ~1.7× a assimetria local no seu anel. No halo isso dá ×1.52 em
montante/jusante — acima dos 30%.

**Responde à pergunta 2 da série pela metade.** O bloqueio químico **existe como perturbação do
campo antes do toque** — mas é **local** (dipolo decai como `R²/r²`, alcance de ~1 raio de pilar) e
não há, no R1, efeito dele sobre a expansão (`dR/dt` −3.1%, seção 10.2).

**Hipótese, NÃO medida:** o acúmulo a montante cria entre o líder que chega e o pilar um `c_s`
crescente em direção ao pilar, ou seja `−∇c_s` apontando **para longe do pilar** — uma repulsão
química que freia o líder antes do contato mecânico. Se existir, contribui para a PARADA da seção 15
e para a queda de `a_mar_bio_p95` da seção 10.4, que nunca foi separada por medição. Teste:
componente radial de `au_mar` dos líderes em função da distância ao pilar, real × falso.

**Para R2/R3:** a amplificação por pilar é a mesma (`R` fixo em todas as doses); o que muda é quantos
pilares o halo atravessa. Rodar esta ferramenta nas duas doses, com o mesmo controle.

### 18.5 A REPULSÃO QUÍMICA — hipótese da 18.4 testada e CONFIRMADA (2026-10-26)

`tools/repulsao_quimica.py runs/rugosidade/R1_lambda28 runs/swarm/P2R23_conduz`. Régua fixada antes
de rodar (docstring): líderes como em `desfecho_contato.py`; só o pilar mais próximo **à frente**
do líder; `cos_g` = cosseno entre `−∇c_s` (direção do flagelo) e a direção líder→pilar; `a_p` =
Marangoni projetada nessa direção; **pareado** com os líderes do liso contra as mesmas posições.

**As direções foram reconstruídas.** O HDF5 grava só as magnitudes (`au_mar`, `au_flag`); `ax_mar`,
`grad_cs_x` etc. estão vazios. A reconstrução usa as fórmulas do solver (fontes = fluido, pilar
fora; filler fora como fonte de `c_s`; KGC com todo o fluido; gate por `|∇ρ_b|`) e foi validada:
`|a_mar|` reconstruída × `au_mar` gravada, **correlação 0.9999 (R1) e 1.0000 (liso)**, razão
mediana 0.9999.

| d à superfície | n R1 | `cos_g` R1 | n liso | `cos_g` liso | diferença | `a_p` R1 | `a_p` liso |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0-1 dx | 4 | **−0.195** | 14 | +0.820 | **−1.015** | **−0.398** | +2.803 |
| 1-2 dx | 14 | **−0.022** | 6 | +0.938 | **−0.960** | **−0.265** | +3.556 |
| 2-3 dx | 18 | +0.135 | 13 | +0.772 | −0.637 | +0.281 | +1.425 |
| 3-5 dx | 37 | +0.415 | 18 | +0.946 | −0.531 | +1.001 | +2.049 |
| 5-7 dx | 33 | +0.771 | 28 | +0.542 | +0.229 | +1.558 | +1.472 |
| 7-10 dx | 74 | +0.726 | 60 | +0.568 | +0.158 | +1.012 | +0.861 |
| 10-15 dx | 34 | +0.636 | 7 | +0.494 | +0.141 | +0.631 | +0.758 |

Erro da mediana de `cos_g` do R1: 0.03-0.11 por faixa.

**A forma por distância é a prova (lição #68), e ela está lá:** a diferença vale −1.0 a 1-2 dx,
encolhe para −0.5 a 3-5 dx e **troca de sinal** a partir de 5 dx, assentando num resíduo pequeno e
constante (+0.14 a +0.23) — que é a diferença global entre as duas realizações, não o pilar. O
critério pré-registrado ("R1 abaixo do liso em d < 3, diferença sumindo em d ≥ 7") passa, com a
ressalva de que o resíduo distante não é zero.

**O que isto quer dizer:**
1. No liso, o motor de um líder que vai em direção àquele ponto aponta PARA ele (`cos_g` +0.8 a
   +0.95). No R1, a menos de ~2 dx do pilar, o motor aponta **de lado ou para trás** (`cos_g` ≈ 0 a
   −0.2) e **a Marangoni passa a empurrar para LONGE do pilar** (`a_p` −0.27 a −0.40, contra +2.8 a
   +3.6 no liso).
2. **O alcance químico (~5 dx) é cinco vezes o mecânico** (o kernel de contato zera a 1 dx). O
   líder começa a ser freado pelo campo antes de qualquer força de contato existir.
3. A PARADA da seção 15 é **químico-mecânica**, não só mecânica: P(parada | contato) a 1 dx (10.8×)
   cai exatamente na faixa em que o motor se inverte. E o "não desvia" da 15.2 fica coerente: a
   repulsão tem o pilar à frente, então ela aponta majoritariamente **para trás**, freando, e não
   para o lado, desviando.
4. Mecanismo é o da 18.2: o pilar isolante acumula `c_s` no lado de montante; entre o líder e o
   pilar, `c_s` passa a crescer em direção ao pilar.

**Ressalva de operador, declarada.** Perto do pilar o suporte do kernel é truncado (o pilar não é
fonte). No gradiente SEM correção, usado pelo flagelo, parte da rotação de `cos_g` pode ser
**viés de truncamento**, não campo. A Marangoni usa KGC, que restaura a consistência de 1ª ordem
sob suporte truncado (Bonet & Lok 1999; teste-ouro do §12 do CLAUDE.md) — logo a inversão de
`a_p`, que sai da Marangoni corrigida, é do **campo**. Exceção: onde `det(M) < 0.25` a KGC cai para
identidade; a fração de líderes nesse fallback a < 2 dx não foi medida. Amostragem: n = 4 a 0-1 dx;
a conclusão se apoia nas faixas 1-5 dx (n = 14-37), não na primeira.

**Consequência para a pergunta 2 da série:** o bloqueio químico não aparece na expansão média do
R1 (−3.1%), mas **existe e age sobre o líder**, antes do toque e com alcance maior que o do contato.
Em φ = 11.6% ele é raro porque o contato é raro (seção 15.4). Para R2/R3: a repulsão por pilar
deve ser a mesma (`R_p` fixo); o que cresce é quantos líderes entram no raio de ~5 dx.

---

## 19. CORREÇÃO DO CONTORNO DE CONTATO — análise preditiva, implementação e pré-registro do R1′ (2026-10-27)

### 19.1 O defeito, contra a referência

Monaghan & Kajtar (2009) provam, pela fórmula de soma de Poisson (eq. 2.4–2.6), que forças radiais
somadas sobre partículas de contorno imitam uma força **normal e uniforme** quando o espaçamento da
parede é menor que o comprimento do kernel (`Δ < h`); recomendam β = 3 para kernels cúbicos (§3,
Tabelas 1–8: `tang/perp` e `dperp/perp` < 10⁻⁵). A implementação do R1 tem **Δ = dx = 2·h_w**, com
o contorno feito do disco de pontos da grade (escada), e não de partículas sobre a curva.

Medido com a fórmula do código (repetição dos testes da §3 do artigo):

| | var. da normal | tang/normal |
|---|---:|---:|
| parede reta, R1, a 0.5 dx | 136% | 47% |
| pilar real do R1, na distância de parada (0.75–0.9 dx) | normal vai de 0 ao máximo | **75% a 390%** |
| artigo (β = 3), a 1 dp | ~10⁻⁵ | ~10⁻⁵ |

E medido no R1: a colônia em contato para a **0.82–0.90 dx** da partícula de pilar mais próxima,
a 5.32 dx do centro na mediana e **4.79 dx no p1** — ela entra nos degraus da escada.

**A parede do R1 é impenetrável (`n_pen` = 0), mas é um cilindro rugoso na escala de 1 dx**, com
empurrões tangenciais que podem agir como atrito numérico. Esse atrito **não foi separado** da
PARADA medida na seção 15.

### 19.2 Restrições e alternativas (análise preditiva offline, sem solver)

| | restrição | fonte |
|---|---|---|
| R-a | Δ < h_w, partículas sobre a curva, força ÷ β | artigo, eq. (2.6), (4.1), §3 |
| R-b | força **zero** no ágar em t = 0 | lição #40 |
| R-c | mesmo raio efetivo (mesma dose) | seção 1 |
| R-d | σ_a e densidade perto do pilar | Violeau §3.6; C2 do §2.5 |
| R-e | a colônia não para em cima de partícula que não a repele | — |

| alternativa | R-a | R-b | R-c | R-d | R-e | veredito |
|---|---|---|---|---|---|---|
| C0 (R1): grade, β=1, h_w=0.5 dx | ✗ | ✓ | ✓ | ✓ | ✓ | não conforme |
| C1: grade com h_w = 1 dx (Δ = h, artigo literal) | ✗ 49–94% de tangencial (escada) | ✗ **48 partículas de ágar/pilar empurradas em t=0** (até 5.7) | — | ✓ | ✓ | rejeitada |
| C2 com vão de 1 dx (anel 4.4, interior ≤ 4.4) | ✓ | ✓ | ✓ | ✗ σ_a 0.71 | ✓ | rejeitada |
| C2 com interior ≤ 5.1 dx | ✓ | ✓ | ✓ | ✓ 0.91 | ✗ colônia a 0.05–0.13 dx de partícula de densidade | rejeitada |
| **C2 adotada** | ✓ | ✓ | ✓ | ~ | ✓ | **adotada** |

**C2 adotada:** anel de contato sobre a curva em **4.6 dx**, espaçamento **dx/4 (β = 4)**,
**h_w = 0.4 dx**, força **÷ β**, `d` = 0; interior da grade **≤ 4.9 dx** só para densidade,
viscosidade e *shifting*; fluido **> 5.4 dx**.

- β = 4 (Δ/h_w = 0.62), não 3 (0.83): a 0.45 dx, var 0.03% / tang 0.06% contra 2.5% / 0.9%.
- `d` = 0 basta: no meio do vão, a 0.05 dx de profundidade, a força normal é 117 (C0: 5.7).
- K ÷ β mantém a força **média** por distância igual à do C0 (é o propósito do fator na eq. 4.1),
  logo a calibração de K (plano §4) continua válida: muda a uniformidade, não a intensidade.
- Raio de parada previsto 5.05–5.23 dx (C0: ~5.1–5.3); φ efetiva cai até ~7%. Anel em 4.6 é o
  mais externo compatível com R-b.
- **Compromisso aceito (R-d):** σ_a do ágar **parado** em volta do pilar cai a 0.75. Ele tem
  pressão nula (não se move) e nunca dispara inserção (gate `ρ_b` > 0.5). A colônia, quando chega,
  ocupa esse vão e fica a 0.58–0.76 dx da partícula de densidade mais próxima (C0: 0.82–0.90).

**Conformidade:** partículas sobre a curva ✓; Δ < h ✓; ÷ β ✓; `d` = 0 ✓; kernel cúbico com a
normalização absorvida em K ✓; K por equilíbrio com a força que empurra o fluido contra a parede ✓
(eq. 4.7, com a força motriz no lugar da hidrostática). **Desvios declarados:** h_w = 0.4 dx em vez
de h ~ dp (lição #40); parede na densidade por interior da grade (o artigo usa continuidade; nós,
`SummationDensity` — Adami 2012); não-deslizamento pela `ViscousForce`, não por Π_aj.

### 19.3 Implementação

- `main.py`: `PILAR_CONTORNO` (`"grade"` reproduz o R1; `"anel"` é o padrão novo), `PILAR_BETA` = 4,
  `PILAR_ANEL_R` = 4.6, `PILAR_DENS_R` = 4.9, `PILAR_FLUIDO_R` = 5.4 (em dx),
  `PILAR_HW` = 0.4 (anel) / 0.5 (grade). `_faz_pilares` cria o array novo `pilar_sup` (o anel) e
  devolve `[pilar, pilar_sup]`; grava em `pilares.json` `contorno`, `hw_dx`, `beta`, `anel_raio`,
  `n_anel_por_pilar`, `dens_raio`, `fluido_raio`. `"raio"` continua o nominal (5 dx): máscaras,
  `n_pen` e ferramentas não mudam.
- `src/scheme.py`: `pilares_contato` (padrão = `pilares`); a `ForcaContornoPilar` usa
  `["pilar_sup"]` no modo anel e recebe K ÷ β.
- Guardas de deposição e `n_contato`: a árvore de pilar inclui interior + anel.
- **`n_pen` mantém a definição pré-registrada** (centro < R_p − 0.5 dx = 4.5 dx, abaixo do anel).

### 19.4 Verificação offline — FEITA (construção das partículas, sem solver)

| verificação | resultado |
|---|---|
| modo `"grade"` × quadro inicial do R1 | fluido 61 127 e pilar 6 994 partículas, **posições bit-a-bit iguais** |
| R-b: menor distância fluido–anel em t=0 | **0.807 dx** ≥ alcance 0.80 dx → força **zero** (margem de 0.007 dx: os centros caem em posições arbitrárias da grade) |
| R-a sobre os 92 pilares construídos | pior pilar: 0.027% / 0.059% a 0.45 dx; 0.84% / 0.22% a 0.63 dx |
| R-d: σ_a do ágar a < 1.5 dx da borda | min 0.754; 39% abaixo de 0.85 (previsto, aceito) |
| partículas | fluido 59 515 (−1 612 em relação à grade), densidade 6 774, anel 10 672 (116/pilar) |

**Rodadas curtas — FEITAS (2026-10-27, autorizadas; 1 thread como o R1, `OMP_NUM_THREADS`=10 sem `--openmp`):**

| teste | resultado |
|---|---|
| **V0**: código atual × `HEAD` (e8034d6), ambos com `use_pilares = False`, até a iteração 600 | as **81 arrays do HDF5 idênticas** nas iterações 0/200/400/600 (`t` = 4.1658 / 7.3000 / 11.2183 iguais) e **`log.csv` idêntico** — o modo anel e o `RASTRO_W_RAZAO` não tocam o baseline liso |
| **modo anel, 200 iterações** (t = 4.22) | compila e roda com `pilar_sup` (79 069 partículas); `ax_rep`/`ay_rep` = **0 exato** em todo o fluido nos quadros 0 e 200; `n_pen` = 0 (menor distância fluido–centro 5.41 dx); anel parado; ágar a < 1 dx da superfície deslocou no máximo 1.3e-7 dx — **o aborto extra da 19.5 não disparou** |

**Limite do teste:** em t = 4.2 a colônia tem `R99` = 0.75 e o pilar mais próximo está a 1.24 da
origem — **nenhum contato aconteceu**. O comportamento da parede lisa sob contato só é medido no R1′.

### 19.5 PRÉ-REGISTRO — R1′ (Λ = 28 dx, modo anel; mesma semente 20260806, 10 threads, t = 50)

O liso (P2R23) não muda. A rede, φ nominal e a força média também não; muda a uniformidade da parede.

| grandeza | predição | leitura |
|---|---|---|
| `dR/dt` t∈[25,50], `R99`(t=50) | dentro do ruído do R1 (0.0733; 3.852) | a parede não mudou de intensidade média |
| fase do m=6 no pico | 27–30° | a rede não mudou |
| **P(parada \| contato) a 3 dx (R1: 24.6% ±5.5, razão 3.0×)** | **razão mantida ≥ ~2×** se a parada é químico-mecânica (§18.5) | **se cair para ~1×, parte da parada do R1 era atrito numérico da parede rugosa** e a §15 tem de ser corrigida |
| varredura da razão por limiar | ainda cresce quando o limiar aperta | forma da curva é a prova (#68) |
| repulsão química (§18.5): `cos_g` a 1–3 dx | negativo/≈0 como no R1 | é do campo, não da parede |
| raio de parada (centro) em contato | 5.05–5.23 dx, sem valores < 4.9 | some a entrada nos degraus |
| `n_pen` | 0 em todos os quadros | aborto 2 |

Critérios de aborto: os da seção 6, mais **`a_rep` > 0 no ágar fora do alcance da colônia nas
primeiras 200 iterações** (R-b violada na prática).

**Consequência para a série:** o R1 deixa de ser o ponto da dose mais baixa; R1′ o substitui, e
R2/R3 rodam no modo anel. O R1 fica arquivado como medição da parede rugosa.

### 19.6 RESULTADO — R1′ (`runs/rugosidade/R1p_lambda28_anel`, 2026-10-27; 1 thread como o R1, 1784 s)

| grandeza | predição (19.5) | R1 | **R1′** | liso (P2R23) | |
|---|---|---:|---:|---:|---|
| `R99`(t=50) (quadros, colônia = quorum ou filler) | ≈ R1 | 3.852 | **3.821** | 4.067 | ✅ |
| `dR/dt` t∈[25,50] | ≈ R1 | 0.0733 | **0.0751** | 0.0756 | ✅ (ver nota) |
| m=6: pico de a₆ / fase no pico | fase 27–30° | 0.223 (t=25.6) / 29° | **0.258 (t=27.8) / 28.7°** | — | ✅ |
| `n_pen` | 0 | 0 | **0** | — | ✅ |
| menor distância fluido móvel–centro (t 39→50) | sem < 4.9 dx | 4.77 dx | **5.06 dx** | — | ✅ some a entrada nos degraus |
| repulsão química `cos_g` a 0–3 dx | ≈ 0 / negativo | ≈ 0 / negativo | **−0.04 / −0.04 / +0.03** (liso +0.77 a +0.94) | | ✅ do campo, não da parede |
| **P(parada \| contato) a 3 dx** | razão ≥ ~2× | 24.6% / 8.2% = **3.0×** | 23.7% ±4.4 / 10.8% = **2.2×** | | ✅ no limite |
| **varredura: a razão cresce quando o limiar aperta** | sim | 10.8× / 7.8× / 5.3× (1 / 1.5 / 2 dx) | **1.6× / 1.2× / 1.1×** | | ❌ |

**Leitura.** A previsão que separava as duas hipóteses era a forma da varredura, e ela mudou. No R1 a
parada se concentrava a ≤ 2 dx da superfície (91% de parada a 1 dx); com a parede lisa, a ≤ 2 dx a
taxa de parada (15–20%) fica dentro do erro do controle (12.8%). **A parada dura de curto alcance do
R1 era atrito numérico da parede rugosa** (§19.1). O que sobrevive tem alcance de 3–4 dx — razão
2.2× a 3 e 4 dx, 1.3× a 7–10 dx —, o mesmo alcance da repulsão química (§18.5, reproduzida acima).

**O que o contato ainda faz:** o avanço radial do líder em contato fica em **0.24–0.29× o do líder
livre em TODOS os limiares de 1 a 4 dx** (R1: 0.12 a 1 dx subindo a 0.71 a 3 dx) — frenagem forte,
plana na distância, sem parada completa. A mediana em contato (0.58 dx) cai em cima do limiar de
PARADA (0.25 × 2.44 = 0.61 dx): **a fração de "parada" é sensível a esse limiar e não deve ser lida
sozinha**; o avanço relativo é a medida robusta. E o líder em contato **deixa de ser líder** em
80% ±10 das observações, contra 23.7% ±1.9 livre (3.4×; R1: 2.9×) — a perda de dedo persiste,
pela frenagem e não pela parada.

**A §15 tem de ser corrigida:** onde ela atribui ao contato mecânico uma parada de curto alcance
(10.8× a 1 dx), o mecanismo medido com a parede conforme é **frenagem de alcance químico**.

**Nota sobre `dR/dt`:** os três valores ficam dentro de 3% um do outro; a forma de medir (ajuste
linear de `R99` dos quadros em t∈[25,50]) não separa R1 de liso nesta dose. O efeito no raio
aparece no valor final (−6%), não na taxa da janela.

**DEFEITO NOVO — partícula pinada derivando para o anel por shifting (não dispara aborto, mas
precisa ser corrigido antes do R2).** `a_rep_max` sobe nos três últimos quadros (5.6 → 14.1 →
17.3). Não é contato dinâmico: as duas partículas com a força alta estão **pinadas** (`|u|` = 0;
um filler da inserção e uma célula com `c_n` < 0.6) e se aproximam do anel 0.63 → 0.44 dx entre
t = 42.8 e 49.8. Como o pin zera a velocidade a cada passo, a força de contato não age sobre elas;
o único termo que as move é o `ParticleShift` (aplicado depois do pin). No modo anel, a grade de
densidade termina em 4.9 dx e o fluido começa em 5.4 dx — o vão de 0.5 dx baixa a concentração
`C` junto ao pilar e `−∇C` aponta para dentro. Contraprova: as 3 pinadas a < 6.5 dx de um centro se
**aproximam** 0.33–0.55 dx em t 39→50 no R1′ e se **afastam** 0.16–0.56 dx no R1 (grade sem vão).
São poucas (3) e lentas (~0.05 dx por unidade de t), por isso o R1′ fecha em t = 50 sem `n_pen`;
mas com mais pilares (R2, R3) e sem nada que as segure, chegam ao anel. Correções candidatas (não
implementadas — exigem análise preditiva e aprovação): (a) incluir o anel nas fontes do
`ParticleShift` (a parede entra na concentração `C`, pela mesma razão que a grade do pilar já entra — comentário do `ParticleShift` em `src/scheme.py`); (b) não deslocar
partícula pinada a menos de ~1.5 dx do anel. `a_rep_max` dos três últimos quadros do R1′ é efeito
desse defeito e não entra em nenhuma leitura.

## 20. CORREÇÃO DO SHIFTING NA PAREDE — análise preditiva, implementação e pré-registro (2026-10-28)

### 20.1 O defeito, medido de novo (reconstrução offline do `ParticleShift` em quadros do R1′)

Fontes como no solver (fluido + grade de densidade do pilar), kernel CubicSpline, h mediano, limite
`SHIFT_CAP`·h; população = partículas no gate do shifting (`rho_b` ∈ [0.1, 0.8)) a < 2h do anel.

| t | n | pinadas | a < 1 dx do anel | componente normal / limite (p50) | fração para dentro | idem a < 1 dx |
|---|---:|---:|---:|---:|---:|---:|
| 27.8 | 498 | 1 | 56 | −0.37 | 62% | −0.62 |
| 42.8 | 1244 | 21 | 204 | −0.31 | 61% | −0.65 |
| 49.8 | 1639 | **53** | 303 | −0.19 | 58% | −0.73 |

**O número de pinadas expostas cresce com a colônia** (1 → 21 → 53): o §19.6 contou só as que
estavam perto nos dois quadros comparados. E **as móveis também são empurradas para dentro**: o
empurrão equivale a ~0.0007 dx por passo, ou seja, a uma velocidade de ~0.0023 (dt ≈ 0.016); o
arrasto γ = 60 vezes essa velocidade dá ~0.14, que é a `a_rep_med` medida no R1′ (0.10–0.16).
**A força de contato mediana do R1′ estava segurando o shifting, não o motor.**

### 20.2 Alternativas

| opção | efeito (offline, mesmos quadros) | veredito |
|---|---|---|
| (a) anel nas fontes de `C` | normal / limite = **+1.00** em 93–96% das partículas a < 2h: empurra para fora no limite, sempre — repulsão artificial que afasta a colônia da parede | rejeitada |
| (b) não deslocar partícula pinada | resolve a deriva das pinadas, mas mantém a carga de ~0.14 nas móveis | parcial |
| **(c) anular só a componente normal para dentro** | 0% para dentro; deslocamento cai a 0.89–0.92 do limite; componente tangencial preservada | **adotada** |
| (d) estender a grade de densidade até o fluido | conta densidade em dobro sob o fluido que encosta (motivo da R-d, §19.2) | rejeitada |

(c) é o tratamento de superfície de Lind et al. (2012) [T8] — suprimir a componente normal do
deslocamento junto ao contorno — e equivale à "rejeição de movimento para dentro do sólido" do
modelo de Cesgrav (2026, §2.7).

### 20.3 Implementação

- `src/equations.py`, `ProjecaoShiftParede` (fontes: `pilar_sup`): normal n = Σ_j W_ij (x_i − x_j)
  normalizada, nula além de 2h do anel; se `shift`·n < 0, subtrai a componente normal. Roda no mesmo
  `Group` do `ParticleShift`, depois dele.
- `src/scheme.py`: parâmetro `shift_parede`; a equação só entra com ele ligado **e** com fontes de contato.
- `main.py`: `PILAR_SHIFT_PAREDE = True`; ligado só com `use_pilares` e `PILAR_CONTORNO == "anel"`.
  Propriedades `wall_nx`, `wall_ny` só com `use_pilares`.
- **Liso e modo `"grade"` intocados por construção** (equação ausente, propriedades ausentes no liso).

### 20.4 PRÉ-REGISTRO — testes curtos (exigem autorização)

1. **V0′**: `use_pilares = False`, até a iteração 600, contra o V0 do §19.4 (`HEAD` e8034d6): 81 arrays e `log.csv` idênticos.
2. **200 iterações no modo anel**: compila com a equação nova; `a_rep` = 0 no ágar; `n_pen` = 0.

### 20.5 PRÉ-REGISTRO — R1″ (Λ = 28 dx, modo anel + projeção; mesma semente, mesmo modo de execução do R1′, t = 50)

O R1′ roda com o defeito e deixa de ser o ponto da série: **R2/R3 rodam com a projeção, então a
dose mais baixa tem de ser refeita com o mesmo código** (pareamento). O R1′ fica arquivado.

| grandeza | R1′ | predição R1″ | leitura se falhar |
|---|---:|---|---|
| pinadas a < 2h do anel: aproximação do centro, t 39→50 | +0.33 a +0.55 dx | **≤ 0.05 dx** | a projeção não pega a deriva |
| `a_rep_med` (t > 30, mediana) | 0.16 | **≤ 0.05** (sem a carga do shifting) | a carga não vinha do shifting (§20.1 errada) |
| `a_rep_max` nos três últimos quadros | 5.6 / 14.1 / 17.3 | **< 3**, sem crescimento | idem |
| avanço do líder em contato / livre, 1–4 dx | 0.24–0.29 | **0.24–0.32** (o empurrão é 3–5% da velocidade do líder) | a frenagem dependia do shifting |
| razão P(parada \| contato) a 3 dx | 2.2× | 1.7–2.7× | idem |
| `R99`(t=50), `dR/dt`[25,50] | 3.821, 0.0751 | dentro de ±3% do R1′ | — |
| `n_pen` | 0 | 0 | aborto |

### 20.6 RESULTADO — R1″ (`runs/rugosidade/R1pp_lambda28_anel_shift`, 2026-10-28; 1 thread, 1295 s)

| grandeza | predição (20.5) | R1 (grade) | R1′ (anel, shifting para dentro) | **R1″ (anel + projeção)** | liso | |
|---|---|---:|---:|---:|---:|---|
| pinadas a < 2h do anel: aproximação, t ≈ 39→50 (n = 8) | ≤ 0.05 dx | — | máx +0.55 dx | **máx 0.000 dx** | — | ✅ |
| `a_rep_med` (t > 30) | ≤ 0.05 | — | 0.162 | **0.038** | — | ✅ |
| `a_rep_max`, três últimos quadros | < 3 | — | 5.6 / 14.1 / 17.3 | **0.36 / 1.48 / 0.65** | — | ✅ |
| `n_pen` | 0 | 0 | 0 | **0** | — | ✅ |
| m=6: pico / fase | fase 27–30° | 0.223 / 29° | 0.258 / 28.7° | **0.254 (t=23.2) / 27.9°** | — | ✅ |
| repulsão química `cos_g` a 0–3 dx | ≈ 0 | ≈ 0 | ≈ 0 | **+0.07 / −0.14 / +0.15** (liso +0.77 a +0.94) | | ✅ |
| razão P(parada \| contato) a 3 dx | 1.7–2.7× | 3.0× | 2.2× | **1.8×** | | ✅ |
| `R99`(t=50) | ±3% do R1′ | 3.852 | 3.821 | **3.682 (−3.6%)** | 4.067 | ❌ no limite |
| `dR/dt` t∈[25,50] | ±3% do R1′ | 0.0733 | 0.0751 | **0.0698 (−7%)** | 0.0756 | ❌ |
| avanço em contato / livre, 1 → 4 dx | 0.24–0.32, plano | 0.12 → 0.90 | 0.24–0.29 (plano) | **0.18 / 0.27 / 0.37 / 0.56 / 0.63** | | ❌ |
| razão P(parada \| contato) a 1 / 1.5 / 2 dx | — | 10.8× / 7.8× / 5.3× | 1.6× / 1.2× / 1.1× | **4.6× / 2.7× / 2.2×** | | |

**A correção funciona:** as pinadas param de derivar, e a força de contato mediana cai a 0.038 — a
carga do shifting (§20.1) sumiu.

**E ela muda a leitura do §19.6.** Eu havia previsto que tirar o empurrão do shifting não mexeria na
frenagem (o empurrão é 3–5% da velocidade do líder). Errado: com a projeção, a frenagem **volta a
depender da distância** (0.18 a 1 dx subindo a 0.63 a 4 dx) e a parada de curto alcance **reaparece**
(4.6× a 1 dx). O shifting do R1′ empurrava o líder **em direção à parede e ao longo dela** justamente
no último 1–2 dx, onde a parada se decide; isso mascarava a parada. Ou seja, o R1′ errou para o lado
oposto ao do R1:

| | parede | shifting | razão a 1 dx |
|---|---|---|---:|
| R1 | rugosa (atrito numérico) | sem vão | 10.8× (inflada) |
| R1′ | lisa | para dentro (empurra o líder contra a parede) | 1.6× (mascarada) |
| **R1″** | **lisa** | **sem componente para dentro** | **4.6×** |

**Leitura corrigida (substitui a do §19.6):** existe parada mecânica de curto alcance no contato,
cerca de **metade** da que o R1 mediu (4.6× contra 10.8× a 1 dx); a outra metade era a parede rugosa.
Sobre ela se soma a frenagem de alcance químico (razão ~1.8× até 5 dx, repulsão reproduzida). A 1 dx
são só 15 observações (±13 pontos), então o valor 4.6× tem incerteza grande — a FORMA da varredura
(cresce quando o limiar aperta) é a evidência robusta (#68).

~~**Efeito na expansão é maior que no R1′:** contra o liso, `R99` −9.5% e `dR/dt` −8% (R1′: −6% e −1%).~~
**RETIRADO no §20.7:** medido num instante só (t = 50), onde o liso salta +6%. Na janela, R1′ e R1″
dão o mesmo raio (−6.0% e −5.7%), e o `dR/dt` do R1″ (−4.4%) fica dentro do ruído (5.8%).

**Consequências:** (1) o R1″ é o ponto de dose mais baixa da série; R1 e R1′ ficam arquivados como
medições de artefatos de contorno. (2) A §15 tem de ser corrigida com estes números (não com os do
R1′). (3) R2 e R3 rodam com o código do R1″.

### 20.7 Ruído entre realizações e morfologia do R1″ (2026-10-28, pós-processamento)

**Não existe segunda realização do P2R23.** A `P2R23_t100_pareamento` é idêntica a ela até t = 48
(mesmo log, mesmos raios): mesmo código e mesmo modo de execução, logo é a mesma realização. O piso
disponível continua sendo o par do P2 (§4): `P2_fillerdonor` (10 threads) × `P2_t100` (5 threads).

**`R99` num instante só é instável.** O liso salta de 3.84 (t = 48.3) para 4.07 (t = 50): +6% em 1.7
unidades de t. Por isso os números abaixo usam a janela, interpolada numa grade comum, e os
"−9.5% em t = 50" e "R1″ atrasa mais que o R1′" do §20.6 ficam **retirados**.

| | `R99` médio, t ∈ [35, 48] | Δ vs liso | `dR/dt`, [25, 48] | Δ vs liso |
|---|---:|---:|---:|---:|
| ruído (P2 × P2, duas realizações) | 3.179 / 3.214 | **1.1%** | 0.0758 / 0.0803 | **5.8%** |
| liso (P2R23) | 3.337 | — | 0.0743 | — |
| R1 (grade) | 3.210 | −3.8% | 0.0735 | −1.1% |
| R1′ | 3.136 | −6.0% | 0.0728 | −2.0% |
| **R1″** | **3.146** | **−5.7% (5× o ruído)** | **0.0710** | **−4.4% (no ruído)** |

**Na dose mais baixa, o raio é menor de forma destacável do ruído; a taxa de expansão não.**
R1′ e R1″ não diferem no raio.

**Réguas de forma, média ± desvio sobre os quadros de t ∈ [35, 50]** (`lit_rank.forma`, máscara da
`fig_tese`). **Correção de régua:** a máscara usada pela forma fecha todos os buracos, e com isso um
pilar engolido pela colônia voltava a contar como colônia (área e largura EDT). Os números abaixo
descontam os pilares da máscara; a área ocupada sai nas duas convenções.

| | núcleo / R | baía / R | dedos | área / disco | área / (disco − pilares) | largura EDT (dx) | AR |
|---|---:|---:|---:|---:|---:|---:|---:|
| P2 (realização 1) | 0.43 ±0.01 | 0.42 ±0.01 | 6.7 | 0.25 | — | 2.80 ±0.30 | 15.2 ±3.0 |
| P2 (realização 2) | 0.41 ±0.02 | 0.40 ±0.03 | 6.6 | 0.25 | — | 4.27 ±0.61 | 9.9 ±1.2 |
| liso (P2R23) | 0.42 ±0.03 | 0.40 ±0.03 | 14.0 ±1.0 | 0.43 ±0.03 | 0.43 | 5.48 ±0.13 | 8.1 ±1.4 |
| R1 (grade) | 0.32 ±0.03 | 0.32 ±0.03 | 13.5 ±0.5 | 0.44 ±0.03 | 0.50 | 5.92 ±0.23 | 8.2 ±1.0 |
| R1′ | 0.40 ±0.06 | 0.39 ±0.06 | 14.3 ±1.0 | — | 0.53 ±0.03 | 5.97 ±0.20 | 7.0 ±0.6 |
| **R1″** | **0.40 ±0.06** | **0.41 ±0.05** | **14.0 ±0.9** | **0.48 ±0.02** | **0.54 ±0.03** | **6.16 ±0.26** | **6.4 ±0.3** |

- **Núcleo, baía e número de dedos: iguais ao liso.** O núcleo menor do R1 (0.32) era da parede em grade.
- **Dedos ~12% mais largos** (5.5 → 6.2 dx) e **colônia mais compacta** (área/disco 0.43 → 0.48).
- **AR 8.1 → 6.4 (−21%)**, ainda dentro da faixa de [T12]/PA14 (3.8–11.2). O AR varia 35% entre as
  duas realizações do P2, então **a queda não é atribuível com segurança**; a largura, com desvio
  interno de 0.13–0.26 dx, é a medida mais firme. (O P2 não serve bem de piso para a largura: lá ela
  não é fixada pelo rastro, como no P2R23.)

**Leitura visual (§11; `plots/fig_R1pp_lambda28_anel_shift.png` × `plots/fig_tese_P2R23.png`,
atenção: eixos ±4.5 × ±5):** a classe continua *fingering* (13–14 dedos radiais a partir de um
núcleo). Os dedos ficam mais largos e de contorno mais irregular, e o corpo central engole o primeiro
anel de 6 pilares (r ≈ 1.5). **O painel (c) mostra o mecanismo do raio menor:** o raio mínimo (fundo
da baía) fica **travado em 1.23 de t = 20 a t = 44** e só então salta para 1.85, enquanto no liso
sobe continuamente de 1.0 a 1.65 — a baía mais funda parou contra o anel de pilares.

**Consequência para o pré-registro do R2:** a resposta da dose mais baixa é **raio −5.7%, taxa no
ruído, dedos mais largos, forma preservada**. A monotonicidade R1″ → R2 → R3 vai ser lida
principalmente no `R99` de janela e na largura, não no `dR/dt`, que não separa do ruído nesta dose.

## 21. PRÉ-REGISTRO DO R2 — Λ = 20 dx (2026-10-28, escrito ANTES de rodar)

### 21.1 Especificação

- **Código:** o do R1″, idêntico (`main.py` e `src/` verificados por `diff` contra
  `runs/rugosidade/R1pp_lambda28_anel_shift` em 2026-10-28). **Única mudança: `PILAR_LAMBDA` 28 → 20.**
  Contorno em anel com projeção do shifting (§19, §20); `PILAR_A` = 10, `PILAR_R_EXCL` = 1.2,
  `PILAR_JITTER` = 0, `SEED` = 20260806, `PILAR_SEED` = 20261006.
- **Execução:** mesmo modo do R1′/R1″ (`python main.py -d ...`, sem `--openmp`, 1 thread), t = 50,
  pasta `runs/rugosidade/R2_lambda20`. Duração estimada 25–40 min.
- **Comparação:** contra o liso P2R23 e contra o R1″, com o par P2 × P2 como piso de ruído (§20.7).
  **Toda métrica em janela** (raio: média em t ∈ [35, 48] interpolada; forma: média ± desvio dos
  quadros em t ∈ [35, 50]), nunca num instante só.

### 21.2 Geometria (calculada offline, varredura de 20 000 raios, braço de 5.4–6.2 dx)

| | R1″ (Λ = 28) | **R2 (Λ = 20)** |
|---|---:|---:|
| contatos por braço, alcance 3.3–3.5 | 1.10–1.16 | **1.60–1.91** |
| braços que interceptam ≥ 1 pilar | 84–88% | **100%** |
| fração sólida dentro de r = 3.4 | 0.113 | **0.188** |
| pilares dentro de r = 3.4 | 18 | 30 |
| primeira coroa presente (centro) / superfície | 1.508 / 1.239 | **1.865 / 1.596** (a coroa 1 cai na exclusão, §14.1) |
| colônia chega à primeira superfície (liso) | t ≈ 21 | **t ≈ 27** |

### 21.3 Predições quantitativas

A âncora é o R1″ (§20.6, §20.7). A hipótese de trabalho é a do §15.4/15.5: **o efeito por contato é
constante; o que cresce com a dose é o número de contatos** (×1.5–1.7 do R1″ para o R2).

| grandeza | liso | R1″ | **predição R2** | o que a refuta |
|---|---:|---:|---|---|
| `R99` médio t ∈ [35, 48], Δ vs liso | — | −5.7% | **−7% a −13%** | Δ ≥ −5.7% (não monotônico: checar primeiro o artefato de largada, §14.1, item 4) |
| `dR/dt` [25, 48], Δ vs liso (ruído 5.8%) | — | −4.4% | **−6% a −12%** | dentro do ruído de novo → a taxa não responde a dose nesta faixa |
| largura EDT do dedo (dx) | 5.48 | 6.16 | **≥ 6.2** (6.2–7.0) | < 6.0 → o alargamento do R1″ não era efeito do pilar |
| área / disco | 0.43 | 0.48 | **≥ 0.48** | — |
| dedos | 14.0 | 14.0 | **11–14** | < 11 → extinção por contato mais forte que linear |
| núcleo / R, baía / R | 0.42 / 0.40 | 0.40 / 0.41 | **0.36–0.44** (sem mudança) | — |
| fundo da baía (raio mínimo, painel (c)) | sobe contínuo | travado em 1.23 de t = 20 a 44 | **travado perto de 1.6 a partir de t ≈ 27** por ≥ 10 unidades de t | sem patamar → o travamento do R1″ era acaso |
| razão P(parada \| contato) a 3 dx | — | 1.8× | **1.3–2.5×** (constante com a dose) | cair abaixo de 1.3× → mecanismo novo (retenção, §15.4) |
| razão a 1 / 1.5 / 2 dx | — | 4.6× / 2.7× / 2.2× | **cresce quando o limiar aperta** | forma plana |
| observações em contato a 3 dx | — | 81 | **120–170** | — |
| deixa de ser líder: contato / livre | — | 73% / 23% | **≥ 60% / 20–30%** | — |
| arco com sinal (quiralidade) | — | +0.35 ±0.24 | sem predição; **vigiar** se fica positivo com erro menor | — |
| repulsão química `cos_g` a 0–3 dx | +0.77–0.94 | ≈ 0 | **≈ 0** (propriedade do contato) | positivo como no liso |
| m=6: fase no pico | — | 27.9° | **≈ 0° (mod 60°)**, pico em `R_med`/coroa ≈ 0.99, ou seja `R_med` ≈ 1.85 (§12, §14.1) | fase ~30° → a rede não está travando a forma |
| `a_mar_bio_p95` (mediana t > 30) | 0.73 | 0.63 | **0.50–0.70** | < 0.5 → canal químico (§3, item 2) |

### 21.4 Critérios de aborto (verificar a cada quadro)

Os da §6, mais os herdados do R1″:
1. `n_pen` > 0 em qualquer quadro.
2. força de contato > 0 no ágar antes de a colônia chegar (primeiros 200 passos).
3. partícula pinada a < 2h do anel se aproximando do centro do pilar (> 0.05 dx entre quadros).
4. `a_rep_max` > 3 sustentado em 2 quadros seguidos.
5. pares < 0.05 dx ou iterações por unidade de t acima de 1.5× o R1″ no mesmo t.

### 21.5 Leitura

- **Se as predições de raio, largura e contato se confirmarem:** a resposta à dose é monotônica e
  mediada pela frequência de contato; segue o R3 com o mesmo desenho.
- **Se o raio não piorar em relação ao R1″:** antes de qualquer conclusão, refazer a leitura em
  `R_med/r_coroa` (§14.1) — o R2 tem a largada mais livre das três doses. Só depois a §1 se aplica
  (réplicas).
- **Ferramentas:** `tools/desfecho_contato.py`, `tools/anisotropia.py`, `tools/repulsao_quimica.py`,
  `tools/sombra_quimica.py`, réguas de forma em janela com os pilares descontados (§20.7),
  `plots/fig_tese.py --t 50`.
