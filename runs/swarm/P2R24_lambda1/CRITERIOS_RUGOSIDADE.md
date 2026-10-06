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
