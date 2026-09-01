# Plano de ação — cicatrização da junção núcleo-braço (série K3)

> **STATUS 2026-08-13 — K3 RODADO, REPROVADO e REVERTIDO.** Registrado como **solução
> insuficiente** na lição #59 do CLAUDE.md. `k_col = 0.0`, massa de `BiomassGrowth` de
> volta a `rate·m`; verificado **bit-idêntico ao C4** (frame 200, `max|diff| = 0` em x, y,
> u, v, m, rho, rho_b_grown, cs, c_n, is_filler, h, p). O código da `BiomassColonization`
> com alvo-doador e gate de `cs` fica **preservado desligado** (§10) — o mecanismo está
> certo, falta a maturação (defeito D, lição #58 e
> [ANALISE_TRES_FRENTES.md](ANALISE_TRES_FRENTES.md)). Religar só depois de D.

Estabelecido 2026-08-13, depois de medir o frame `main_01600.hdf5` (t≈21.6 s) do run
C4 em curso. Sucede o [PLANO_JUNCAO.md](PLANO_JUNCAO.md), cujo desfecho (série J,
`k_col` com relaxação-Shepard) foi **revertido** — ver lições #48 a #50 do CLAUDE.md.

Baseline: **`runs/C4`**. Janela de avaliação: **t = 48 s** (lição #47 — além de t≈55
o nutriente esgota e a colônia congela; medir fora da janela mede fome, não morfologia).

Regras válidas em todos os passos: **uma alavanca por passo** (§2.3), **predição
quantitativa ex-ante antes de rodar** (§2.3), **veredito só com frames + log cruzados**
(§11), **`SEED = 20260806` fixa até o fim da série** (§2.5), commit isolado por passo.

---

## 1. O que a medição mostra (frame `main_01600.hdf5`, t≈21.6 s)

Fração de partículas com biomassa, por anel radial:

| r | 0.2 | 0.6 | 1.0 | 1.4 | 1.8 |
|---|---:|---:|---:|---:|---:|
| `frac(ρ_b>0.1)` | 0.921 | 0.385 | **0.079** | 0.063 | 0.017 |
| `frac(ρ_b == 0)` | 0.000 | 0.064 | **0.775** | 0.937 | 0.983 |
| `cs` mediano | 0.409 | 0.359 | 0.294 | 0.216 | 0.154 |
| `c_n` mediano | 0.057 | 0.705 | 0.810 | 0.882 | 0.941 |

Dois fatos decidem o plano:

1. **`c_n` na junção é 0.705**, acima do piso 0.4 do gate de crescimento. O bloqueio
   ali **não é metabólico** — é o estado absorvente `ρ_b = 0` puro (lição #48).
2. **A junção não está quimicamente saturada**: `frac(cs > 0.7·cs_max)` cai de 0.60
   (r=0.6) para 0.003 (r=1.4). Ainda existe gradiente radial exatamente onde qualquer
   preenchimento vai atuar — logo o custo em `cs` é real e tem de ser medido, não
   assumido.

---

## 2. Rotas refutadas (decisão registrada — não repetir)

### 2.1 Ativação por limiar (`ρ_b = 0 → 0.01` se vizinhança densa) — REFUTADA

- **O gatilho não dispara em ninguém.** A média Shepard de `ρ_b` nos buracos vale
  0.025 (r=0.6), 0.012 (r=1.0), 0.006 (r=1.4). Um critério "média da vizinhança > 0.3"
  recruta **zero** partículas.
- **0.01 é sub-quórum → mecanicamente invisível** (lição #53): `BiomassEOS` zera
  `fade_rep`/`fade_att` abaixo de 0.1 ([equations.py:655-657](../src/equations.py#L655-L657)),
  gate flagelar é [0.1, 0.6], `ParticleShift` exige ≥ 0.1.
- **Sair de 0.01 leva 131 s** (`ln(10)/(r_growth·c_n_f)`), contra janela de 50 s.
- **Dispara o runaway da lição #50**: ao cruzar `ρ_b > 1e-12` a partícula passa no gate
  do `BiomassGrowth` e ganha massa a taxa cheia, porque `d_am = rate·m` não é
  proporcional a `ρ_b` ([equations.py:42](../src/equations.py#L42)).

### 2.2 Exclusão mecânica de fases (bactéria empurra ágar) — REFUTADA

- O ágar tem `p ≡ 0`. Na `MomentumEquation` de Monaghan ele **recebe** o empurrão pelo
  termo `p_j/ρ_j²` da bactéria, mas não tem termo próprio contra compressão: **empilha
  em vez de escoar**. `frac_clump` já está em 0.31–0.36 sem empurrão nenhum.
- **Lição #40**: vazio aberto no ágar nunca cicatriza (sem força restauradora) — a rota
  de reciclagem abriu ~300 furos e nenhum fechou.
- Quebra o **C2 do §2.5** onde mais dói: o ágar dentro dos braços **é** o suporte de
  kernel que sustenta `mean_sig_all = 0.929` (Liu §6.5, Violeau §3.6).
- Feito corretamente exigiria multifase real (força de superfície contínua + pressão da
  fase passiva, [T2] van't Hoff) — é a Frente 4/6, não um empurrão unilateral.

---

## 3. Rota escolhida — recrutamento com alvo-doador

### 3.1 O defeito do `BiomassColonization` atual

[equations.py:101](../src/equations.py#L101) relaxa para a média Shepard
(`rho_b_smooth / sigma_a`). Medido sobre os buracos (`ρ_b < 0.01`):

| r | Shepard (alvo atual) | doador `Σ V ρ_b² W / Σ V ρ_b W` | razão |
|---|---:|---:|---:|
| 0.6 | 0.025 | **0.185** | 7× |
| 1.0 | 0.012 | **0.276** | 24× |
| 1.4 | 0.006 | **0.381** | 60× |

A média entre a mãe e o vácuo é o vácuo: o alvo atual recruta para **abaixo do quórum**,
direto na armadilha da lição #53. O alvo correto é a densidade da **mãe** — a filha nasce
com a densidade de quem a gerou (corolário da lição #48).

Um gate `doador ≥ 0.5` recruta **zero** partículas (o doador máximo na junção é 0.44).
O limiar utilizável é ~0.3, ou nenhum limiar, deixando a graduação radial natural.

### 3.2 As quatro mudanças de K3 (revisado após o Passo 0)

1. **Alvo-doador**: em `BiomassColonization`, trocar o denominador `d_sigma_a` por
   `Σ_j V_j ρ_b_j W_ij`. Continua **unilateral** (só soma → o núcleo nunca drena, que foi
   o que matou o D1: `ρ_b` 1.0→0.48, `n_pinned` 43→0) e continua **auto-gateada** (numa
   baia as duas somas → 0). Guarda de 0/0 obrigatória.
2. **Filler fora da soma de doadores** (`s_is_filler < 0.5` no `loop`) — ver §5. Sem
   isso, 93% dos doadores são fantasmas congelados (Passo 0.0).
3. **Gate de saturação química** `cs > 0.6·cs_max` no `post_loop`. Medido: sem ele,
   **570 dos 1260** recrutamentos caem em r > 2, ou seja **na frente**, onde vive o
   gradiente que move a colônia — é o mecanismo exato do afogamento de K2/J5. Com ele
   sobram 286, **todos** em r < 1.12 (mediana 0.95 = a junção), e por construção não
   podem afogar: onde `cs` já está a 60% do teto, o próprio fator `(1 − cs/cs_max)`
   limita o que uma fonte nova acrescenta. É o guard que satisfaz "não silenciar o motor
   de Marangoni" por construção, não por calibração.
4. **Consistência de massa**, em `BiomassGrowth` **e** na colonização:
   `d_am = m · d_a_rho_b_grown / rho_max`. O incremento de massa SPH passa a ser
   proporcional à biomassa efetivamente criada, não à taxa. Desarma a lição #50 **antes**
   que o recrutamento a acione.

São uma alavanca só: (4) sozinho não faz nada, (1) sem (4) dispara o runaway, (1) sem
(2) recruta a partir de fantasmas e (1) sem (3) recruta na frente.

**`k_col` a decidir — 0.03 (herdado da série J) não serve aqui.** A relaxação é
`dρ_b/dt = k_col·c_n_f·(alvo − ρ_b)` e o gate de `c_n` fecha a janela: `c_n_junc` cai de
0.997 (t=0) a **0.370** (t=48), cruzando o piso 0.4 do gate em t≈40. Com janela efetiva
de ~35 s e `c_n_f ≈ 0.8`, `k_col = 0.03` entrega 57% do alvo — **0.091 sobre um alvo de
0.16, abaixo do quórum**, direto na lição #53. Para 80% do alvo na janela é preciso
`k_col ≳ 0.06`. Opções: (i) `k_col = 0.06` na janela atual; (ii) `k_col = 0.03` com
`k_src > 0` para manter o gate aberto — mas aí são duas alavancas, e a série N mostrou
que `k_src` muda o regime inteiro.

### 3.3 A restrição que matou todas as tentativas anteriores (§3.3.6)

`cs_∞ = σ·q·c_n_f / (λ_eff + σ·q·c_n_f/cs_max)`, com os `c_n` medidos por anel:

| zona | ρ_b | q (Hill) | `cs_∞` previsto | `max_cs` medido |
|---|---:|---:|---:|---:|
| núcleo | 1.0 | 0.990 | 0.480 | 0.479 (r=0.2) |
| junção | 0.3 | 0.900 | 0.490 | 0.472 (r=0.6) |
| mid-arm | 0.4 | 0.941 | 0.491 | 0.488 (r=1.4) |
| **sub-quórum** | **0.05** | **0.200** | **0.462** | — |
| ágar | 0 | 0 | 0 (só halo difusivo) | — |

O modelo bate com o medido em três zonas. E diz o essencial: **uma partícula com 5% da
biomassa atinge 94% do `cs` de uma com 100%** — produção por presença, não por conteúdo.
É por isso que K2/J5/J7 afogaram (`|∇cs|` 0.684→0.251, lição #48).

Com produção linear `σ·ρ_b`, `σ ≈ 3`, **mantendo o teto**: núcleo 0.44 / junção 0.42 /
mid-arm 0.44 / sub-quórum 0.234 — discriminação 1.9× contra 1.06× do Hill, perfil ainda
decrescente para fora (§3.3.2 ✓). O teto fica: Y3/lição #56 provou que removê-lo explode
a faixa dinâmica (núcleo `cs` = 14.3, `a_mar` = 116).

---

## 4. Passos

### Passo 0 — Baseline C4 até t=50 ✅ (2026-08-13)

Run completo em t = 50.000 (frame `main_03129`; o `log.csv` para em t=47.8 porque
`print_freq=200` e a última iteração é 3129). Wall time 2378 s. Arquivado em `runs/C4`.

**Valores de referência — toda comparação da série sai daqui:**

| medida | valor | fonte |
|---|---:|---|
| `a_mar_front` (r > 0.75·R99) | **8.605** | `diag_juncao` |
| `a_mar_bio_med` / `p95` | 2.97 / 8.41 | `compare_runs` |
| `constrast_cs` / `mean_cs` / `max_cs` | 12.93 / 0.0379 / 0.4902 | `log.csv` |
| `mass_total` / aceleração | 205.5 / 2.48 | `compare_runs` |
| `n_pinned` | 43 | `log.csv` |
| `a_pressure` (mediana / freq. picos) | 3.00 / **0%** | `compare_runs` |
| `mean_v` | 5.06e-04 | `log.csv` |
| C1 vazio >1.5dx (disco / envelope) | 0.17% (39 dx²) / 0.090% | `compare_runs` |
| C2 `frac(σ_a<0.85)` (disco / **envelope**) | 8.0% / **3.3%** | `compare_runs` |
| componentes conexas da biomassa | **17** (76.7% na maior) | `compare_runs` |
| `rho_b_junc` (p90) / `rho_b_dip` | 0.313 / 0.340 | `log.csv` |
| `c_n_junc` | 0.370 | `log.csv` |
| `iter / t_final` | 62.8 | `log.csv` |
| `n_bio_arms` / `n_ins_arms` | **63 / 1733** | `log.csv` |

**Colapso de dt investigado e descartado:** 62.8 iter/s contra ~52 do C4 histórico
(lição #51) — 1.2×, dentro do guardrail E3 (≤2×). A pausa em t≈17.5–18 foi transitória.

**Ressalva sobre `m.acel = 2.48`:** o guardrail de convergência de massa (≤1.3) só se
aplica a `t ≥ 90` (lição #46c). Em t=50 nenhuma rodada saiu da fase de crescimento —
para K3 o critério vale como **comparação contra este 2.48**, não contra 1.3.

### Passo 0.0 — A colônia visível é 93% filler congelado (medido no Passo 0)

O `diag_juncao` reportou `n_bio = 168` e a primeira versão do `diag_recrut` reportou
2496. A diferença é a definição de população: `diag_juncao` usa `(ρ_b>0.1) & (is_filler<0.5)`.
**Dos 2496 com `ρ_b > 0.1`, apenas 168 são biomassa real; 2328 são filler.** Tanto
`INSERT` ([main.py:886](../main.py#L886)) quanto `WAKE` ([main.py:1052](../main.py#L1052))
inserem com `is_filler = 1.0` e `rho_b_grown` **copiado da mãe** — valor histórico,
congelado, que não cresce, não produz `cs` e não se move.

Distribuição em t=48, separando as populações:

| região | partículas reais | **bio real** (ρ_b>0.1) | filler bio | buracos reais |
|---|---:|---:|---:|---:|
| núcleo r<0.4 | 93 | **93** | 63 | 0 |
| junção r∈[0.4,1.2) | 930 | **38** | 465 | 823 |
| anel médio r∈[1.2,2.0) | 2919 | **9** | 507 | 2910 |
| braços r∈[2.0,3.5) | 9234 | **13** | 1036 | 9221 |

Fora do núcleo existem **60 partículas de biomassa viva** contra **2008 de filler acima
do quórum**. Os braços que aparecem no painel `rho_b` são, em massa, o andaime congelado
— não colônia crescendo. É a lição #41 (dendritos = cauda da gaussiana advectada) levada
ao limite pela inserção: o `n_bio_arms = 63` contra `n_ins_arms = 1733` do log diz o
mesmo, na razão de 27:1.

**Consequência instrumental:** `rho_b_junc` e `rho_b_dip` no `log.csv` são calculados
sobre `fluid.rho_b_grown` **sem filtrar filler** ([main.py:474-493](../main.py#L474-L493)).
O `rho_b_junc = 0.313` do baseline é ~12× otimista em relação ao estado biótico real da
junção (38 vivas contra 465 fantasmas naquele anel). Corrigir antes de usar essas
colunas como critério — hoje elas medem o andaime, não a colônia.

**Censo completo → lição #57 do CLAUDE.md.** São três populações, não duas: o filler do
`WAKE` (2379 dos 2739) **não é pinado** ([scheme.py:64-68](../src/scheme.py#L64-L68)),
então em t=50 há 16 vivas móveis no domínio contra 638 fantasmas em movimento; 90.8% de
`Σ ρ_b·V` está em fantasmas; e `BiomassGradient` não filtra filler, de modo que o `ρ_b`
congelado define o gate de interface da `MarangoniForce`. Também registrada ali a
divergência do pin químico (`c_n<0.6` no código contra 0.4 no §7).

### Passo 0.1 — Instrumento de medição compartilhado ✅

[`tools/diag_recrut.py`](../tools/diag_recrut.py) — mede o **alvo do recrutamento**
(Shepard vs doador por anel, quantos buracos são recrutáveis por limiar, quanto disso
cruza o quórum). Complementa `diag_juncao.py`, que já mede o **defeito** (vale `V`) e o
**motor na frente** (`a_mar_front`, guardrail 4.0). Rodar no C4 **antes** de qualquer
mudança, para que baseline e K3 sejam medidos pelo mesmo instrumento.

### Passo 1 — K3 ❌ REPROVADO em M2 + topologia (rodado 2026-08-13, `runs/K3`)

**Resultado (t=47.7 contra C4 em t=47.8):**

| critério | limiar | C4 | K3 | |
|---|---|---:|---:|---|
| Q1 zeros no disco | <5% | 77.7% | 69.8% | ❌ (mas ver nota) |
| M1 `a_mar_front` | ≥6.88 | 8.605 | **7.236** | ✅ |
| M2 `contrast_cs` | ≥12 | 12.93 | **10.23** | ❌ |
| M3 `a_pressure` med / picos | ≤4 / ≤15% | 3.0 / 0% | 2.9 / 0% | ✅ |
| E1 massa | ≤222 | 205.5 | **203.8** | ✅ |
| E2 `n_pinned` | ≥40 | 43 | 43 | ✅ |
| E3 `iter/t` | ≤126 | 62.8 | 50.3 | ✅ |
| C2 `frac(σ_a<0.85)` geométrica | ≤15% | 3.3% | 3.5% | ✅ |
| — componentes conexas | não piorar | 17 (76.7% na maior) | **59 (42.8%)** | ❌ |

**Q1 falhou pelo critério estar mal escopado, não pelo mecanismo.** Na junção o zero
absorvente foi **eliminado**: `frac(ρ_b==0)` em r=1.0 vai de 0.552 para **0.000** e em
r=1.4 de 0.831 para 0.004. O disco inteiro não melhora porque além de r≈1.6 não há
doador vivo — exatamente o escopo que o próprio plano declarou ("cura a junção, não os
braços"). O critério Q1 deveria ter sido escrito sobre o anel.

**Predição vs resultado (§2.3):** previ o anel indo a 300–330 vivas acima do quórum;
mediu-se ~80. A causa não é o mecanismo — é o **alvo**: `k_col=0.06` entregou `ρ_b` p90
de **0.217** (r=0.6) e **0.097** (r=1.0), contra os 0.134 previstos. A previsão acertou
a magnitude e **errou o lado do limiar**: o quórum é 0.10. O recrutamento tirou todo
mundo do zero e depositou a maioria em (0, 0.1) — lição #53 outra vez, agora por
calibração apertada demais.

**Causa da reprovação em M2:** `cs/cs_max` na junção chega a **0.993 (r=0.6) e 0.983
(r=1.0)** — saturado no teto. O gate `cs>0.6·cs_max` limitava o NÍVEL, mas as recrutadas
sub-quórum produzem por PRESENÇA sob o Hill (`qs(0.1)=0.5`) e cravam o anel inteiro no
teto, achatando `∇cs` no interior. Consequência morfológica medida: `R99` 4.36 → **3.86**
(−11%) e braços mais curtos e atarracados (frames t=50), com a junção visivelmente
preenchida no painel `rho_b` log. É a lição #43 — mais biomassa achata o gradiente.

**Decisão: seguir para o Passo 2A** (produção linear em conteúdo), que era a ramificação
prevista para falha em M2 e agora está duplamente justificada — sob produção linear uma
partícula em `ρ_b=0.1` produziria 10% do máximo em vez de 50%.

<details><summary>implementação (mantida)</summary>

Aplicado: [equations.py](../src/equations.py) `BiomassColonization` (alvo-doador via
`rho_b_w2 / rho_b_smooth`, filler fora da soma, gate `cs > cs_min`, massa ∝ `d_a_rho_b`)
e `BiomassGrowth` (mesma correção de massa); [scheme.py](../src/scheme.py) (`cs_max` e
`col_cs_frac = 0.6` como parâmetros, eliminando o literal duplicado); `rho_b_w2`
registrada em [main.py](../main.py); **`k_col = 0.06`**.

`k_col` derivado do log completo do baseline: `c_n_junc` cruza o piso 0.4 em **t≈45**
(0.997 → 0.370), logo `∫ k_col·c_n_f dt ≈ 0.06 × 31 = 1.86` → **84% do alvo** = 0.134,
acima do quórum 0.10. Com 0.03 seriam 60% → 0.096, sub-quórum (lição #53).

Smoke test de 2 passos: transpila, `mass_total` inicial 197.51 idêntica ao baseline.
</details>

Aplicar §3.2. **Predição ex-ante (t=48, contra C4):**

Escopo revisado pelo Passo 0: com o gate `cs > 0.6·cs_max` o recrutamento atinge **286
buracos, todos em r < 1.12** — cura a junção, **não** os braços. Isso é deliberado: os
braços não têm doadores reais (13 partículas vivas em r∈[2.0,3.5)), e recrutar ali é
recrutar na frente. A predição é sobre o anel r∈[0.4,1.2), não sobre a colônia toda.

| métrica (anel da junção) | C4 | K3 previsto |
|---|---:|---:|
| bio real (ρ_b>0.1) em r∈[0.4,1.2) | **38 / 930** | **300–330 / 930** |
| `frac(ρ_b>0.1)` no anel | 4.1% | **32–36%** |
| zeros reais no anel | 88.5% | **55–62%** |
| `a_mar_front` | 8.605 | **≥ 6.9** (−20% tolerado) |
| `constrast_cs` | 12.93 | **≥ 12** |
| `mass_total` | 205.5 | **≤ 222** (+8%) |
| componentes conexas | 17 | **≤ 17** (não piorar) |

Quantidade que decide o veredito: o alvo mediano no anel é **0.16**, e o quórum é 0.10.
A margem é pequena — se o `k_col` escolhido não entregar ~80% do alvo dentro da janela
de nutriente, o recrutamento pousa **abaixo** do quórum e K3 vira lição #53 de novo.

**Validação — aprova só se TODOS passarem:**

| # | critério | limiar |
|---|---|---|
| Q1 | zeros no disco | < 5% |
| Q2 | `V` na crista | ≤ 0.20 |
| M1 | `a_mar_front` | ≥ 0.8× C4 |
| M2 | `contrast_cs` | ≥ 12 |
| M3 | `a_pressure` mediana / freq. picos >4 | ≤ 4 / ≤ 15% (lição #46b) |
| E1 | massa | ≤ +8%; aceleração 2ª/1ª metade ≤ 1.3 |
| E2 | `n_pinned` | ≥ 40 (guarda contra o modo D1) |
| E3 | `iter / t_final` | ≤ 2× C4 (guardrail invariante de config, lição #51) |
| C2 | `frac(σ_a<0.85)` na população **geométrica** | ≤ 15% |

> C2 na população geométrica, **não** em `ρ_b>0.1`: a lição #56 mostrou que normalizar
> por limiar de biomassa penaliza estruturalmente qualquer rota que crie biomassa
> (18.3% reportado vs 7.4% real).

**Ramificação:** tudo passa → Passo 3. Falha M1/M2 → Passo 2A. Falha Q1/Q2 com M1/M2 OK
→ Passo 2B. Falha E1/E2 → reverter e reabrir (sintoma de que (2) não fechou o balanço,
não de calibração).

### Passo 2A (condicional) — produção linear em conteúdo

Só se K3 afogar o `cs`. Alavanca única: `qs = ρ_b²/(ρ_b²+0.01)` → `ρ_b` em
`SurfactantEquation.post_loop`, com `σ: 10 → 3.0`, **mantendo o teto**.

**Obrigatório antes de editar:** rodar `/protocolo-cs-zonas` (§3.3.6) e registrar a
tabela `cs_∞`. Se o pico migrar para o rim externo ou para o ágar, **não aplicar**
(§3.3.4). Validação: critérios do Passo 1 + `max_cs ≤ 1.1·cs_max` (guarda contra a
explosão de faixa dinâmica da lição #56).

### Passo 2B (condicional) — subir a taxa de recrutamento

Só se K3 preservar o motor mas encher pouco. `k_col: 0.03 → 0.06`, nada mais. Atenção
redobrada a M1/M2: foi entre `k_col` 0.03 e 0.1 que o motor caiu 49% na série J (#49).
Se M1 reprovar, o ponto de operação é 0.03 e a limitação é a lei de produção → 2A.

### Passo 3 — Validação completa §2.5 + §11

```bash
tools/archive_run.sh K3
python tools/compare_runs.py runs/C4 runs/K3
python tools/compare_frames.py runs/C4 runs/K3
python tools/diag_juncao.py runs/C4 runs/K3 --t 48
python tools/diag_recrut.py runs/C4 runs/K3 --t 48
python tools/validate_model.py runs/K3
```

Protocolo §11 nas três etapas (inspeção visual isolada é falha analítica): descrição
quantitativa dos frames → cruzamento com `log.csv` → síntese com veredito explícito.

**Critério de sucesso da série:** junção contínua **e** AR ≥ 5 preservado **e** C1/C2 não
regredidos. Cicatrizar a junção às custas de braços atarracados é **reprovação** — foi o
que a série N fez (AR 5.3 → 3.0).

### Passo 4 — Extensão a t=100 (só depois de passar em t=50)

`total_sim_time: 50 → 100`, nada mais. A lição #34 é explícita: a validação de t=50
**enganou** na C3.2 (auto-limitante até t≈57, runaway depois). Procurar especificamente:
inserções cravando no cap, massa acelerando na 2ª metade, `contrast_cs` caindo **com**
`mean_cs` subindo junto (só as duas juntas são afogamento — lição #47).

Ressalva: com `k_src = 0` o nutriente esgota e a colônia congela em t≈55 (lição #47, com
a atribuição corrigida pela #51 — o gargalo é contagem de portadores). Se K3 recrutar
portadores, essa trava pode se mover: medir `c_n` nos braços e `dR/dt` nas janelas
[33,55] e [55,95] para ver se o expoente α saiu de 0.40.

### Passo 5 — Registro

`/registrar-pass` → CLAUDE.md §9: predição ex-ante vs resultado real, causa-raiz, lição
numerada. **Registrar as refutações também** (§2.1 e §2.2 deste plano) — poupa rodadas
futuras, como a rota do inóculo fechada por predição na lição #50.

---

## 5. Filler como doador — MEDIDO, decisão tomada (2026-08-13)

`BiomassColonization.loop` soma **todas** as fontes, sem checar `s_is_filler`
([equations.py:80-83](../src/equations.py#L80-L83)), enquanto o `post_loop` exclui filler
como destino. O filler é inserido com `ρ_b ≥ 0.5` herdado e congelado
(`INSERT_RHO_B_MIN = 0.5`), então sob o alvo-doador ele viraria um **doador forte e
imóvel** — valor histórico que não decai, a armadilha da lição #39 (onde o filler
congelado com `cs` histórico inflava `mean_cs`).

Medido (`tools/diag_recrut.py . --t 21.6`), alvo mediano nos buracos:

| r | doador (todas as fontes) | doador **sem filler** |
|---|---:|---:|
| 0.6 | 0.185 | 0.102 |
| 1.0 | 0.278 | 0.107 |
| 1.4 | 0.381 | **0.000** |

A diferença é material, e no anel r=1.4 é total: a mediana cai a zero, ou seja **os
doadores ali são quase inteiramente filler**. Confere com o log (`n_bio_arms = 75`
contra `n_ins_arms = 391` nos braços) — a biomassa real nos braços é um esqueleto de
~75 partículas cercado por ~390 inseridas congeladas.

**Decisão: excluir filler da soma de doadores** (`s_is_filler < 0.5` no `loop`). Custo
medido: recrutáveis caem de 2631 para **1406** (gate 0.1), com alvo mediano **0.362** —
ainda bem acima do quórum. O que se perde são exatamente os recrutamentos que seriam
fundados num valor histórico congelado, não em biomassa viva.

**Consequência para a predição do Passo 1:** o recrutamento fica confinado à vizinhança
da biomassa real, então `frac(ρ_b>0.1)` em r≈1.4 sobe menos do que o previsto em §4. Se
Q1/Q2 reprovarem **por falta de doadores reais** (e não por motor), o diagnóstico não é
`k_col` baixo — é que o esqueleto de 75 portadores não tem alcance, o que remete à
lição #51 (contagem fixa de portadores) e ao Passo 2B/2A, não a afrouxar o gate.
