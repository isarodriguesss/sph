# Plano de ação — fluxo quimiotático de biomassa

> **STATUS 2026-08-12 — X1 REPROVADO, rota encerrada na forma proposta.** O fluxo
> transporta biomassa a território novo (60% dos zeros preenchidos em 5 s) e mantém
> `cs` estruturado (faixa 3,7×, passou no discriminante do X0), mas **dilui**: a banda
> do gate flagelar caiu de 2041 para **34** partículas e a colônia congelou em
> `R99 = 0.49` contra 4.39 do C4. Ver lição #54.
>
> **C4 permanece o baseline.** `chi = 0.0`; `ChemotacticFlux` preservada desligada.
>
> **X2 e X3 perdem o sentido** com este resultado: `NutrientSource` não ataca diluição,
> e o `K` do Hill mexe na produção de `cs`, que não foi o modo de falha. **X4** (`dx`
> menor) segue válido como paliativo independente. **X5/X6** inalterados.
>
> **X1b (limite de capacidade, `ρ_target=0.4`) TAMBÉM REPROVADO** — banda flagelar 52,
> `R99` 0.418. O limite restringe o **nível** do receptor, não o **número** deles: com
> receptores em `ρ_b≈0` a condição quase nunca ativa. **Isso fecha a classe de
> transporte**: com orçamento fixo, `ρ_b ~ Σρ_b·V/(N·V)`, e nenhuma regra par-a-par
> muda `N`.
>
> **X5 analisado preditivamente e NÃO resolve o problema** (custo zero, sem rodada):
> o influxo van't Hoff traz **solvente**, não células. Vira `dm/dt`, não toca `rho_b`,
> então a colônia ganha área com o mesmo orçamento e a concentração **cai** — espalhado
> em 1×/2×/4× a área, o `ρ_b` médio dá 0.0042/0.0021/0.0011, todos abaixo do quórum.
> Além disso `V₀` é máximo no núcleo (pinado) e **diverge em `φ→1`**, onde o núcleo está.
> X5 segue valendo como **motor de expansão** ([T5]: osmótica é dominante), mas para
> outro problema.
>
> ## Reorganização — o que resta
>
> A única coisa que aumenta o orçamento de biomassa é o `BiomassGrowth`. N1 já provou
> que a fonte de nutriente o sustenta (2,15× em 100 s, `c_n` nunca esgota); em ~300 s
> daria ~10×, suficiente para popular os braços. O bloqueio é a restrição #1 — mais
> biomassa afoga o `cs` — e a alavanca contra ela **nunca foi testada**.
>
> **Caminho proposto:** `k_src = 0.3` (validado em N1) + domínio `[-12,12]` + `t ≈ 120`
> + **`K` do Hill 0.1 → 0.4 com `σ` recalibrado** (o X3 deste plano, que era
> contingência e passa a ser o eixo). Custo ~2,3 h por rodada.
>
> **Alternativa estrutural:** gates **relativos** em vez de absolutos — tornaria a
> diluição irrelevante por construção, mas toca todas as equações (patamar X6).
>
> Permanece válido da análise abaixo: o diagnóstico do problema, os critérios de
> sucesso (§ Critérios), o protocolo de validação e o X0 — cujo saldo foi positivo
> apesar do desfecho (pegou o gate do núcleo, que evitaria repetir o D1, e o `χ`
> calibrado com o gradiente errado). **Duas correções ao X0 estão na lição #54:** o
> discriminante precisava incluir uma medida de expansão (`R99`), e curvatura medida
> no estado final não prediz a dinâmica a partir do estado inicial.

Estabelecido 2026-08-12. Baseline: **`runs/C4`**. Janela de leitura: **t ≤ 50 s**.

Ataca a descontinuidade documentada nas lições #48–#53: as partículas da borda do
inóculo nascem em `ρ_b = 0` exato (underflow da cauda quártica além de r≈1.25–1.88)
e o crescimento multiplicativo (`rate·ρ_b`) torna zero um **estado absorvente**.
Provado causalmente: 100,000% das 65 400 partículas em zero em t=50 nasceram em zero;
**0 de 65 400** adquiriram biomassa; o piso `1e-300` (teste F1) deu rodada
bit-idêntica, logo não é artefato numérico.

---

## Por que esta rota e não as outras

O §3.0 registra que **[T3] Giverso, Verani & Ciarletta 2016** compara crescimento
**volumétrico** (`Γ = K_γ·ρ·n`) contra **fluxo quimiotático** (`m = χ·ρ·∇n`), e que
o segundo produz *"padrões mais simétricos com múltiplos dendritos"* — o nosso alvo.

O modelo atual implementa os dois pela metade: o crescimento é volumétrico, e a
quimiotaxia existe apenas como **força sobre partículas** (`FlagellarForce`), nunca
como **fluxo do campo de biomassa**. Como `ρ_b` é escalar passivo advectado, a
quimiotaxia move sempre os mesmos ~150 portadores e nunca leva biomassa a território
novo. **É essa metade faltante que este plano implementa.**

Vantagem estrutural sobre tudo que já foi reprovado — o fluxo **conserva massa**:

| restrição estrutural | colonização (K2/J5/N2) | fluxo quimiotático |
|---|---|---|
| #1 Hill produz por PRESENÇA | encher multiplica produtores **10,9×** | conserva `∫ρ_b`; enche em `ρ_b` baixo → **~+34%** |
| #6 `d_am ∝ m` (runaway de massa) | abre o gate de massa p/ o domínio | **não cria massa** |
| alargamento dos braços (AR 5.3→3.0 em N2) | relaxa isotropicamente | segue `∇c` → **anisotrópico por construção** |
| pontas descolando (35% dos fragmentos além do p90) | não ataca | move biomassa **para** a ponta |

---

## Critérios de sucesso — definidos ANTES de calibrar (lição #46)

População e estatística fixadas aqui; nenhuma rodada muda a definição depois.

| | métrica | população | baseline C4 | alvo |
|---|---|---|---:|---:|
| **P1** | `frac(ρ_b>0.1)` na junção `r∈[0.8,1.4]`, dentro dos cones dos braços | geométrica | 20–30% | **≥ 50%** |
| **P2** | componentes conexas da biomassa / fração na maior | `ρ_b>0.1`, arestas <2h | 32 / 71% | **≤ 15 / ≥ 85%** |
| **M1** | `a_mar_frente` (mediana em `r>0.75·R99`) | biomassa | 7.23 | **≥ 6.0** |
| **M2** | `contrast_cs` | — | 12.9 | **≥ 10** |
| **F1** | AR (`validate_model.py`) | — | 5.3 | **≥ 5.0** (§2.2) |
| **F2** | nº de dendritos | — | 17 | **15–20** |
| **V1** | vazio >1.5dx, **envelope** | geométrica | 0.126% | **≤ 0.15%** |
| **V2** | `frac(σ_a<0.85)`, **envelope** | geométrica | 3.2% | **≤ 6%** |
| **V3** | aceleração de massa (2ª/1ª metade) | — | 2.61 | **≤ 3.0** |

**P1 e P2 são o problema.** M1/M2 e F1/F2 são o que não se pode perder. V1–V3 são o
seu requisito de não regredir em vácuo (reportados sob envelope, que é a definição
geométrica — o disco R99 infla a área 3,2× contando baias).

**Reprovação automática:** qualquer rodada que satisfaça P1/P2 violando F1 repete o
N2 (encheu e engordou) e não conta como avanço.

---

## X0 — Análise preditiva (custo zero, OBRIGATÓRIA antes de X1)

Sem rodar nada, produzir:

1. **Forma discreta conservativa.** `∂ρ_b/∂t = −∇·(χ·ρ_b·∇c)`. Exigir forma
   antissimétrica par-a-par para que `Σ_i ρ_b,i·V_i` seja conservado exatamente —
   citar **Liu §3.4** e **Violeau §5.3** (conservação sob formulação antissimétrica).
   Verificar numericamente a conservação num teste isolado antes de acoplar.
2. **Escolha do atrator.** `∇c_n` (fidelidade a [T3]) ou `−∇cs` (consistência com
   `FlagellarForce`). Ambos apontam para fora; decidir e justificar.
3. **Estimativa de `χ`.** Alvo: transportar biomassa na escala do braço em ~30 s.
   Com `|∇c_n| ≈ 0.22` medido no C4 em r∈[1.2, 3.0] e velocidade de ponta 0.09/s:
   **`χ ≈ 0.4`**. Reportar a faixa `[0.2, 0.8]` como bracket.
4. **Protocolo §3.3.6** — tabela `cs_∞` em 4–5 zonas com o fluxo ligado. É proibição
   explícita do §10 mexer em produção/sumidouro de `cs` sem isso; aqui a produção
   muda indiretamente (mais partículas acima do quórum), então vale.
5. **Tempo de drenagem do núcleo.** `τ ≈ L_núcleo/(χ|∇c_n|)`. Se `τ < 30 s`, X2
   (nutriente) passa a ser pré-requisito e não contingência.

**Sem esses cinco itens escritos, não rodar X1.**

### X0 — EXECUTADO (2026-08-12)

**1. Forma discreta conservativa.** Velocidade de deriva `u_i = χ·(−∇cs_i)` (já
disponível: `FlagellarForce` acumula `grad_cs_x/y`). Fluxo `J = ρ_b·u`, e

```
dρ_b_i/dt = −Σ_j V_j (ρ_b_i·u_i + ρ_b_j·u_j)·∇W_ij
```

A forma é **antissimétrica par-a-par** — o termo do par (i,j) em i é o negativo do
termo em j, porque `∇_j W_ji = −∇_i W_ij`. Logo `Σ_i d(ρ_b_i V_i)/dt = 0` exatamente
enquanto `V_i` for constante no passo (Liu §3.4; Violeau §5.3 prova a conservação sob
formulação antissimétrica). **Verificar `Σ ρ_b·V` numericamente antes de acoplar** —
é a propriedade da qual todo o argumento depende.

Propriedade-chave: uma partícula com `ρ_b = 0` tem fluxo de saída **zero** mas recebe
do vizinho com biomassa. O estado absorvente deixa de existir por transporte, sem
termo aditivo arbitrário.

**2. Atrator: `−∇cs`** (não `∇c_n`, apesar de [T3] usar o nutriente). Três razões,
medidas no C4 em t=48: (a) na junção `r∈[0.9,1.5]` — que é o alvo P1 — `|dcs/dr|`
vale **0.29** contra `|dc_n/dr| ≈ 0.19`; (b) `−∇cs` sobrevive ao X2, enquanto o
`NutrientSource` achata `∇c_n` (medido em N1: gradiente cai ~2×); (c) é a direção que
o §3.2 Frente 5 já declara como a quimiotaxia do modelo ("quimiotaxia para agar
fresco"). Desvio de [T3] documentado: a classe do mecanismo é a mesma, o atrator não.

**3. `χ = 0.5`**, bracket `[0.15, 1.0]`. Calibrado para a deriva acompanhar a ponta
(`dR/dt = 0.099`) na banda de gradiente máximo: `χ = 0.099/0.19 = 0.51`.

**4. Gate `ρ_b < 0.8`** — exclui o núcleo pinado. Sem ele o fluxo **atravessa o pin**
(o pin zera velocidade da partícula, não o fluxo do escalar) e drena o núcleo em
`τ = 0.35/(χ·|dcs/dr|) = 41 s` com `χ=0.5` — dentro da janela. Justificativa
biológica: §2.4, "núcleo denso imóvel (matriz EPS madura)"; célula madura não
quimiotaxa. Com o gate, o modo de falha do D1 fica fechado por construção.

**5. Protocolo §3.3.6 — AQUI ESTÁ O RISCO PRINCIPAL.**

| zona | `ρ_b` | `c_n` | **C4 medido** | **se o fluxo espalhar** |
|---|---:|---:|---:|---:|
| núcleo (r<0.5) | 0.80 | 0.10 | 0.78 | 0.97 |
| junção (r≈1.0) | 0.15 | 0.40 | **0.76** | 0.96 |
| braço (r≈2.0) | 0.35 | 0.55 | **0.28** | 0.97 |
| ponta (r≈3.0) | 0.45 | 0.78 | **0.18** | 0.96 |

*(valores de `cs/cs_max`)*

O perfil do C4 varia **4,3×** entre núcleo e ponta; se o fluxo dobrar a fração de
produtores, vira **plano em 0.96–0.97** — top-hat, gradiente radial zero.

A razão é aritmética: `cs/cs_max = 1/(1 + λ_eff·cs_max/P)` com `P = σ·qs·c_n_f`. Com
`σ=10`, `P ≈ 8·qs`, e `λ_eff·cs_max = 0.15`, então **satura sempre que `ρ_b > 0.014`**.
O campo de `cs` é top-hat por construção; o que o mantém graduado hoje é que a maioria
das partículas **não produz nada** (o esqueleto esparso).

**Consequência para o plano: X3 pode ser pré-requisito, não contingência.** Para `cs`
ficar graduado seria preciso `σ` entre 0.2 e 1.5 — território do J5, onde o motor
morreu.

**6. Contra-argumento medido — a curvatura de `cs` favorece a rota.** `div(u) = −χ∇²cs`:

| r | 0.6 | **1.2** | **1.8** | **2.4** | 3.0 | 3.6 |
|---|---|---|---|---|---|---|
| | diverge | **converge** | **converge** | **converge** | diverge | converge |

O fluxo **acumula** biomassa exatamente em `r ∈ [1.2, 2.4]` — a região que precisa
ser preenchida — e espalha junto ao núcleo, de onde ela deve sair. Se a concentração
prevalecer sobre o espalhamento, os produtores ficam **mais** localizados e `cs`
**afia** em vez de achatar, invertendo o risco do item 5.

**Conclusão de X0: X1 é um teste discriminante, não uma aposta.** Os dois desfechos
são distinguíveis por uma única medição — o perfil radial de `cs/cs_max`:

- **achatou** (junção e ponta ambas > 0.9, faixa < 1.5×) → restrição #1 pela quinta
  vez; a rota exige X3 antes de qualquer coisa, e o par `σ`/`K` passa a ser o
  problema central;
- **afiou ou preservou** (faixa ≥ 3×, ponta < 0.4) → a convergência venceu, e a rota
  segue para P1/P2 sem precisar de X3.

Autorizado rodar X1 com `χ=0.5`, gate `ρ_b<0.8`, atrator `−∇cs`.

---

## X1 — Fluxo quimiotático isolado *(risco médio, 1 rodada)*

`χ` da estimativa de X0, sobre o C4, tudo o mais inalterado, t=50.

**Predições a registrar antes:** P1 20–30% → **≥45%**; P2 32 comp → **≤20**;
`∫ρ_b` conservado (variação < 1% descontado o `BiomassGrowth`); `mass_total` ≈ C4
(+4,3%, o fluxo não cria massa); M1 ≥ 6.

**Modos de falha esperados e o que cada um significa:**
- **núcleo drena** (`n_pinned` 43 → <20, `ρ_b` máximo < 0.8) → é o modo do D1; vai
  para X2.
- **`cs` infla** (`contrast_cs` < 8, `cs_bio_arms` > 0.47) → a convexidade de `qs`
  mordeu mesmo com massa conservada; vai para X3.
- **braços alargam** (AR < 4.5) → o fluxo não está sendo anisotrópico na prática;
  revisar a forma discreta antes de qualquer calibração.
- **nada muda** (P1 inalterado) → `χ` baixo demais; subir dentro do bracket [0.2, 0.8].

---

## X2 — Fluxo + `NutrientSource` *(risco baixo, 1 rodada)*

Só se X1 drenar o núcleo. `k_src = 0.3`, já validado em N1 (lição #51: `min_c_n`
nunca abaixo de 0.376, biomassa 2,15×, vácuo melhor que o baseline, frente +23%).

**Atenção — dois efeitos colaterais medidos em N1**, que precisam entrar na leitura:
(a) o pin químico `c_n < 0.6` afrouxa, despinando tudo com `ρ_b < 0.40`; (b) `c_n`
uniforme achata o gradiente que estruturava a produção de `cs` e `a_mar_bio_med` cai
(medir pela FRENTE, nunca pela mediana).

---

## X3 — `K` do Hill: `qs = ρ_b²/(ρ_b² + K²)`, `K` 0.1 → 0.4 *(risco médio-alto)*

Só se X1/X2 inflarem `cs`. Torna `qs` **graduado** na faixa de operação em vez de
degrau: `qs(0.05)` cai de 0.20 para 0.015, `qs(0.3)` de 0.90 para 0.36.

**Não é neutro:** a produção nos braços cai ~2,5×, exigindo `σ` 10 → ~25 no mesmo
passo. É par, não alavanca — declarar isso e prever pelo protocolo §3.3.6.

Ataca a restrição #1, que reprovou K2, J5, J6 e N2 — quatro rotas pelo mesmo
mecanismo. É a alavanca de maior alcance entre as baratas.

---

## X4 — Resolução global mais fina *(risco baixo, custo alto)*

`dx` 0.0538 → 0.036 (grade 261² → ~390², 2,2× partículas). `dt` ~0.6×, wall ~3,4×
(≈70 min por rodada).

Não resolve o estado absorvente, mas os braços deixam de ser esqueletos de ~44
partículas, e todas as métricas de cobertura melhoram por orçamento de partículas.
**Paliativo honesto** — vale como validação final da melhor configuração, não como
tentativa de conserto.

---

## X5 — Duas fases [T2] Srinivasan *(risco muito alto, reescrita)*

Substituir `ρ_b` por fração volumétrica `φ` com balanço de massa próprio e influxo
van't Hoff `V₀ = Q₀·(φ/(1−φ) − φ₀/(1−φ₀))`. É o **Objetivo 2 do §1** (pressão
osmótica); a `OsmoticForce` está no repositório desativada esperando isso, e o §10
proíbe removê-la justamente por ser roadmap.

Resolveria o problema **por construção** — `φ` teria transporte próprio em vez de ser
escalar advectado. Invalida a calibração das séries K/M/T/C inteira.

## X6 — Separar ágar e colônia em fases distintas *(risco muito alto)*

A raiz conceitual: num corte 2D de camada única não existe "embaixo". Medido: o ágar
não é empurrado (deslocamento **+0.0004** em `r₀∈[2.5,3.5]`) porque tem pressão zero
(`BiomassEOS`, `ρ_b<0.1`) e Marangoni zero (gate `|∇ρ_b|≥0.15`). Num swarm real as
bactérias se espalham **sobre** o ágar. Reescrita maior que X5.

---

## Protocolo de validação — idêntico em todo passo

```bash
make run
tools/archive_run.sh X1
python tools/validate_model.py runs/X1                 # ancoras [T1][T2][T3][T6][T7]
python tools/compare_runs.py runs/C4 runs/X1           # C1/C2 nas DUAS definicoes
python tools/diag_juncao.py runs/C4 runs/X1 --t 48 --out runs/juncao_X1.png
python tools/plot_envelope.py runs/C4 runs/X1 --times 15 30 48 --out runs/env_X1.png
python tools/plot_fields.py runs/X1 --times 15 30 48 --out runs/campos_X1.png
```

Cada passo compara contra **o passo anterior E contra o C4**. Registro na §9 com
predição *ex-ante* × resultado medido; discrepância > 2× é investigada antes de
avançar (§2.3).

**A leitura visual é o árbitro (§11).** Três vezes nesta sessão uma métrica escalar
disse o oposto do painel: `a_mar_bio_med` contaminado por composição, `frac(σ_a<0.85)`
penalizando quem recruta biomassa, e o halo de `cs` do `validate_model.py` reportando
melhora quando o campo havia saturado.

---

## Regras de parada

- **Avançar** só com o critério do passo satisfeito e sem violar M1/M2/F1.
- **Parar e reportar** se dois passos consecutivos falharem pelo mesmo mecanismo —
  foi o padrão de K2/J5/J6/N2 (quatro rotas, sempre a restrição #1).
- **Não escalar para X5/X6** sem decisão explícita: são projeto de tese, não
  calibração, e invalidam a calibração existente.
- **C4 permanece o baseline** até um passo satisfazer P1+P2 sem violar F1.
