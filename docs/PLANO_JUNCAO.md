# Plano de ação — vale de biomassa na junção núcleo-braço

> **STATUS 2026-08-11 — RESOLVIDO por `k_col = 0.03` (J7).** O plano abaixo (J1–J4) foi
> escrito antes de medir o campo `rho_b` e está **superado**: ele mirava o vale no perfil
> da crista, quando o defeito real era o **estado absorvente `rho_b = 0`** em 77.6% do
> disco. J7 elimina 100% dos zeros custando 1.2% do motor. Ver §9 lições #48 e #49 do
> CLAUDE.md para a série completa (K2, J5, J6, J7, J8) e o bracket de `k_col`.
>
> O que da análise abaixo permanece válido: o diagnóstico da cadeia causal (§1), a regra
> de projeto (§2), o protocolo de validação (§5) e a tabela do que já foi reprovado (§6).
> J1 (pin `c_n<0.6`→`0.4`) e J2 (`D_n`) **nunca foram rodados** e seguem em aberto —
> nenhum dos dois toca o zero absorvente, mas J1 continua valendo como higiene da §7.

Estabelecido 2026-08-11. Baseline: **`runs/C4`** (§7, config `SHIFT_CAP=0.0006`,
`D_n=0.05`, domínio `[-7,7]` 261²). Janela de avaliação: **t = 48 s** (§ lição #47 —
além de t≈55 o nutriente esgota e a colônia congela; comparar fora dessa janela mede
fome, não morfologia).

---

## 0. Correção de registro (feita antes de qualquer coisa)

`main_output/` e `log.csv` na raiz **não eram C4 — eram a rodada K2** (`k_col=0.3`,
reprovada). Confirmado por igualdade até o último decimal com `runs/K2` em `t`,
`mean_v`, `a_marangoni`, `mass_total` e `biomass_total`. Toda medição feita sobre
esse diretório descrevia K2, não o baseline.

Consequência prática: **antes de qualquer rodada nova, limpar `main_output/`.**
`make run` já limpa, mas o arquivamento (`tools/archive_run.sh`) deve ser imediato.

---

## 1. O que o baseline C4 realmente mostra

Medido por `tools/diag_juncao.py` sobre `runs/C4/main_output` em t=48 s
(figura: `runs/juncao_baseline.png`).

| | C4 (baseline) | K2 (`k_col=0.3`, reprovada) |
|---|---:|---:|
| `n_bio` (partículas com `rho_b>0.1`) | **167** | 606 |
| `V` — profundidade do vale | **0.390** | 0.360 |
| `rho_b` mínimo na crista | 0.222 | 0.270 |
| `rho_b` no braço (r>2) | 0.364 | 0.422 |
| `c_n` na junção | 0.361 | 0.373 |
| `cs/cs_max` na junção | **0.729** | **0.973** |
| fração congelada da biomassa | 0.904 | 0.705 |
| **`a_mar_bio_med`** | **3.145** | **0.126** |
| vazio >1.5dx / `frac(σ_a<0.85)` | 0.16% / 7.4% | 0.44% / 16.2% |
| `compare_runs.py` | aprova | **REPROVA** (C2 16.2%, contrast 4.2) |

**O defeito não é vazio de partícula.** C1 e C2 de §2.5 estão satisfeitos (vazio
0.16%, `σ_a` 7.4%). A colônia é materialmente contínua. O defeito está no **campo
escalar `rho_b`**: a crista cai de 0.80 (r=0.3) para 0.31 (r=0.5) e fica em 0.28–0.45
até r=2.5. É um **degrau na borda do núcleo**, exatamente onde a cauda quártica do
inóculo desaba.

**Cadeia causal medida** (painel `c_n` da figura):

| fronteira | raio | efeito |
|---|---:|---|
| `c_n` cruza 0.4 (gate de `BiomassGrowth`) | r ≈ **0.95** | dentro disso, crescimento **exatamente zero** |
| `c_n` cruza 0.6 (pin cinemático, [scheme.py:64](../src/scheme.py#L64)) | r ≈ **2.2** | dentro disso, `u=v=0` — **nenhuma advecção** traz biomassa |
| degrau do inóculo `exp(-(r/R_θ)⁴)` | r ≈ **0.4–0.5** | semeia `rho_b ~ 5e-4`; crescimento é multiplicativo (`rate·rho_b`) → **nunca recupera** |

Os três se sobrepõem na zona da junção. É por isso que nada a preenche sozinho.

---

## 2. Regra de projeto (derivada da medição — vale para toda alavanca)

`qs = rho_b²/(rho_b²+0.01)` satura em `rho_b ≈ 0.1`: vale 0.89 em `rho_b`=0.28 e
0.99 em `rho_b`=1.0. Logo:

> **Subir `rho_b` em partícula que JÁ é biomassa é quimicamente quase grátis
> (`qs` já ≈1). Recrutar partícula VAZIA para `rho_b > 0.1` é o que satura `cs`.**

K2 recrutou: `n_bio` 167 → 606 (3.6×) → `cs/cs_max` na junção 0.73 → 0.97 → o
gradiente sumiu → `a_mar_bio_med` 3.15 → 0.13. Foi assim que ela conseguiu melhorar
`V` em 8% destruindo o motor em 25×.

**Métrica de controle obrigatória em todo passo: `n_bio` ≤ 200** (baseline 167).
Rodada que multiplique `n_bio` está reproduzindo K2, independente do que `V` diga.

---

## 3. Critério de sucesso — definido ANTES de calibrar (lição #46)

| | métrica | alvo |
|---|---|---|
| **A** | `V = 1 − min_r(ρ_b crista) / mediana(ρ_b crista, r>2)` | **≤ 0.20** (C4 = 0.39) |
| **B** | `a_mar_bio_med` (§2.5.1 C3) | **≥ 2.5** (C4 = 3.15; perder ≤20%) |
| **C** | `n_bio` | **≤ 200** (C4 = 167) |
| **D** | C1/C2 + guardrails de `compare_runs.py` | aprovar |

Os quatro são obrigatórios e indissociáveis. **A sem B** é K2. **A sem C** é K2 por
outro caminho. `V` normaliza pelo nível do braço justamente para que inflar biomassa
global não "melhore" o vale sem reconectar nada.

---

## 4. Passos — uma alavanca por vez (§2.3)

### J1 — pin químico `c_n < 0.6` → `< 0.4`

[scheme.py:64](../src/scheme.py#L64). É **revert ao valor documentado** na §7
("hard cutoff < 0.4", M-B.8) — o código divergiu sem registro na §9.

*Por que primeiro:* `c_n` tem dupla função (gate de crescimento em 0.4, pin em 0.6).
Mexer no nutriente (J2) antes de fixar o pin move os dois eixos de uma vez e destrói
a atribuição — é o erro da lição #12 (D_ext/λ_ext simultâneos).

**Predição:** fração congelada da biomassa 0.904 → ~0.83; a faixa r∈[1.2, 2.2]
descongela (100% → ~10%); `V` praticamente inalterado (0.39 ± 0.04 — o pin não cria
biomassa); `n_bio` inalterado.
**Risco:** lição #22 — afrouxar pin sem coesão proporcional dá fragmentação
distribuída. Vigiar `frac_clump` (C4 = 0.40), `a_pressure` (mediana 3.0, picos 13%)
e `max_v`.
**Reverte se:** picos de `a_pressure` > 25% das amostras, ou `frac_clump` > 0.55, ou
partículas soltas nos frames.

### J2 — `D_n: 0.05 → 0.15`

[main.py:130](../main.py#L130), com o pin já fixado em 0.4.

*Por que é a alavanca certa:* o gate de `c_n` já está **saturado em 1.0 na frente**
(r>3.2) e **em 0 na junção** — subir `c_n` não pode amplificar a frente, só religa a
junção. E religa por **crescimento multiplicativo**, que não recruta partícula vazia
(quem tem `rho_b=0` continua em 0 para sempre) — satisfaz a regra da §2 por
construção. Precedente: C1 (lição #45) subiu `D_n` 0.02→0.05 e levou `c_n` nos braços
de 0.47 → 0.61.

**Predição:** `L_D_n` na junção `√(D_n/(k_n·ρ_b))` 1.0 → 1.73; `c_n_junc` 0.361 →
**0.55–0.65**; gate de crescimento em r∈[0.5,1.2] de 0 → **0.3–0.5**; crista sobe de
0.22 → **0.40–0.50** → `V` 0.39 → **0.15–0.25**; `n_bio` 167 → 170–200;
`a_mar_bio_med` ≥ 2.6 (a junção sobe `qs` de 0.89 → 0.96, +8% de produção — longe do
salto que saturou K2).
**Reverte se:** `n_bio` > 220 ou `cs/cs_max` na junção > 0.92.

### J3 — inóculo com platô (só se J2 fechar B/C mas não A)

[particles.py:39-40](../src/particles.py#L39-L40):

```python
R0 = 0.45 + 0.09 * np.cos(8 * theta)
rho_b = np.exp(-((np.maximum(0.0, dist - R0) / 0.15) ** 4))
```

**Não confundir com I1** (`p=4→2`, R maior), que foi reprovado: I1 suavizou a
*inclinação* (biomassa +2.5×, `cs` uniforme, `a_mar` 2.30→0.75). Aqui a inclinação
fica **1.3× mais afiada** (borda 0.9→0.1 em Δr=0.10 contra 0.20 hoje) e só o
**platô** é alargado — preenche o degrau em r∈[0.3, 0.45] sem espalhar produção.
Custo estimado: biomassa inicial ~1.9× (contra 2.5× do I1 e 3.6× da K2).

**Predição:** `rho_b_junc` (p90) 0.32 → 0.7–0.9; `V` ≤ 0.15; `n_bio` 167 → 190–210
(limite de C).
**Reverte se:** `a_mar_bio_med` < 2.5 ou `n_bio` > 220 — é o sinal do modo I1.

### J4 — fonte de nutriente (estrutural; só se J2/J3 não sustentarem além de t≈55)

Termo em `NutrientConsumption.post_loop`:

```
dc_n/dt += k_src * (1 - c_n)
```

Fundamentação [T2] Srinivasan 2019 (§3.0): swarming de *P. aeruginosa* é regime
**nutrient-rich** com `c ≈ c0` constante; um `c_n` que esgota monotonicamente é o
regime de **biofilme**, que a própria §3.0 identifica como o errado para este
organismo. Fisicamente: a placa de ágar é um reservatório 3D e a simulação é um corte
2D — a reposição vertical não está representada. Resolve também a lição #47 (janela
útil t ≲ 55 s) e estende a validação a t=100.

`k_src` inicial sugerido: 0.02 (τ = 50 s, mesma escala da janela útil).

---

## 5. Protocolo de validação — idêntico em todo passo

```bash
# 1. rodar (t=50, ~25 min)
make run

# 2. arquivar IMEDIATAMENTE (evita o erro do §0)
tools/archive_run.sh J1

# 3. métrica primária do problema + figura comparativa
python tools/diag_juncao.py runs/C4 runs/J1 --t 48 --out runs/juncao_J1.png

# 4. critérios C1/C2 e guardrails do projeto (§2.5)
python tools/compare_runs.py runs/C4 runs/J1

# 5. leitura visual — árbitro final (§11)
python tools/compare_frames.py --t 48 runs/C4 runs/J1
python tools/plot_fields.py runs/J1
```

Passo 5 não é opcional: §11 proíbe declarar sucesso por métrica escalar. `V` mede a
profundidade do vale, não a morfologia — um `V` baixo com dendritos borrados é falha.

**Registro:** cada passo entra na §9 do CLAUDE.md com predição *ex-ante* × resultado
medido (§2.3). Discrepância > 2× é investigada antes do passo seguinte.

---

## 6. O que NÃO fazer — já medido

| proposta | por que não | evidência |
|---|---|---|
| `k_col` 0.03 → 0.3 | na junção `deficit = local − ρ_b ≤ 0` **e** gate `c_n` = 0 → contribuição zero; só age na frente, onde achata a crista | `runs/K2`: `a_mar` 3.15→0.13, `n_bio` 167→606, REPROVA em C2 |
| `BiomassDiffusion` (`D_b`) | único termo bilateral; o maior gradiente com gate de produto ativo é núcleo(1.0)↔vale(0.12) → o fluxo só pode drenar o núcleo | [main.py:166](../main.py#L166): `rho_b` 1.0→0.48, `n_pinned` 43→**0** |
| inóculo `p=4 → 2` com R maior | suaviza a inclinação → biomassa +2.5× → `cs` uniforme → gradiente some | [particles.py:34-38](../src/particles.py#L34-L38): `a_mar` 2.30→0.75, morfologia Circular |
| as três juntas | perde atribuição (§2.3/§10) | — |

---

## 7. Divergências código × documentação a corrigir na §7 (só doc, risco zero)

Levantadas em 2026-08-11; todas do mesmo tipo — o Pass mudou o código, a §9 registrou
a narrativa, a tabela da §7 nunca acompanhou.

| item | código | §7 | impacto medido |
|---|---|---|---|
| pinning químico | `c_n < 0.6` | `< 0.4` | **alto** — 90.4% vs 83.2% da biomassa congelada; fronteira r 2.2 vs 1.1 |
| `k_consume` | `0.0` | `1.0` | nulo — T1 já o substituiu (`7041a6f`); restaurar mudaria `cs` em −1.2% (junção), −2.9% (braço) |
| gate Marangoni | `[0.15, 0.5]` | `[0.05, 0.6]` | **nulo** — `\|∇ρ_b\|` tem mediana 0.685 e p10 0.221; gate médio 0.779 sob ambos, razão 1.00× |
| `c0` | `0.35` | `0.8` | `B = ρ0c0²/7` 5.2× menor que o documentado (K16d, nunca revertido) |
| perturbação azimutal | `R_θ = 0.30+0.06cos8θ` | `0.8·cos(8θ)` | formulações diferentes; `noise` confere |

Nota: o gate `[0.15, 0.5]` é exatamente o valor que a §9 lista como "K.1 (refutado)";
voltou em `d6378e3` (2026-04-22, "pass k.16") sem registro. O veredito "sem efeito
mensurável" de K.1 estava **correto pelo motivo certo** — os dois gates são
indistinguíveis sobre a distribuição real de `|∇ρ_b|` — e não por cegueira do `max()`,
como a lição registra.
