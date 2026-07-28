# Consultando propriedades onde não há partículas

**Pergunta do cliente:** *"No último frame do movie, os braços do swarm são regiões
que literalmente não têm partícula nenhuma. Se eu quiser saber as propriedades
(surfactante, biomassa, nutriente) dessas regiões, como isso deve ser feito?"*

Este documento explica **a problemática** (os braços são vazios reais no campo de
partículas) e **a solução** (como extrair um valor honesto de um ponto qualquer,
incluindo dentro de um vazio, sem inventar dado nem corromper o solver).

---

## 1. A problemática — os braços são vazios reais

O SPH não guarda campos numa grade; guarda **partículas**, e cada uma carrega `cs`
(surfactante), `rho_b` (biomassa), `c_n` (nutriente). O "campo" num ponto `x` só existe
como interpolação ponderada pelos vizinhos:

```
A(x) ≈ Σ_j (m_j/ρ_j) · A_j · W(x − x_j, h)
```

Se o ponto não tem vizinhos próximos, a soma perde suporte e o valor fica indefinido.

Plotando **todas as partículas como um pontinho** (branco = vazio) do frame real
(`main_02273.hdf5`, N=35 252):

- O **ágar** é uma rede densa de partículas preenchendo o fundo.
- Os **braços dendríticos são canais brancos — literalmente sem partícula no interior.**
  As ~462 partículas de biomassa ficam nas **bordas** desses canais; o miolo é vazio.

Medindo o tamanho desse vazio (distância ao vizinho mais próximo, numa grade fina):

| Métrica | Valor medido |
|---|---|
| Domínio com vizinho mais próximo além de `1·dx` | **7,1%** |
| Domínio com vizinho mais próximo além de `1,5·dx` | **2,9%** |
| No **centro de um braço**: partículas dentro de `h` | **0** |
| No centro de um braço: partículas dentro de `2h` (só as paredes) | ~32 |
| **`σ_a` no centro do braço** | **≈ 0,03** |

> **Correção de um erro meu na 1ª versão deste documento.** Eu havia medido cobertura
> com raio `2h ≈ 0,19`, que é **maior que a largura do braço**. Assim um ponto no centro
> do braço "encontrava" as partículas das paredes e eu concluí, erradamente, que o
> domínio estava todo preenchido. **Está errado: os braços são vazios reais.** O ponto
> tem 0 partículas dentro de `h` e ~32 apenas na parede, a `~2h` de distância — onde o
> peso do kernel já é minúsculo. Por isso `σ_a ≈ 0,03`, não ~1.

### Por que esses vazios existem (e não são preenchidos)

A frente de biomassa avança para fora (Marangoni + flagelar) e, com o núcleo pinado e
**sem inserção de partículas na frontier**, o braço estica mais rápido do que o
suprimento de partículas — abrindo um tubo evacuado. A inserção que existe no baseline
(`use_insert=True`) tem `INSERT_RHO_B_MIN=0.5`: preenche só o vazio **estrutural** do
núcleo (`rho_b>0.5`), e **exclui de propósito a frontier dos braços** (`rho_b<0.5`).
Isso é intencional — preencher os braços quebra a física (ver §4).

---

## 2. `σ_a` é o juiz — e ele acerta no vazio

A grandeza que mede "quanta partícula o kernel enxerga em `x`" é a **partição da
unidade discreta**, já computada a cada passo pela equação
[`KernelSum`](src/equations.py#L94):

```
σ_a(x) = Σ_j (m_j/ρ_j) · W(x − x_j, h)
```

O ponto forte de `σ_a`: ele **pondera pelo kernel**, então não se deixa enganar pelas
32 partículas distantes na parede do braço — ele lê **≈ 0,03** no miolo do vazio,
exatamente o que se espera. Validação em [test_kernel_sum.py](test_kernel_sum.py):

| Cenário | `σ_a` | Significado |
|---|:---:|---|
| Suporte completo (bulk) | `1,00` | Dado plenamente confiável |
| Superfície livre (borda) | `≈ 0,69` | Kernel truncado, mas dado válido |
| Partícula isolada (auto-termo) | `≈ 0,14` | Piso — não chega a 0 |
| **Vazio (nenhum vizinho próximo)** | `→ 0` | **Sem material → sem dado** |

Consultando três tipos de ponto no frame real:

| Ponto | `σ_a` | vizinhos em `h` / `2h` | Diagnóstico |
|---|:---:|:---:|---|
| **Centro de braço** `(0.98, -0.02)` | **0,04** | 0 / 34 | **Vazio real — não há material, não há dado** |
| Baía de ágar entre braços `(1.5, 0.3)` | 0,66 | — / 38 | Ágar preenchido: **tem dado** (`rho_b=0`, `cs≈0,07`) |
| Junção do núcleo `(0.3, 0.0)` | 0,79 | — / 20 | Esqueleto denso: **tem dado** (`rho_b=0,57`) |

Ou seja: **nem tudo que parece vazio é vazio.** As baías entre os braços são ágar
preenchido (`σ_a` alto → há dado real: "ágar virgem"). Já o **miolo dos braços** é vazio
de verdade (`σ_a ≈ 0,04` → não há material a reportar). `σ_a` distingue os dois
automaticamente — é essa a grandeza que a consulta deve usar.

---

## 3. A solução

### 3.1 Onde há suporte (`σ_a ≥ σ_min`): interpolação de Shepard

Usa-se a soma SPH **normalizada pelo próprio `σ_a`** (Price 2007), o que restaura a
consistência de ordem zero e entrega valor válido mesmo no esqueleto ralo:

```
                Σ_j (m_j/ρ_j) · A_j · W(x − x_j, h)
A_Shepard(x) = ─────────────────────────────────────
                          σ_a(x)
```

Já implementado em [reconstruct_field (plots/plot.py:71)](plots/plot.py#L71) para os
painéis físicos.

### 3.2 Onde não há suporte (`σ_a < σ_min`): reporte "sem cobertura" — **não invente**

Este é o ponto crítico para o miolo do braço, onde `σ_a ≈ 0,03`. **Aqui não se pode
usar Shepard.** Dividir por `σ_a = 0,03` **amplifica ~30×** o valor das partículas da
parede e o projeta através do vazio — um número **fabricado**, não uma medição (não há
material ali para medir). A regra correta usa `σ_a` como máscara:

```python
if sigma_a(x) < SIGMA_MIN:      # ex.: 0.25
    return float('nan')          # "sem cobertura" — não há material, não há dado
else:
    return num / sigma_a(x)      # Shepard, válido
```

Um `σ_min ≈ 0,25` é defensável: fica acima do `σ_a ≈ 0,03` dos vazios e do piso de
auto-suporte (0,14) de uma partícula isolada, e abaixo das baías de ágar reais (0,66).
A resposta física honesta dentro de um braço é: **"este ponto é um vazio — não há
biomassa nem surfactante *aqui*; o material está nas paredes do braço."**

O gap atual do código: `reconstruct_field` faz `num / (den + 1e-12)`
([plot.py:80](plots/plot.py#L80)) — no vazio isso devolve `≈ 0` **silenciosamente**,
que se confunde com "medi zero". Bom para renderizar imagem; **enganoso para entregar
dado ao cliente.** A máscara por `σ_a` corrige isso.

### 3.3 Se o cliente quer o gradiente (∇cs), não o valor

Para gradiente numa zona rala (avaliar Marangoni), Shepard não basta — usa-se **Kernel
Gradient Correction / Bonet-Lok** ([`KernelGradientCorrection`](src/equations.py#L110),
disponível via `use_kgc`; **está `False` neste baseline** — ligar se o gradiente for
requisito). Torna `∇cs` exato para campo linear mesmo com vizinhança incompleta, sem
mover nem criar partícula. Validação em [test_kgc_consistency.py](test_kgc_consistency.py).

---

## 4. O que **NÃO** fazer: preencher o vazio do braço com partículas no solver

A tentação óbvia — "criar partículas no braço para o cliente poder consultar" — está
**proibida** pelo histórico do projeto, porque destrói a física que gerou o dado:

- Filler **congelado** na frontier → **mata o fingering** (lição #31).
- Filler **ativo** na frontier → **reacende o runaway** de massa/cs (lição #34/#36).
- Refinar 7× (Pass N) na frontier → **colapsa o `dt`** 35–97× (lições #29/#38).

É exatamente por isso que a inserção do baseline exclui a frontier (`INSERT_RHO_B_MIN=0.5`,
lição #35). **O conjunto de partículas do solver é sagrado.** "Obter dado num vazio" é
problema de **consulta / pós-processamento**, nunca de mexer no solver.

---

## 5. Receita prática (resumo)

Para responder "quais as propriedades em `(x, y)`?":

1. Somar sobre os vizinhos dentro de `2h` (via `cKDTree`), acumulando
   `num_A = Σ (m_j/ρ_j) A_j W` e `σ_a = Σ (m_j/ρ_j) W` — mesmo kernel e parâmetros do
   solver (`CubicSpline(dim=2)`, `h = 1.8·dx`).
2. Se `σ_a < σ_min` (≈ 0,25) → retornar `NaN` / **"vazio: sem material a reportar"**.
   É o caso do miolo dos braços.
3. Senão → retornar `num_A / σ_a` para cada propriedade (`cs`, `rho_b`, `c_n`). É o caso
   das baías de ágar e do núcleo.
4. Para gradiente → aplicar a correção KGC em vez da soma bruta.

Reaproveita integralmente a `KernelSum` (o `σ_a` já é calculado a cada passo). O único
acréscimo indispensável é a **máscara honesta por `σ_a`**, que separa "vazio real"
(braço, `σ_a≈0,03`) de "ágar ralo mas válido" (baía, `σ_a≈0,66`).

---

## Referências

- **Price 2007**, *splash* (PASA 24:159) — reconstrução de campo SPH via Shepard.
- **Bonet & Lok 1999** (CMAME 180:97) — correção de gradiente (KGC), consistência de
  1ª ordem com distribuição irregular.
- **Liu & Liu 2003** [T6] §3.3 — consistência do operador SPH e deficiência de suporte.
- **Violeau 2012** [T7] §3.4–3.6 — partição da unidade `σ_a` como medida do erro.
- Testes locais: [test_kernel_sum.py](test_kernel_sum.py) (σ_a: →0 no vazio, 1 no bulk),
  [test_kgc_consistency.py](test_kgc_consistency.py) (gradiente corrigido).
- Contexto do projeto: `CLAUDE.md` §12 "Rotas para o vácuo"; lições #31, #34, #35, #38.
