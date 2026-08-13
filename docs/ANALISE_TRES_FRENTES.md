# Análise das frentes estruturais

> **Estabelecido 2026-08-13.** Escrito depois de K3 reprovar (`docs/PLANO_K3_JUNCAO.md`) e
> de o censo de partículas (lição #57) mostrar que a colônia visível é 94% filler
> congelado. Não são calibrações — são defeitos de modelagem que nenhuma alavanca
> paramétrica alcança.
>
> **REVISADO no mesmo dia.** A leitura inicial tinha três frentes e uma taxonomia errada
> (tratava o filler como poluente). A medição do ágar invadido mostrou que são **quatro
> defeitos** e que o filler é fase passiva legítima. As Frentes 1–3 abaixo permanecem
> válidas como diagnóstico; a **Frente 0** e a seção de **acoplamento** as reordenam e
> corrigem duas afirmações — leia essas duas primeiro se estiver decidindo o que fazer.

Todas as medições vêm de `runs/C4` (baseline validado, t=0→50) e `runs/K3`, lidas
direto do HDF5. Referências da literatura conforme catalogadas em CLAUDE.md §3.0.

**Índice.** [Frente 0 — ágar invadido](#frente-0--o-ágar-invadido-a-população-que-polui-tudo)
(escrita por último, leia primeiro) · [Frente 1 — inóculo](#frente-1--o-inóculo-semeado-e-o-halo-de-ρ_b--0) ·
[Frente 2 — rastro](#frente-2--o-rastro-furado-velocidade-da-ponta-vs-cadência-do-wake) ·
[Frente 3 — `ρ_b`/`ρ_max`](#frente-3--ρ_b-nunca-chega-a-ρ_max) ·
[Acoplamento](#acoplamento--quatro-defeitos-não-três)

---

## Frente 1 — O inóculo semeado e o halo de `ρ_b = 0`

### 1.1 A problemática, quantitativa

A condição inicial em [particles.py:38-39](../src/particles.py#L38-L39) é

```
R_θ  = 0.30 + 0.06·cos(8θ)
ρ_b  = exp(−(dist/R_θ)⁴)
```

Medido em t=0:

| | valor |
|---|---:|
| partículas com `ρ_b > 0.1` | **152** |
| partículas com `ρ_b > 0.01` | 213 |
| partículas com `ρ_b == 0` **exato** | **67 596** de 68 121 (99.2%) |
| raio máximo com `ρ_b>0.1` | **0.431** |
| biomassa total `Σ ρ_b·V` | 0.2544 |

A 4ª potência faz um quase top-hat: `ρ_b` cai de 1.0 a 0.15 entre r=0.30 e 0.45 e chega
a **zero exato** em r≈0.7. Como `BiomassGrowth` é multiplicativo (`dρ_b = rate·ρ_b`), zero
é **estado absorvente**: as 67 596 partículas fora desse raio nunca entram na biologia.
Validado por rastreamento de identidade (lição #48): das que estão em zero em t=50,
**100% nasceram em zero e nenhuma saiu**.

E a modulação `cos(8θ)` faz com que essa fronteira de zero seja **lobada** — a colônia
nasce com 8 setores que já contêm biomassa e 8 que já não contêm.

### 1.2 A problemática, morfológica

A colônia cresce em área ~25× (R de 0.43 para 4.4) enquanto o número de portadoras vai
de 152 para 168 (+11%). A densidade areal de biomassa cai ~19×. **O que se vê como
dendrito é a cauda da gaussiana inicial sendo advectada, não colônia colonizando espaço**
(lição #41). O halo roxo entre núcleo e braços é a região que a colônia *ocupa
mecanicamente* mas onde `ρ_b` nunca foi diferente de zero.

### 1.3 Um resultado que muda a recomendação

Mediu-se o espectro azimutal do resultado final (FFT do histograma angular, r>1.5, t=50):

| população | modos dominantes |
|---|---|
| biomassa **real** (n=32) | **11, 12, 10** |
| todas com `ρ_b>0.1` (n=1942, majoritariamente filler) | 16, **8**, 12 |

**O modo semeado (m=8) não é o que aparece na biomassa viva.** E a contagem visual nos
frames é de ~18–20 dendritos, não 8.

Isso bate com a escala de seleção física. O alcance lateral do surfactante no ágar é
`L_D_ext = √(D_ext/λ_ext) = √(0.08/0.75) = 0.327`; dedos separados por ~2·L_D_ext dão

| R | N previsto = 2πR / (2·L_D_ext) |
|---:|---:|
| 1.0 | 9.6 |
| 1.5 | 14.4 |
| **2.0** | **19.2** |
| 3.0 | 28.9 |

**Conclusão: a física já seleciona o comprimento de onda (10–20 dedos em R≈1.5–2), e a
semeadura m=8 não é o que fixa a contagem.** Isso é evidência direta de que remover a
semeadura geométrica é seguro — o que ela faz hoje é sobretudo criar a fronteira lobada
de `ρ_b=0`, ou seja, o custo sem o benefício.

### 1.4 O que a literatura diz

- **[T3] Giverso, Verani, Ciarletta 2016** — análise de estabilidade linear com *dispersion
  curves* e número de onda característico. O padrão tem comprimento de onda **selecionado
  pela dinâmica**; expansão volumétrica dá k=1 (assimetria), quimiotática dá múltiplos
  dendritos simétricos. Nada é imposto na condição inicial.
- **[T1] Trinschek, John, Thiele 2018** — os 7–9 dedos do painel (b) emergem do par
  (`W`, `Γ_max`). A morfologia é resultado dos parâmetros, não de semeadura.
- **[T5] Bru et al. 2023** — PA14 forma camada de surfactante e tendrils; MPAO1, que não
  produz surfactante na superfície, **não forma tendrils**. O discriminante biológico é a
  produção de surfactante, não a geometria do inóculo.
- **Xavier et al. 2011** (citado em CLAUDE.md §3.2) — a expressão de `rhlAB` é regulada
  por repressão catabólica e é heterogênea entre células.

A literatura é convergente: **o dedo nasce de heterogeneidade química amplificada por uma
instabilidade, não de um inóculo pré-lobado.**

### 1.5 Soluções propostas

**S1.1 — inóculo circular + heterogeneidade na produção de surfactante (recomendada).**
`R_θ = R₀` constante; a modulação sai da geometria e vai para o campo `noise`, que
multiplica a produção de `cs`. Em vez de `noise = 1 + 0.6·sin(8θ) + 0.01·rand`, usar um
**campo aleatório espacialmente correlacionado** com comprimento de correlação
`ℓ ≈ L_D_ext = 0.33` e sem modo azimutal imposto — banda larga, deixando a instabilidade
escolher. Biologicamente é a expressão heterogênea de `rhlAB` ([T5], Xavier 2011); numericamente
é o ruído de banda larga que [T3] pressupõe.
*Custo:* zero (condição inicial). *Risco:* se a instabilidade for subcrítica, o resultado
vira disco circular — mas §1.3 mostra que ela **não** é subcrítica hoje.
*Critério de aceitação:* contagem de dendritos entre 10 e 20 e espectro azimutal sem pico
em m=8; `a_mar_front` dentro de 20% do baseline.

**S1.2 — não suavizar o perfil.** Registrado para não se repetir: alargar a gaussiana
(p=2, R=0.45) **cura a borda e mata o motor** — +2.5× biomassa → `cs` uniforme →
`a_mar` 2.30→0.75, morfologia Circular (lição #48, e o próprio comentário em
[particles.py:34-38](../src/particles.py#L34-L38)). O perfil agudo é o preço da
seletividade; o problema não é a nitidez, é a ausência de colonização (Frente 3).

**S1.3 — o halo não se resolve na condição inicial.** O inóculo com `ρ_b>0.1` vai a
r=0.431 e a colônia chega a R99=4.36: **99.0% da área que ela ocupa em t=48 era ágar em
t=0**. Semear mais largo não alcança essa região, e semear o domínio inteiro afoga o `cs`
(custo medido de 172% a 3484% na produção — lição #50). **A Frente 1 elimina a fronteira
lobada; quem elimina o halo é a conversão (defeito C), e quem faz a conversão pegar é o
crescimento (defeito D) — ver Frente 0 e a seção de acoplamento.**

---

## Frente 2 — O rastro furado: velocidade da ponta vs cadência do wake

### 2.1 A problemática, quantitativa

O `WAKE` ([main.py:900-995](../main.py#L900-L995)) roda a cada `WAKE_FREQ = 100`
iterações. Para cada partícula com `ρ_b>0.05` que acumulou deslocamento `≥ 1·dx` desde a
última deposição, ele deposita **um cluster** (1 partícula na posição *antiga* `x_dep`,
mais até 6 num anel de raio 0.75·dx) e zera o acumulador.

Medido em `runs/C4`, ao longo do run:

| t | dt | intervalo do wake (100 it) | desloc. p99 | desloc. máx |
|---:|---:|---:|---:|---:|
| 12.9 | 0.0210 | 2.10 s | 2.3 dx | 3.4 dx |
| 25.5 | 0.0195 | 1.95 s | 2.1 dx | **10.2 dx** |
| 41.2 | 0.0185 | 1.85 s | 2.2 dx | **9.5 dx** |
| 44.6 | 0.0170 | 1.70 s | 2.0 dx | **10.3 dx** |
| 50.0 | 0.0170 | 1.70 s | 2.2 dx | 5.0 dx |

**A ponta típica percorre 2.0–2.4 dx entre duas chamadas do wake; a mais rápida percorre
5–10 dx.** O cluster depositado cobre ~1.5 dx. O restante do rastro fica vazio até a
chamada seguinte — e é depositado **no lugar onde a ponta estava, não onde ela está**.

Consequência medida no baseline: `void_07` salta para **26.4%** em t≈4 e a área sem
partícula só volta a 7% depois de ~10 s de cicatrização por shifting.

### 2.2 A problemática, morfológica

O rastro é preenchido **depois**, por partículas congeladas (`is_filler=1`) inseridas na
posição antiga. Isso produz exatamente o que o censo (lição #57) mediu: os dendritos são
90 fantasmas para cada partícula viva, e 638 fantasmas WAKE estão em movimento contra 16
vivas. **O braço não é construído pela bactéria que avança; é remendado atrás dela.**

### 2.3 O que a literatura diz

Não é questão biológica, é de discretização. O enquadramento correto é
**[T6] Liu §3.3 / [T7] Violeau §3.4-3.6**: a representação SPH exige cobertura espacial
contínua; um vazio na região que o fluido ocupa é erro de consistência, não estética
(§2.5 do CLAUDE.md). E **[T8] Lind et al. 2012** trata o caso de suporte truncado — mas
por redistribuição, que já está em uso (`ParticleShift`) e é justamente o mecanismo que
cicatriza tarde.

### 2.4 Soluções propostas

**S2.1 — depositar AO LONGO do segmento, não no ponto (recomendada).** Quando a partícula
acumula deslocamento `d` desde `x_dep`, depositar `n = round(d/dx)` partículas igualmente
espaçadas no segmento `x_dep → x`, em vez de uma na origem. A trajetória entre chamadas é
quase retilínea (ponta balística sob arrasto), então a interpolação linear é fiel.
*Custo:* multiplica a deposição por ~2 na p99 e ~5–10 nas mais rápidas — que são poucas.
O orçamento comporta: `WAKE_MASS_BUDGET=0.12` equivale a ~8200 partículas e hoje só
2379 são usadas (**29% do teto**).
*Critério:* `void_07` no pico (t≈4) abaixo de 15%, e razão fantasma/viva nos dendritos
caindo abaixo de 50.

**S2.2 — aumentar a cadência.** `WAKE_FREQ 100 → 20` põe o intervalo em ~0.35 s e o
deslocamento típico em ~0.45 dx, dentro de um cluster. *Custo:* 5× mais construções de
`cKDTree` (≈50 ms cada, 150 no run inteiro — desprezível contra 2378 s). *Limitação:* não
resolve a cauda rápida (2–3 dx por chamada nas pontas de `max_v`).

**S2.3 — as duas juntas.** S2.1 é robusta a qualquer velocidade; S2.2 reduz a extrapolação
que S2.1 precisa fazer. Recomendo aplicar S2.1 primeiro (alavanca única) e medir.

**Ressalva de escopo:** isto melhora a *cobertura* do rastro, não o fato de ele ser
preenchido por partícula inerte. Enquanto o depositado for `is_filler=1`, o braço continua
mecanicamente presente e biologicamente vazio. A conversão do rastro em biomassa viva é
Frente 3.

---

## Frente 3 — `ρ_b` nunca chega a `ρ_max`

### 3.1 A problemática, quantitativa

Contagem de partículas **vivas** por faixa de `ρ_b`, ao longo de todo o run:

| t | `max ρ_b` | `n ≥ 0.8` | `n ∈ [0.5,0.8)` | `n ∈ [0.1,0.5)` | biomassa |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 1.0000 | **43** | 42 | 67 | 0.2544 |
| 12.9 | 1.0000 | **43** | 42 | 73 | 0.2607 |
| 25.5 | 1.0000 | **43** | 43 | 77 | 0.2914 |
| 37.5 | 1.0000 | **43** | 51 | 73 | 0.3120 |
| 50.0 | 1.0000 | **43** | 63 | 62 | 0.3505 |

**`n(ρ_b ≥ 0.8) = 43` do primeiro ao último instante** — e são exatamente as 43 que
*nasceram* acima de 0.8. Em 50 s, **nenhuma partícula atravessou esse limiar**. O mesmo
valor aparece em K3, e é a coluna `n_pinned` do log.

Três causas independentes:

**(a) Teto por construção.** [equations.py:24-28](../src/equations.py#L24-L28) só permite
crescimento se `ρ_b < 0.8`. Como `rho_max = 1.0`, **o valor nominal é inalcançável por
crescimento** — o teto real da dinâmica é 0.8.

**(b) Escala de tempo 2–6× maior que a janela.** Solução logística de `ρ_b=0.3` a 0.8:

| `c_n_factor` | `r_eff` | tempo |
|---:|---:|---:|
| 1.0 | 0.0200 | **112 s** (2.2× a janela) |
| 0.7 | 0.0140 | **160 s** (3.2×) |
| 0.4 | 0.0080 | **279 s** (5.6×) |

**(c) O termo logístico nunca satura.** Com a população ativa em `ρ_b ~ 0.1–0.5`, o fator
`(1 − ρ_b/ρ_max)` vale 0.5–0.9 o tempo todo: a equação opera em regime **exponencial**, e
`ρ_max` é um parâmetro que não exerce função alguma.

### 3.2 A problemática, morfológica

A colônia é um núcleo congelado de 43 partículas herdadas da condição inicial, cercado de
material que nunca amadurece. Não há a estrutura de três zonas que a biologia exige
(núcleo denso maduro → swarmers → baias). O que existe é: 43 partículas em `ρ_max`,
~125 vivas entre 0.1 e 0.8, e 2739 fantasmas carregando 91% do campo `ρ_b`.

### 3.3 O que a literatura diz

**Não existe valor de `ρ_b` na literatura para comparar** — e isso é um ponto do próprio
modelo, não uma lacuna bibliográfica. `ρ_b` é campo **normalizado** (0 a `ρ_max`,
capacidade de suporte do logístico); [T2] usa fração volumétrica `φ`, [T1] usa altura de
filme e `Γ`. Registrado na lição #53.

O que a literatura **constrange**:

- **[T2] Srinivasan, Kaplan, Mahadevan 2019** — swarming é regime **nutrient-rich**
  (`c ≈ c₀` constante), com duas fases: ativa (células) e passiva (fluido/matriz). O corpo
  maduro está na capacidade de suporte; o gradiente vive na frente que avança. **O freio
  por nutriente que o modelo usa hoje (série M-B) é física de biofilme, não de swarm.**
- **[T1] Trinschek 2018** — crescimento bioativo com saturação; o fingering ocorre com a
  biomassa **em/perto da saturação** no corpo da colônia.

Ou seja: a literatura diz que o corpo deveria estar em `ρ_max` e a ação deveria estar na
borda. O modelo produz o oposto — corpo em 0.1–0.5 e 43 pontos congelados em 1.0.

### 3.4 Por que as tentativas anteriores falharam

| tentativa | resultado | lição |
|---|---|---|
| `r_growth` 0.02 → 0.04 / 0.08 | colônia colapsa (R 4.92→3.46), `a_mar` 3.37→2.47 | #43 |
| `k_col` (colonização) 0.03–0.3 | preenche mas afoga `cs`, motor cai 49–96% | #48, #49 |
| K3 (alvo-doador + gate de `cs`) | zeros eliminados na junção, mas `contrast_cs` 12.9→10.2 e 59 componentes | K3 |
| `k_src` (fonte de nutriente) | `c_n` nunca esgota, biomassa 2.15× — **portadoras continuam 167→178** | #51 |
| fluxo quimiotático de `ρ_b` | dilui: banda flagelar 2041→34 | #54 |

O padrão: toda rota que **cria ou move** biomassa esbarra na produção de `cs` por
PRESENÇA (Hill `qs` satura em `ρ_b≈0.1`) ou dilui abaixo do quórum. A série Y mostrou que
sob produção **linear** o preenchimento deixa de matar o motor (lição #55) — mas a faixa
dinâmica de `cs` explode sem o teto (#56).

### 3.5 Soluções propostas

**S3.1 — remover o teto artificial (pré-requisito trivial).** Trocar o gate
`ρ_b < 0.8` de `BiomassGrowth` por `ρ_b < ρ_max`; o próprio fator `(1 − ρ_b/ρ_max)` já
zera a taxa na saturação. Custo zero. **Sozinho não resolve** — sem (b) nada chega lá.

**S3.2 — separar maturidade de densidade: campo `φ_m` (estrutural, recomendada).**
Hoje `ρ_b` acumula três funções incompatíveis: densidade de biomassa viva, chave da EOS
(`fade` liga em 0.1) e critério do pin (≥0.8). Por isso (i) o fantasma precisa herdar
`ρ_b` para ter coesão, contaminando o campo (lição #57), e (ii) qualquer aumento de `ρ_b`
congela a colônia via pin (armadilha K.17, medida em K.21).

Introduzir um campo de **maturidade/EPS** `φ_m` que:
- cresça com o tempo passado em alta densidade (`dφ_m/dt = k_m·ρ_b·(1−φ_m)`);
- **controle o pin** (`φ_m > φ_pin`), em vez de `ρ_b ≥ 0.8`;
- **ative a EOS** — o que permite ao filler carregar `φ_m` sem carregar `ρ_b`, devolvendo
  `ρ_b` ao papel de biomassa viva pura.

Ancoragem: é a fase **passiva** de [T2] (matriz/fluido) separada da fase **ativa**
(células) — exatamente a decomposição que o paper propõe e que o modelo hoje colapsa numa
variável só. Resolve simultaneamente a contaminação de métricas (#57), o painel enganoso,
o gate espúrio da Marangoni e a armadilha do pin.
*Custo:* um campo novo, um termo de evolução, e revisão do pin e da EOS. É a maior das
três frentes.

**S3.3 — regime nutrient-rich + produção proporcional ao conteúdo.** Com `φ_m` no lugar do
pin, `r_growth` pode subir sem congelar a colônia — mas precisa de (i) `k_src > 0` para o
nutriente não esgotar (validado em N1: `min c_n` nunca abaixo de 0.376, biomassa 2.15×) e
(ii) produção de `cs` linear em `ρ_b` **com** o teto mantido, para que encher não afogue
(lição #55) nem exploda a faixa dinâmica (#56).
*Alvo quantitativo:* `r_eff` tal que 0.3→0.8 leve ~20 s, i.e. `r_eff ≈ 0.11` →
`r_growth ≈ 0.15` com `c_n_factor ≈ 0.7`.

---

## Frente 0 — O ágar invadido: a população que polui tudo

> Esta seção foi escrita depois das outras (2026-08-13) e **corrige** o enquadramento
> delas. Ela vem primeiro porque é o defeito de que os outros três dependem.

### 0.1 Taxonomia correta das populações

Havia uma confusão de conceito nas versões anteriores deste documento: filler foi tratado
como "vazio biótico". **Está errado** — o filler carrega as propriedades herdadas da mãe e
é justificável como a **fase passiva** de [T2] (matriz/EPS). Ele não é o poluente.

Contando dentro do **envelope** da colônia (a menos de 1.5h do corpo `ρ_b>0.1` — e **não**
do disco `R99`, que inclui as baias e infla a área por 3.3×):

| população | definição | C4 (t=50) | K3 (t=50) |
|---|---|---:|---:|
| **fase ativa** | `ρ_b>0.1`, `is_filler=0` | 168 — 1.4% | 242 — 2.8% |
| **fase passiva (matriz/EPS)** | `is_filler=1` (wake + insert) | 2739 — **23.1%** | 2283 — 26.0% |
| limbo sub-quórum | `0 < ρ_b ≤ 0.1`, `is_filler=0` | 357 — 3.0% | 1675 — 19.1% |
| **ágar invadido** | `ρ_b = 0`, `is_filler=0` | **8608 — 72.5%** | 4565 — 52.1% |

**O filler é 23% e é legítimo. O poluente é o ágar invadido: 72.5% do corpo da colônia no
baseline** — meio externo que a colônia engoliu sem converter, e que entra em toda métrica
normalizada por "colônia" e em toda leitura morfológica.

### 0.2 O mecanismo do rasgo, medido

**(a) O anel de ágar existe e é comprimido.** Em t=3.65, 203 partículas de ágar formam um
anel em r=0.774 com `ρ/ρ0 = 2.03`.

**(b) Ele é empurrado de forma desigual, com assinatura m=8.** Deslocamento do anel entre
t=3.6 e t=15.7, por setor azimutal (24 setores):

| | valor |
|---|---:|
| mediana global | 3.17 dx |
| setor mínimo | 1.08 dx |
| setor máximo | 5.31 dx |
| razão máx/mín | **4.9×** |
| modo dominante do espectro azimutal | **m = 8** |

É a impressão digital do inóculo lobado: os 8 lóbulos furam o anel em 8 pontos e mal o
tocam nos 8 setores entre eles.

**(c) A frente ultrapassa o anel.** No mesmo intervalo o raio da colônia vai de 0.75 a
1.47 (**+0.72**) e a mediana do anel vai de 0.774 a 0.94 (**+0.17**) — a frente avança
**3.6×** mais rápido. E mesmo no setor *mais* empurrado (5.31 dx = 0.286), ainda avança
**2.5×** mais rápido.

**(d) O que não é empurrado é engolido.** Partículas de ágar dentro do raio da colônia:
61 (t=3.6) → 416 → 1093 → **1787** (t=15.7). **43% delas nunca se moveram**
(deslocamento < 0.5 dx); deslocamento mediano das que estão dentro: 0.83 dx — menos de um
espaçamento de rede.

**Leitura:** os lóbulos decidem **onde** o anel rasga; quem determina **que** ele rasga é
outra coisa — §0.3.

### 0.3 A EOS do ágar é alavanca de primeira ordem

[equations.py:655-657](../src/equations.py#L655-L657) zera `fade_rep` **e** `fade_att`
para `ρ_b < 0.1`. Consequência medida em t=3.65: o ágar tem `|p| = 0` **exato** e
`|v|` mediano = 6.5e-76 — não é fluido, é poeira estática que só serve de suporte de
kernel.

**Um meio sem pressão não transmite empurrão.** A colônia comprime a única camada em
contato direto de kernel — daí o `ρ/ρ0 = 2.03` — e essa camada não repassa a compressão
para a seguinte. Sem onda de compressão, **qualquer** frente acaba ultrapassando o meio; a
geometria do inóculo só escolhe o azimute onde isso acontece primeiro.

Isso é fato estrutural, não calibração: é o mesmo `p = 0` da lição #40 (vazio aberto no
ágar nunca cicatriza) visto pelo outro lado — o ágar não fecha buraco **e** não sai da
frente.

### 0.4 As duas arquiteturas coerentes — e por que o modelo não é nenhuma

Se a frente tem que deixar de ser porosa, há exatamente dois caminhos autoconsistentes:

| | **A — Deslocamento** | **B — Incorporação** |
|---|---|---|
| o que faz | a colônia empurra o meio para fora | a colônia converte o meio que alcança |
| exige | EOS no ágar **+** inserção massiva | termo **aditivo** em `ρ_b` + maturação |
| partículas a inserir no run | **~20 400** (área 59.1 / dx²) | 0 |
| massa | **+30%** da massa do domínio | **zero** |
| contra o orçamento atual (`WAKE_MASS_BUDGET=0.12` = 8 190 part.) | **2.5×** | — |
| cadência necessária | 191 part/s em t=15.7 (≈11× o wake atual em R=4) | — |
| ancoragem na literatura | nenhuma referência descreve isso | [T1] filme fino espalhando sobre substrato; [T2] influxo osmótico `V₀` puxa fluido **do ágar para dentro da colônia** |

**Sobre a rota A — inserir para empurrar não desloca, acumula.** Uma inserida com `ρ_b>0.1`
herdado tem pressão e empurra o vizinho, mas o vizinho tem `p=0` e não repassa: cada
inserção desloca só o que está dentro de um raio de kernel, e o material se empilha ali.
O anel iria de `ρ/ρ0 = 2.03` para 3–4×, produzindo uma casca de matéria sem pressão a
várias vezes a densidade nominal — nem ágar nem colônia. Levada até o fim, a rota A
**exige** dar pressão ao ágar; e se o ágar tiver pressão, a onda de compressão se propaga
sozinha e a inserção massiva deixa de ser necessária.

**O modelo hoje não é A nem B**: engole e deixa inerte. É a pior combinação — paga o custo
geométrico do engolimento sem colher nem o deslocamento nem a conversão. Daí os 72.5%.

**A literatura aponta para B.** Em [T1] o domínio é um filme sobre substrato: a colônia
espalha **sobre** o ágar, não o empurra de lado. Em [T2] o influxo osmótico traz fluido
**do ágar para dentro** da colônia. Nas duas, o meio é **incorporado**. Corolário
desconfortável e importante: **o ágar engolido não é um defeito a evitar — é a
matéria-prima da colônia.** O defeito é ele nunca ser incorporado.

**Ressalva sobre B:** converter ágar em biomassa não é gratuito biologicamente — célula
nova consome nutriente. O consumo existe (`OxigenConsumption`), mas se a conversão ficar
vigorosa o balanço de `c_n` vira restrição de primeira ordem e a fonte `k_src` (validada em
N1: `min c_n` nunca abaixo de 0.376) deixa de ser opcional.

---

## Acoplamento — quatro defeitos, não três

A versão anterior desta seção dizia que "a Frente 3 é a única que transforma ágar em
colônia". **Está errado.** `BiomassGrowth` é **multiplicativo** (`dρ_b = rate·ρ_b`):
corrigir o teto de 0.8 e acelerar `r_growth` deixa as portadoras existentes mais densas e
não tira ninguém do zero. Converter exige termo **aditivo** — que é a
`BiomassColonization` (o Pass K3), e não a Frente 3.

| | defeito | mecanismo necessário | estado |
|---|---|---|---|
| **A** | inóculo lobado → rasgo em m=8 | condição inicial (Frente 1) | proposto |
| **B** | rastro furado → vazio geométrico | deposição ao longo do rastro (Frente 2) | proposto |
| **C** | ágar engolido não é convertido | termo **aditivo** (colonização) | K3 rodado |
| **D** | convertido não amadurece → `ρ_b` nunca chega a `ρ_max` | crescimento (Frente 3) | proposto |

### C e D são um encanamento — e é isso que explica a falha do K3

O K3 **converteu**: ágar invadido 72.5% → 52.1%, ou seja 16–20 pontos percentuais. Mas
desses, ~92% pararam no limbo sub-quórum (3.0% → 19.1%) e só ~1.4 ponto virou fase ativa
(1.4% → 2.8%).

A causa é aritmética. Uma recrutada nasce em `ρ_b ≈ 0.15`; para chegar a 0.5 (swarmer de
verdade), a solução logística dá:

| `r_eff` | tempo |
|---:|---:|
| 0.014 (hoje: `r_growth=0.02`, `c_n_f≈0.7`) | **124 s** |
| 0.11 (alvo da Frente 3) | **16 s** |

Numa janela de 50 s o primeiro nunca chega: a partícula fica em (0, 0.1), **mecanicamente
invisível e quimicamente ativa**, que é o que afogou o `cs` e derrubou `contrast_cs` para
10.2. **A Frente 3 não converte, mas é ela que faz a conversão pegar.** Sem D, qualquer
mecanismo aditivo deposita no limbo.

### O que cada frente alcança, e o que não alcança

- **A (inóculo)** muda o *padrão* do engolimento, não a *quantidade*: uniformizar o
  empurrão levaria o anel a r≈1.0 uniforme em vez de 0.80–1.07 lobado, e o vão entre anel
  e frente cairia de 0.53 para ~0.47. Marginal. Também não alcança o halo: o inóculo com
  `ρ_b>0.1` vai a r=0.431 e a colônia chega a R99=4.36, ou seja **99.0% da área final era
  ágar em t=0** — qualquer coisa feita no perfil inicial redistribui biomassa dentro de 1%
  da área.
- **B (rastro)** é vazio geométrico de verdade — a partícula saiu e não há ninguém para
  converter. Alvo pequeno e localizado, ao contrário dos 20 400 do domínio inteiro. Mas
  enquanto o depositado for `is_filler=1`, ele engrossa a fase passiva (já 23%) sem
  aumentar a ativa.
- **C (conversão)** é a única que transforma meio em colônia, e sai de graça em massa.
- **D (crescimento)** é o gargalo de C **e** o que dá sentido a `ρ_max`.

### Ordem revisada

```
D (crescimento/φ_m)  →  destrava C
C (conversão)        →  elimina o ágar invadido
B (rastro)           →  fecha o vazio geométrico residual
A (inóculo)          →  tira a assinatura m=8, por último
```

`D` primeiro porque sem ela `C` deposita no limbo — foi exatamente o que K3 mediu. E
dentro de `D`, `φ_m` primeiro: sem separar maturidade de densidade, subir `r_growth` faz
as partículas cruzarem 0.8 e o pin congela a colônia (armadilha K.17, medida em K.21).

`A` por último e não por primeiro, ao contrário do que a versão anterior deste documento
recomendava: os dendritos de hoje são a cauda semeada sendo advectada (lição #41), então
remover a semente **antes** de haver conversão capaz de amplificar ruído de banda larga
arrisca produzir colônia circular — o modo de falha clássico do projeto.

**O que NÃO fazer**, já refutado por medição: suavizar o perfil do inóculo (#48), subir
`r_growth` sem fonte de nutriente (#43), colonizar sem produção linear (#48/#49/K3),
transportar `ρ_b` por fluxo (#54), inserir partícula biologicamente ativa na frente (#36),
e — novo — **inserção massiva para empurrar o ágar**, que custa +30% de massa, produz
acumulação em vez de deslocamento enquanto `p=0`, e compra a arquitetura que a literatura
não descreve.
