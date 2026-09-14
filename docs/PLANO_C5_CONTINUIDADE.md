# Plano C5 — continuidade da expansão

> **Estabelecido 2026-09-01.** Escrito depois de a usuária apontar, na trajetória do E5
> ([plots/traj_E5.png](../plots/traj_E5.png)), que em t≈13–17 os braços já estão avançados **sem nada que os
> conecte ao centro** — e que isso não acontece num swarm real, onde a expansão é gradual
> e as bactérias estão unidas umas às outras.
>
> A medição confirmou e foi além do que a observação visual sugeria: **mais de 90% da
> colônia não tem caminho de volta ao próprio centro**, em qualquer instante depois de
> t≈8. Isso virou o critério **C5** do CLAUDE.md §2.2 e é anterior a todos os alvos
> morfológicos — uma colônia que não expande de forma contígua não é um swarm,
> independente do AR que ela exiba.

Medições em `runs/E5_hillK025` e `runs/C4` (t=0→50), lidas direto do HDF5. Literatura
conforme catalogada em CLAUDE.md §3.0.

**Índice.** [O critério](#o-criterio-c5) · [Diagnóstico](#diagnostico--a-colonia-nao-expande-ela-e-construida) ·
[Rota A](#rota-a--portadoras-via-crescimento) · [Rota B](#rota-b--pressao-no-agar) ·
[Rota C](#rota-c--wake-continuo-recomendada) · [O que já foi tentado](#o-que-ja-foi-tentado-e-por-que-nao-resolveu)

---

## O critério C5

Um swarm expande como população **contígua**: cada bactéria em raio R chegou lá dividindo
e empurrando pelos raios intermediários, junto das vizinhas. Não existe mecanismo pelo
qual material apareça em R=4 enquanto o corpo conexo termina em R=0.8. Confirmado nas
imagens de PA14 em ágar de swarming (`reference.jpg`; ver também a série A–F de Morris
et al. em várias concentrações de ágar): a frente é contínua com a colônia em **todos**
os estágios, nunca há anel de material destacado.

| | métrica | alvo |
|---|---|---|
| **C5a** | `R_conn / R99` — raio alcançado pela componente conexa que contém o centro do inóculo, sobre o raio da colônia | → 1.0 |
| **C5b** | `frac_conn` — fração do material nessa componente | → 100% |

**Baseline em t=50, escala de contato (1.05 dx):**

| ligação | C4 `C5a` / `C5b` | E5 `C5a` / `C5b` |
|---|---|---|
| 1.05 dx | 0.12 / 4.6% | 0.20 / 8.2% |
| 1.20 dx | 0.18 / 12.7% | 0.49 / 36.0% |
| 1.40 dx | 0.20 / 17.2% | 0.49 / 49.0% |
| 2.00 dx | 0.66 / 27.7% | 1.02 / 59.0% |
| 2.70 dx | 1.02 / 49.5% | 1.02 / 72.9% |

A varredura de escala é obrigatória (lição #68). Mesmo na escala mais frouxa, **metade
do material do C4 continua desligado do centro**.

**É critério de TRAJETÓRIA, não de estado final — é o que o separa do C4.** Uma rodada
pode terminar aceitável em C4 e violar C5 o tempo todo, que é o caso do baseline. A
assinatura é o **raio conexo estagnar**: no E5 fica em 0.66–0.82 de t=8 a t=50 enquanto
`R99` vai de 0.95 a 4.13. As duas métricas decaem monotonicamente de 1.0 / 100% e nunca
se recuperam.

Validar C5 exige `python tools/plot_trajetoria.py <run>` e olhar a série, nunca só o
último frame.

---

## Diagnóstico — a colônia não expande, ela é construída

Dos 1343 pontos de material nos braços (`r > 2`), **1312 (98%) foram inseridos ali**.
Apenas 31 são partículas originais, e essas viajaram 52 dx de mediana partindo de
`r`=0.34.

| mecanismo | partículas | natureza |
|---|---:|---|
| advecção (transporte real) | **31** | contígua, e desprezível |
| **wake** (rastro da ponta) | **2410** | criação pontual |
| insert (vácuo) | 698 | criação pontual |

### (A) A advecção não funciona porque o ágar não é empurrado

Deslocamento mediano das 20641 partículas de ágar dentro de `r`<4.5: **0.005 dx**; só
26% se moveram mais de 1 dx. A causa está na `BiomassEOS`: abaixo de `rho_b = 0.1`,
`fade_rep = fade_att = 0`, então o ágar tem pressão **exatamente zero** e não transmite
empurrão. É a lição #58 vista pelo eixo temporal — a colônia não desloca o meio, ela o
atravessa.

### (B) O crescimento não funciona porque não há portadoras

Vivas por dx de perímetro: **3.02 (t=0) → 0.57 (t=50)**. A população viva vai de 152 a
276 enquanto o perímetro vai de 50 para 482 dx — cerca de **uma célula viva a cada 2 dx
de frente**. E `BiomassGrowth` é multiplicativo, então zero é absorvente: só a
colonização (aditiva) recruta, e o gate `cs > 0.6·cs_max` atrasa a partida em ~20 s.

### (C) O wake agrava porque deposita em pulsos

Aglomerados de até 7 partículas a cada `WAKE_FREQ`=100 iterações, enquanto a ponta
percorre 2–7 dx no intervalo. Medido: as 2410 partículas de wake formam **612 componentes
na escala de contato, com mediana de 2 partículas cada**. O braço parece um braço e **é
um colar de grumos desconectados**.

### Conclusão

**C5 é violado por construção, não por calibração.** Nenhum ajuste de parâmetro corrige,
porque o mecanismo que constrói os braços é deposição, e deposição é descontínua por
natureza.

---

## Rota A — portadoras via crescimento

**BLOQUEADA por uma desigualdade.** Protocolo §3.3.6 com `sigma`=11.1, `HILL_K`=0.25,
`cs_max`=0.5, três formas testadas:

| forma | perfil de `cs` | veredito |
|---|---|---|
| **com teto** (atual) | mid-arm a **0.98** do teto | sem espaço para mais biomassa |
| **sem teto + sumidouro ∝ `rho_b`** | núcleo 0.13, rim 0.43 — `núcleo/rim` = **0.31** | **patológico** (§3.3.4: push para dentro) |
| **sem teto, sem sumidouro** | faixa dinâmica explode | Y3 mediu `a_mar` = 116 (lição #56) |

**Por que restaurar `k_consume` não resolve:** no denominador
`lam_eff + k·rho_b + P/cs_max`, o termo do teto vale **~14** contra `lam_eff` = 0.30 e
`k·rho_b` ≤ 3. **O teto já é o sumidouro dominante** — `k_consume`=3 leva o mid-arm de
0.98 para apenas 0.89 do teto.

**Por que o sumidouro proporcional a `rho_b` é estruturalmente errado:**
`cs_inf = P/(lam + k·rho_b)` dá **menos** `cs` onde há **mais** biomassa, invertendo o
perfil. A razão `núcleo/rim` fica em 0.31 para qualquer `k`, ou seja o pico migra para o
rim externo — a linha patológica da §3.3.4, em que a maior parte da colônia recebe push
para dentro.

Para o sumidouro importar seria preciso `sigma·qs ≤ cs_max·(lam + k·rho_b)`, o que dá
**`sigma` ≤ 1.0** contra 11.1 — e o E4 já mediu que `sigma`=1.93 mata o motor
(`a_mar` 3.71, `R99` 1.85).

**Desbloqueio requer resolver `cs_max` primeiro**, que é o problema aberto da lição #64 e
é anterior a esta rota.

---

## Rota B — pressão no ágar

**NÃO TESTADA NA FORMA CERTA.** A alavanca é `fade_rep > 0` com `fade_att = 0` abaixo de
`rho_b = 0.1` — ramo **só repulsivo**, que resiste à compressão sem puxar a colônia para
dentro.

A lição #66-C mediu que dar coesão a região sub-densa **contrai** a colônia (P3: `R99`
4.62 → 1.93), mas isso é o **ramo atrativo**, o errado para este fim. E
`runs/A2b_agarfade` testou `agar_fade=1.0` com os **dois** ramos, contra critérios
anteriores ao C5 (a nota da rodada diz "melhora numérica, zero efeito na desjunção",
medido pela ocupação em t≈29 — não por `C5a`/`C5b`).

**Predição se implementada:** o ágar passa a resistir e a ser empurrado; o anel comprimido
que a lição #58 mediu em `rho/rho0` = 2.03 com `|p|` = 0 passaria a empurrar de volta.
Risco conhecido: a lição sobre `AGAR_DRAG_RATIO` (`runs/A1`) mostra que pressão e arrasto
no ágar não são independentes — baixar o arrasto sem pressão levou `rho/rho0` a 16.7 e
`R99` a 3.06.

**Custo:** uma rodada. **Risco:** alto — muda a dinâmica de ~40000 partículas de ágar.

---

## Rota C — wake contínuo (RECOMENDADA)

**A mais barata, e ataca C5 diretamente.** Medido no E5:

| | |
|---|---|
| componentes do material a 1.05 dx | 663, **todos a ≤ 8 dx de um vizinho** |
| ligações necessárias (árvore geradora mínima) | 662 |
| comprimento das ligações | mediana **1.22 dx**, p90 2.30 dx, total 994 dx |
| ligações que cabem em ≤ 2 partículas | **74%** (491 em ≤1.5 dx, 610 em ≤2.5 dx) |
| **custo** | **~1467 partículas** (+61% do wake atual, **+2% de massa**) |

**Predição:** `C5b` de 8.2% para ~100%; `C5a` de 0.20 para ~1.0; `R99` e contagem de dedos
preservados (a deposição é ao longo do rastro que já existe, não uma nova frente).

### Diferença crucial para o B1, que foi reprovado

O **B1** (lição #65, `runs/_logs_reprovados/B1_REPROVADO`) depositava
`n = round(d/dx)` pontos ao longo do segmento **inteiro** a cada chamada, houvesse vão ou
não, e perdeu **25% de AR** (10.48 → 7.83) por engrossamento, com a razão fantasma/viva
subindo de 91.8 para 122.4.

A proposta aqui é **depositar só onde o rastro está quebrado** — ~1467 partículas
dirigidas contra a deposição cega do B1. E o B1 foi julgado contra critérios que **não
incluíam C5**; sob C5 ele ataca o defeito certo pelo método errado.

### Ponto de partida: E5, não C4

Medido o custo da própria rota C nos dois baselines:

| | material | % wake | componentes | custo MST | **custo relativo** |
|---|---:|---:|---:|---:|---:|
| C4 t=50 | 3028 | 82% | 742 | 1634 | **54%** |
| C4 t=95 | 5035 | 77% | 1190 | 2436 | 48% |
| **E5 t=50** | 3384 | 71% | **663** | **1467** | **43%** |
| **E11 t=100** | 9591 | 66% | 1299 | 2634 | **27%** |

1. **O defeito é menor no E5** — menos componentes (663 contra 742) com *mais* material.
2. **C4 congela em t≈55** (`mean_v` cai 25×, `min_c_n` → 0, lição #47). C5 é critério de
   trajetória; num baseline morto a melhora é inmensurável. O E11 em t=100 está
   **acelerando** (`mean_v` 0.00125 contra 0.00088 em t=50).
3. **A rota C adiciona material, que precisa de nutriente para amadurecer.** No C4 o
   nutriente acabou e o depositado ficaria inerte — é a mesma razão pela qual o gate da
   colonização estava fechado no C4 (lição #71-A: 36 elegíveis contra 712).
4. **O custo cai com o tempo no E5** — 43% em t=50, **27%** em t=100, com 89% das ligações
   em ≤1.5 dx. A colônia se preenche sozinha e sobra menos para a rota C.

O único eixo em que o C4 seria preferível é `frac(cs>0.45)` (9.7% contra 52.4%), que é
preocupação da rota A e não toca C5.

**Instante de referência morfológica: t=50.** Entre t=50 e t=100 o AR cai de 7.84 para
**4.16** — a biomassa dobra e os dedos engordam (lei (L), lição #67-L). O t=100 vale como
prova de que a janela existe e para medir C5, não como a figura da morfologia.

### Implementação

No bloco do wake em `main.py`, antes de depositar: medir o vão entre a posição atual da
ponta e o material já existente no rastro; depositar `ceil(vao/dx)` pontos **apenas se o
vão exceder a escala de contato** (1.05 dx). Um parâmetro novo, `WAKE_GAP_MAX`, limitando
quantos pontos por chamada para não mascarar um salto anômalo da ponta.

### Critérios de aceitação (pré-registrados, §2.3)

1. **`C5a` ≥ 0.8 e `C5b` ≥ 70%** a 1.05 dx em t=50 (baseline: 0.20 / 8.2%)
2. **`C5a` e `C5b` não decaem monotonicamente** na trajetória — é o ponto do critério
3. **AR ≥ 6.8** (baseline E5: 7.61; tolerância de 10%, contra os 25% que o B1 perdeu)
4. **dedos ≥ 22** (baseline: 24)
5. **massa ≤ 211** (baseline: 206.1; +2% previsto)
6. `a_pressure` mediana ≤ 3.0 e `a_mar` ≥ 10 (baseline: 2.58 e 12.31)

**Se falhar em 3 ou 4** (engrossamento), o modo provável é o mesmo do B1 e a resposta é
apertar `WAKE_GAP_MAX`. **Se falhar em 1 ou 2** com 3–6 passando, o vão não estava no
rastro e o diagnóstico precisa ser refeito.

---

## O que já foi tentado, e por que não resolveu

Cerca de quinze rodadas de preenchimento — séries **J**, **K3**, **N**, **P**, **V**,
**W** e **E7–E10** — atacaram o sintoma de **estado final** de um defeito de
**trajetória**. Promover material em t=50 não faz a colônia ter se expandido de forma
contígua, e por isso nenhuma variante resolveu:

| série | mecanismo | por que falhou |
|---|---|---|
| J, K3, N2 | colonização (`k_col`) | afoga `cs` ou deposita sub-quórum (lições #48–#53, #59) |
| P | promoção do limbo | coesão em região sub-densa **contrai** (lição #66-C) |
| X | fluxo quimiotático de `rho_b` | dilui: espalhar por 25× mais partículas divide o valor por 25 (lição #54) |
| E7–E9 | piso por enclausuramento | come a baia; e desloca os dedos por competir com a colonização (lição #72) |
| E10 | piso topológico | **reloca** os buracos: fecha 609 dx² e abre 333, só 37 coincidem (lição #72-J) |

A conclusão metodológica está na lição #72-K: o preenchimento **visual** pertence à
renderização (`tools/plot_piso.py --fill`), com perturbação zero. O que pertence ao
solver é o mecanismo de **expansão**, que é o objeto deste plano.

---

## Resultado do E11 (t=100) — encher a junção não resolve C5

> **Executado 2026-09-01**, `runs/E11_t100` — E5 sem piso, `total_sim_time`=100. Teste
> pré-registrado: se `C5a` saltasse para perto de 1.0 quando a banda [0.6, 1.8) enchesse,
> a desjunção seria o gargalo. **A predição de que não era está confirmada.**

**A junção fechou, e por crescimento próprio:**

| | t=50 | t=100 |
|---|---:|---:|
| `rho_b_p90` da banda (só vivas) | 0.100 | **0.194** |
| portadoras na banda | 151 | **444** |
| biomassa real total | 0.455 | 0.666 (+46%) |
| `c_n` nos braços | 0.552 | 0.758 |
| `a_mar` líquida | 12.3 | **13.65** |
| massa | 206.1 | 214.1 (+3.9%) |

**Mas C5 quase não se moveu**, medido pelo maior componente a 1.05 dx:

| | partículas no maior | alcance | `C5a` | `C5b` |
|---|---:|---|---:|---:|
| t=50 | 237 | r ≤ 0.82 | 0.21 | 7.5% |
| t=75 | 1228 | r ≤ 1.67 | 0.34 | 19.3% |
| **t=100** | **1684** | r ≤ **1.92** | **0.36** | **17.6%** |

O corpo conexo cresceu **7× em partículas** e mais que dobrou em alcance — e **parou em
r=1.92 com `R99` = 5.31**. Os braços além disso continuam grumos depositados, exatamente
como o diagnóstico previu. **A junção era uma das quebras, não a quebra.**

Ver [plots/E11_C5.png](../plots/E11_C5.png) e a trajetória completa em
[plots/traj_E11.png](../plots/traj_E11.png): o maior componente (verde) é um disco que termina no círculo
verde; todos os braços ficam cinza.

### Defeito novo — o núcleo pinado se descola

| | comp. do CENTRO | maior comp. | centro está nele? |
|---|---:|---:|:---:|
| t=50 | 237 | 237 | sim |
| t=75 | 1228 | 1228 | sim |
| **t=100** | **88** | 1684 (r 0.30 → 1.92) | **não** |

Entre t=75 e t=100 o núcleo (`rho_b ≥ 0.8`, pinado com `u=v=0` por K.17) fica parado
enquanto a casca em volta cresce e se move, abrindo um vão em **r ≈ 0.30**. É o vácuo
cinemático da lição #32, agora no centro em vez da junção.

**Não aparece em t=50** — só surge quando a colônia amadurece. Mais um caso de critério de
trajetória pegando o que o frame final esconde. Deve ser tratado junto da rota C, ou o
C5a fica limitado por um vão de 1 dx no centro por mais conectados que os braços fiquem.

### O que o E11 valida como subproduto

Motor vivo e **crescente** em t=100 (`a_mar` 13.65 contra 12.31 em t=50), nutriente
sustentado (`c_n` = 0.758 nos braços, contra 0.086 do C4 em t=89 — lição #47), massa
controlada. **A janela útil do E5 vai bem além de t=50**, o que não era verdade para o C4.

---

## C1 — rota C v1 (deposição no segmento): REPROVADA

> **Executado 2026-09-01**, `runs/C1_wakeseg`, branch `feat/continuidade`. `WAKE_SEG=True`,
> `WAKE_SEG_MAX=6`: o wake percorre o segmento entre `x_dep` e a posição atual da ponta,
> depositando a cada ~1 dx **onde o teste de proximidade permite**. O anel de adensamento
> (`n_extra`) fica só no ponto de origem — a diferença para o B1.

**Critérios pré-registrados: 4 de 6 passam, os 2 que importam falham.**

| | E5 | **C1** | veredito |
|---|---:|---:|---|
| **`C5a` @1.05dx** | 0.20 | **0.35** | **REPROVA** (≥ 0.8) |
| **`C5b` @1.05dx** | 8.2% | **11.4%** | **REPROVA** (≥ 70%) |
| AR | 7.61 | 6.97 | PASSA (≥ 6.8) |
| dedos | 24 | 22 | PASSA (≥ 22) |
| massa | 206.1 | 205.6 | PASSA (≤ 211) |
| `a_press` mediana | 2.58 | 2.65 | PASSA (≤ 3.0) |
| `a_mar` | 12.31 | 13.09 | PASSA (≥ 10) |

### Por que falhou

**Não adicionou material — redistribuiu.** 2185 partículas contra 2119 do E5 nas mesmas
28 chamadas (o teto `WAKE_MAX` travou em só 7% delas, e o orçamento de massa usou 2185 de
8174). O que aconteceu: os pontos do segmento consomem a mesma cota que iria para o anel
de adensamento, e o cálculo de déficit do `WAKE_MODE=2` reduz `n_extra` quando a
vizinhança já ficou densa.

**E tirou material dos braços.** Em `r > 2`: 1122 partículas contra **1343** do E5, com
`nn` mediano 0.79 dx contra 0.75 e a fragmentação inalterada (287 contra 293 componentes,
maior 32 contra 29). O ganho todo foi na **junção** — o maior componente vai de 278
partículas até r=0.82, para 356 até r=1.41.

**O alvo estava errado.** O segmento entre `x_dep` e a posição atual não é onde está a
quebra; preencher *dentro* de um segmento não emenda um segmento no **seguinte**. Com
`WAKE_FREQ=100` a ponta anda ~2 dx entre chamadas e cada chamada deposita uma corrente de
1–2 pontos: medido nos depósitos mais recentes, 40 componentes para 73 partículas, o maior
com 5. Elas nascem 2.5× mais próximas **entre si** (`nn` 0.83 dx contra 2.12 do E5) mas
**menos ancoradas** ao material preexistente (62% contra 82% a menos de 1.05 dx), porque
marcham na direção da ponta, entrando em ágar virgem.

### O que sobrevive

**A deposição dirigida é barata:** AR perdeu 8% contra os 25% do B1, e a massa ficou
igual. O princípio de só depositar onde há vão funciona — o alvo é que estava errado.

### Correções de método nesta rodada

1. **A predição de massa do plano (+2%, "massa ≤ 211") era do custo de ligar o estado
   final UMA vez** (árvore geradora mínima). O wake contínuo paga a cada chamada. O
   critério passou por acidente — a rota não adicionou massa por outro motivo.
2. **O primeiro plot comparava braços diferentes** — o zoom escolhia "o mais povoado" em
   cada run independentemente, e sugeria que o C1 tinha braços contínuos. Com a mesma
   janela angular (θ=35°) o C1 tem **198** partículas contra 256 do E5, `nn` 0.77 contra
   0.64: os braços estão mais **esparsos**. Ver `plots/C1_rotaC.png`.

### Próximo: cadência

A hipótese é que o defeito é de **frequência**, não de interpolação: com `WAKE_FREQ=100` a
ponta anda 2 dx entre chamadas e os depósitos saem em ilhas. Com `WAKE_FREQ=50` ela andaria
~1 dx e o rastro sairia contínuo por construção. Custo: 2× mais chamadas do wake, cada uma
O(N log N) numa KDTree.

---

## C2 — cadência (`WAKE_FREQ` 100 → 50): REPROVADA, mas com tendência

> **Executado 2026-09-01**, `runs/C2_freq50`. Hipótese: o defeito do C1 é de **frequência**,
> não de interpolação — com `WAKE_FREQ=100` a ponta anda ~2 dx entre chamadas e os depósitos
> saem em ilhas que não emendam.

| | E5 | C1 (freq 100) | **C2 (freq 50)** |
|---|---:|---:|---:|
| **maior cadeia do wake** | 99 | 215 | **467** |
| **maior componente** | 278 | 356 | **575** |
| **`C5b`** | 8.2% | 11.4% | **16.6%** |
| `C5a` | 0.20 | 0.35 | 0.33 |
| componentes @1.05 | 663 | 651 | 641 |
| material total | 3384 | 3134 | **3458** |
| AR / dedos | 7.61 / 24 | 6.97 / 22 | **6.93 / 24** |
| massa / `a_press` / `a_mar` | 206.1 / 2.58 / 12.31 | 205.6 / 2.65 / 13.09 | 206.2 / 2.65 / 12.01 |

**A hipótese está parcialmente certa.** A maior cadeia do wake mais que dobra a cada
divisão da frequência (99 → 215 → **467**, ~2.2× por passo) e o `C5b` dobra contra o E5.
A rota volta a **adicionar** material (3458 contra 3384) em vez de redistribuir, com AR e
dedos recuperados e massa inalterada. Custo real: wall time +21% (1339s → 1616s), não 2×.

**Mas reprova:** `C5a` = 0.33 contra 0.8; `C5b` = 16.6% contra 70%.

**Sinal de alerta — a contagem de componentes não cai** (663 → 651 → 641). Se as correntes
estivessem se fundindo, ela despencaria. Em vez disso a maior cresce enquanto o total fica:
cada chamada forma cadeia longa **e** ilhas novas na mesma taxa.

**As duas extrapolações discordam**, e é isso que o C3 (`WAKE_FREQ=25`) decide:

- pela **maior cadeia** (2.2× por divisão): freq 25 daria ~1000 e freq 12 daria ~2000 de
  2569 — quase tudo numa cadeia só
- pelo **`C5b`** (1.46× por divisão): chegar a 70% exigiria `WAKE_FREQ ≈ 4`, inviável

**Critério de decisão pré-registrado para o C3:** se `C5b` ≈ 30% e a maior cadeia ≈ 1000,
a tendência é boa e vale continuar. Se `C5b` ficar em ~20% com a cadeia crescendo, a rota C
**satura** e o caminho passa a ser a rota B (pressão no ágar).

---

## C3 e o fim da rota C — a série de frequência era RUÍDO

> **Executado 2026-09-02**, `runs/C3_freq25`. O critério pré-registrado do C2 pedia
> `C5b` ≈ 30% e cadeia ≈ 1000 para continuar. Deu `C5b` = 11.8% e cadeia 195 — a
> tendência **inverteu**. E a checagem seguinte mostrou que nem a tendência nem a
> inversão existiam.

### `C5b` varia ±7 pontos dentro de um mesmo run

| run | t=35 | t=40 | t=44 | t=47 | t=50 | média | **desvio** |
|---|---:|---:|---:|---:|---:|---:|---:|
| E5 | 10.6 | 14.6 | 15.5 | 11.5 | 8.2 | 12.1 | 2.6 |
| C1 (freq 100) | 26.2 | 30.3 | 26.1 | 17.3 | 11.4 | 22.3 | **6.9** |
| C2 (freq 50) | 7.6 | 28.6 | 16.9 | 13.7 | 16.6 | 16.7 | **6.8** |
| C3 (freq 25) | 11.5 | 17.9 | 29.2 | 24.9 | 11.8 | 19.1 | **7.0** |

As diferenças que eu interpretei entre configurações (11.4 → 16.6 → 11.8) são de ~5
pontos, **menores que o ruído de ~7**. Li uma tendência em três amostras ruidosas
(99 → 215 → 467 na cadeia do wake) e depois li uma inversão nas **mesmas** amostras — as
duas leituras eram infundadas.

**`SEED` fixo não resolve.** Ele garante a condição inicial; a morfologia diverge
caoticamente assim que os parâmetros mudam. É a mesma classe da lição #46: decidir a
ESTATÍSTICA antes de ranquear. Corrigido no CLAUDE.md §2.2 e em
`tools/plot_trajetoria.py`, que agora imprime média ± desvio sobre t ∈ [35, 50].

### O que a rota C entregou, medido corretamente

| | `C5b` médio (t=35–50) |
|---|---:|
| E5 (sem rota C) | **12.1** |
| rota C (média de C1/C2/C3) | **19.4** |

**+7 pontos, cerca de um desvio** — real mas modesto. E **`WAKE_FREQ` não tem efeito
atribuível** entre 25 e 100: os três valores são indistinguíveis.

### Veredito: rota C REPROVADA

Contra o alvo de `C5b` ≥ 70%, um ganho de 12 → 19 satura muito abaixo. O mecanismo de
deposição dirigida é barato (AR −8% contra os −25% do B1, massa inalterada) e vale
preservar, mas **não resolve C5**.

Isso é consistente com o diagnóstico da lição #73: o wake é um mecanismo de **construção**,
e refinar como ele constrói não o transforma em expansão. Enquanto o material chegar aos
braços por deposição — 98% em `r > 2` — a colônia não vai ser contígua.

**Próximo: rota B** (pressão no ágar, ramo só repulsivo abaixo de `rho_b = 0.1`), que ataca
a causa em vez do sintoma: com um meio que resiste, a colônia passa a **deslocar** em vez de
atravessar, e a advecção — hoje 31 partículas de 1343 — vira o mecanismo de transporte.

---

## B1 e B2 — rota B (pressão no ágar): REPROVADA

> **Executado 2026-09-02**, `runs/B1_agarfade` e `runs/B2_agar_press_drag`, isolados
> (`WAKE_SEG=False`, `WAKE_FREQ=100`). O mecanismo **já existia** no código: `agar_fade`
> faz exatamente o ramo só repulsivo (`fade_rep = agar_fade`, `fade_att = 0` para
> `rho_b < 0.1`), apenas desligado.
>
> Os dois testes anteriores dele eram inválidos para esta pergunta: **A2** rodou com o bug
> de `Group` da `BiomassEOS` (o ramo repulsivo nunca executava), e **A2b** mediu por
> ocupação em t≈29, antes do C5 existir.

| | E5 | B1 (pressão) | **B2 (pressão + arrasto 0.1)** |
|---|---:|---:|---:|
| `\|p\|` máx do ágar | 0.000 | **0.121** | 0.078 |
| ágar com `p > 0` | 0 | **21268** | — |
| `rho/rho0` do ágar máx | 5.86 | 3.63 | **3.12** |
| deslocamento do ágar p50 | 0.011 dx | 0.025 dx | **0.291 dx** |
| % do ágar que move > 1 dx | 29.6% | 30.2% | **42.2%** |
| **ADVECTADAS em r > 2** | **31** | **31** | **34** |
| **`C5b`** (t 35–50) | 11.3% ± 2.9 | 9.4% ± 1.1 | **2.6% ± 0.8** |
| AR / dedos | 7.61 / 24 | 7.03 / 25 | **3.13 / 15** |
| `R99` | 4.15 | 4.07 | 3.57 |
| `a_mar` | 12.31 | 11.73 | **4.02** |

**B1 — o mecanismo engajou e não adiantou.** O ágar passou a resistir (21268 partículas
com pressão, deixou de ser esmagado a 5.86) mas **não passou a ser empurrado**: advecção
idêntica, `C5b` indistinguível do baseline. Dar resistência fez o meio parar de ser
comprimido, não passar a escoar.

**B2 — mover o meio destrói a colônia.** Com o arrasto em 0.1 o ágar se move 26× mais, e
se move **junto** com a colônia em vez de sair da frente: AR 7.61 → **3.13**, dedos 24 →
15, motor a um terço (`a_mar` 4.02), `C5b` a 2.6%. É o modo de falha do `runs/A1`
reproduzido, agora com pressão.

### O achado: a advecção é invariante

**31 → 31 → 34.** Três regimes radicalmente diferentes do meio — rígido sem pressão,
rígido com pressão, móvel com pressão — e o número de partículas que a colônia
efetivamente **transporta** além de `r = 2` não muda.

**Logo o gargalo não é o meio.** Não é que a colônia não consiga empurrar o ágar; é que ela
não tem força motriz para carregar um corpo contíguo para fora, faça o meio o que fizer.
Com ~276 portadoras vivas e `a_mar` ~12, o que ela consegue é mover as pontas e depositar
atrás.

---

## Convergência: as três rotas apontam para `cs_max`

| rota | ataca | status |
|---|---|---|
| **C** — wake contínuo | como o material é **depositado** | **REPROVADA** — `C5b` 12 → 19, satura (C1/C2/C3) |
| **B** — pressão no ágar | o **meio** que deveria ser deslocado | **REPROVADA** — advecção invariante (B1/B2) |
| **A** — portadoras via crescimento | a **força motriz** | bloqueada pelo teto de `cs` |

C e B estão refutadas empiricamente. A que resta é justamente a que ataca o que a invariância
da advecção identificou como limitante — **a força motriz** — e ela está bloqueada por uma
desigualdade, não por calibração (§3.3.6 acima, e lição #64).

**As três rotas convergem no mesmo ponto: `cs_max`.** Enquanto o teto existir, mais biomassa
satura `cs` por área e achata `∇cs`; sem o teto, ou a faixa dinâmica explode (Y3, `a_mar`
116) ou o perfil inverte e o push vira para dentro. **Esse é o problema aberto do projeto, e
ele é anterior a C5.**

---

# REFORMULAÇÃO (2026-09-02) — o enquadramento estava errado

> **Correção da usuária, e ela invalida a premissa das rotas B e C.** Escrevi que "a
> colônia não expande, é construída — 98% do material em `r>2` foi inserido lá". **Isso
> não é defeito.** Colônia bacteriana expande principalmente por **divisão celular**: a
> bactéria da frente se divide e a filha ocupa o espaço adjacente. **Criar partícula na
> frente É o mecanismo correto.**
>
> Eu tratei criação como sintoma e advecção como o "certo", e por isso gastei duas rotas
> inteiras (B — o meio; C — a deposição) atacando **transporte**, que nunca foi o problema.
> A invariância da advecção (31 → 31 → 34) permanece como medição válida, mas **não é o
> diagnóstico** — advecção não deveria ser o mecanismo.
>
> A usuária também estabeleceu que **engrossar os braços não é problema**, contanto que não
> fiquem grossos demais. Isso relaxa a lei (L) (lição #67-L) como objeção a mais biomassa.

## O defeito reformulado: escala de comprimento do movimento

Não é *de onde vem* o material — é **quão longe cada partícula viaja antes de parar**:

| | medido em `runs/E5_hillK025`, t=50 |
|---|---:|
| deslocamento da viva mediana | 7.1 dx |
| p90 | 35.6 dx |
| **máximo** | **77.5 dx** |
| diâmetro de uma partícula | 1 dx |

Uma bactéria que se divide põe a filha a **~1 diâmetro**. No modelo, **13% das vivas viajam
50–77 diâmetros** — elas não se dividem, **cavalgam**, e o que fica para trás é vazio que só
o wake preenche, em grumos (612 componentes de mediana 2).

**A métrica do defeito é a razão frente/corpo:**

```
v_frente = 69 dx / 50 s = 0.074       (a colônia vai de R99 0.43 a 4.15)
v_corpo  = 0.0025                     (|v| mediano das vivas moveis)
razao    = 30x                        (alvo: <= 3)
```

Num swarm real essa razão é ~1: todo mundo avança junto, devagar, dividindo. Aqui uma
vanguarda de **35 partículas** (13%) corre de `r`=0.33 a 2.91 e o resto fica parado.

**Isso explica por que B e C falharam:** as duas mexiam em *como o material chega na frente*,
quando o problema é que **o corpo não acompanha**.

## Causa direta: `cs` é plano dentro da colônia

| zona | `cs_∞` (E5) |
|---|---:|
| núcleo | 0.486 |
| junção | 0.490 |
| mid-arm | 0.490 |
| rim | 0.481 |
| ágar | 0.316 |

Como a força é `−β∇cs`, um campo plano dá **força zero no interior**. Só o rim, onde `cs`
cai para o ágar, é empurrado — e por isso só ele anda. O gate da Marangoni **não** é o
culpado: medido pelo `au_mar` real, **265 das 276 vivas já recebem força**. Falta gradiente,
não permissão.

*(Nota de medição: `grad_rho_b_mag` e `grad_cs_x` são acumuladores zerados no `initialize`
e aparecem como 0 em todo dump — não servem para diagnóstico pós-hoc. Usar `au_mar`.)*

## Propostas reformuladas

### Proposta 1 — perfil monótono de `cs` (principal)

Três mudanças **acopladas**, que só funcionam juntas:

| | de → para | por quê |
|---|---|---|
| teto `(1 − cs/cs_max)` | **remover** | é o sumidouro dominante (`P/cs_max ≈ 14` contra `λ_eff` = 0.30) e achata tudo |
| `c_n_factor` na produção de `cs` | **remover** | `c_n` é 0.1 no núcleo e 0.8 nos braços, então suprime justo onde a produção deveria ser máxima — cria o bump da junção (17.4 → 24.6) |
| `β·σ` | **÷ 4.2** (β 5 → 1.2) | recalibra `a_mar` do corpo para o orçamento §8 |

**Perfil resultante (§3.3.6):** núcleo **34.8** → junção 32.8 → mid-arm 29.6 → rim 14.4 →
ágar 0.94. **Único monótono decrescente** das quatro configurações avaliadas: push outward
em **toda** a colônia (§3.3.4, linha correta).

**Predição:** `|∇cs|` no corpo ≈ 5.1 (hoje ~0); `a_mar` do corpo ≈ 6 e da fronteira ≈ 15;
razão frente/corpo de 30× para **~2.5×**.

**Riscos:** o Y3 removeu o teto e teve `a_mar` = 116 — mas sem recalibrar `β` nem remover
`c_n`. Remover `c_n` da produção **aproxima** o modelo do Xavier/CCR (a tensão com o
acoplamento direto está documentada em §3.2). Guardrails de `mean_cs`/`contrast_cs` precisam
ser reescalados: o campo fica ~70× maior.

### Proposta 2 — crescimento na frente (desbloqueada pela reformulação)

Com "engrossar não é problema", a objeção da lei (L) cai. O que restava contra `r_growth`
alto era o achatamento de `∇cs` (lições #43/#64) — que é exatamente o que a Proposta 1
resolve. **As duas são a mesma solução:** um perfil de `cs` que sobrevive a mais biomassa, e
mais biomassa para **dividir** em vez de cavalgar.

### Proposta 3 — distribuir a força, em vez de concentrá-la na vanguarda

O gate flagelar (`rho_b ∈ [0.2, 0.6]`) é o que **seleciona** a vanguarda: só **72
partículas** têm `au_flag > 0`. A `FlagellarForce` tem magnitude constante (`f0·gate·n̂`),
então é imune à saturação de `cs` — precisa só da direção. Alargar para `[0.1, 0.8]`
distribui a força por ~230 das 276 vivas.

**Nunca testado sob o E5** — o K.2 só testou **estreitar**. **Risco:** lição #15 (brittle
neck) — mais tração exige reavaliar `tension_ratio` no mesmo passo.

## Ordem e métrica

**Proposta 3 primeiro** (alavanca única, barata, mede a razão frente/corpo diretamente),
depois **1+2 juntas**, que são acopladas.

**Métrica nova:** razão `v_frente / v_corpo`, medida como **p90/p50 do deslocamento das
vivas** em t=50. Hoje **35.6/7.1 = 5.0** pela distribuição, ou **30×** pela velocidade
característica. Alvo ≤ 3. Ela captura o defeito melhor que `C5a`/`C5b`, que medem a
consequência.

---

## P3 — alargar o gate flagelar: REPROVADA, e o motivo argumenta a favor da Proposta 1

> **Executado 2026-09-02**, `runs/P3_flaggate`. Alavanca única sobre o E5: gate flagelar
> `[0.2, 0.6] → [0.1, 0.8]`, parametrizado em `FLAG_GATE_LO`/`FLAG_GATE_HI`.

| | E5 | **P3** |
|---|---:|---:|
| partículas com `au_flag > 0` | 72 | **169** (2.3×) |
| deslocamento **p50** (corpo) | 7.11 dx | **6.92** |
| deslocamento **p90** (frente) | 35.6 dx | **47.3** |
| **RAZÃO p90/p50** | 5.00 | **6.84** |
| % que viaja > 30 dx | 12.7% | **19.2%** |
| `a_mar` | 12.31 | 8.95 |
| `C5a` / `C5b` | 0.24 / 11.3% | 0.21 / 10.2% |
| AR / dedos / `R99` | 7.61 / 24 / 4.15 | **7.60 / 26 / 4.14** |
| `a_press` med / picos>4 | 2.58 / 7% | 2.64 / 7% |

**A alavanca engajou e piorou o alvo.** O corpo não se moveu (p50 6.92 contra 7.11) e as 97
partículas recrutadas viraram **mais cavaleiros**: p90 de 35.6 para 47.3, fração que viaja
> 30 dx de 12.7% para 19.2%. A razão p90/p50 foi de 5.00 para **6.84** (alvo ≤ 3).

O risco pré-registrado (brittle neck, lição #15) **não se materializou** — AR, dedos,
`R99` e `a_press` todos preservados.

### A lição: força constante não move corpo, cria cavaleiros

A `FlagellarForce` tem magnitude **constante** (`f0·gate·n̂`): toda partícula no gate recebe
3.0 e adquire velocidade terminal `f0/gamma` **independentemente das vizinhas**. Adicionar
partículas ao gate não faz o corpo andar junto — faz **mais partículas cavalgarem sozinhas**.

**Para mover um corpo coerentemente a força precisa ser GRADIENTE DE CAMPO, não constante
por partícula.** É o gradiente que faz a força variar suavemente no espaço e o corpo se
deslocar como bloco; uma constante numa banda produz N cavaleiros independentes.

**Isso argumenta A FAVOR da Proposta 1.** A P3 testou "mais força, mesma natureza" e mostrou
que a natureza é que está errada. O que falta é um `∇cs` que exista no interior — que é
exatamente o que o perfil monótono entrega.

Gate revertido para `[0.2, 0.6]`; código parametrizado fica (§10).

---

## P1 — perfil monótono de `cs`: REPROVADA, e o que ela revelou

> **Executado 2026-09-02**, `runs/P1_monotono`. `CS_CEILING`=0 (sem teto), `CS_PROD_CN`=0
> (produção sem `c_n`), `sigma` 11.1 → 2.64 (β·σ ÷ 4.2).

| | E5 | **P1** |
|---|---:|---:|
| razão p90/p50 | 5.00 | **4.67** |
| desloc p50 (corpo) | 7.11 | **7.81** |
| `C5a` / `C5b` | 0.24 / 11.3% | 0.29 / 12.2% |
| AR / dedos | 7.61 / 24 | **7.82 / 25** |
| `a_press` med / picos>4 | 2.58 / 7% | 2.59 / **0%** |
| **`a_mar`** | 12.31 | **75.90** |
| **`R99`** | 4.15 | **3.30** |
| VIVAS / biomassa | 276 / 0.472 | 178 / 0.289 |

**Direção certa, magnitude incontrolável.** O corpo acelerou (p50 7.81), a razão caiu, AR e
dedos preservados, picos de pressão zeraram — mas `a_mar` foi a **75.9** e a colônia
encolheu 20%.

**O perfil saiu diferente do previsto:**

| banda | previsto | medido |
|---|---:|---:|
| núcleo | 8.3 | **8.25** ✓ |
| junção | 7.8 | **6.84** ✓ |
| **mid-arm** | **7.0** | **0.284** ✗ (25×) |
| rim | 3.4 | 0.333 |

Errei por usar `cs_∞ = P/λ_eff`, que **ignora a perda difusiva** — dominante onde a biomassa
é rala e vizinha do ágar (`D_ext`=0.08). Mesma armadilha da lição #55.

**E isso revela o que o teto realmente fazia:** `P/cs_max ≈ 14` é um sumidouro local tão
grande que domina a difusão, tornando o campo localmente determinado e portanto **uniforme**.
Era ele que sustentava `cs` nos braços. Sem ele, a estrutura difusiva real aparece — e é
dominada pelo núcleo, com o gradiente inteiro concentrado na transição núcleo/braço
(razão de gradientes junção:frente ≈ 25:1). Como o sistema é **linear em σ**, mudar σ não
altera essa razão.

**Efeito colateral não previsto:** `CS_GATE_REF = 8.3` pôs o gate da colonização em `cs > 5.0`
e os braços têm `cs` = 0.28 — **desliguei o recrutamento nos braços sem querer**. É o que
explica VIVAS 178 contra 276, e o limbo com `ρ_b` plano em zero na trajetória.

---

# SÍNTESE — o que é necessário para a continuidade

## Os três requisitos, com número medido

**1. A população precisa crescer com a ÁREA.** Área cresce **93×** (R 0.43 → 4.15), população
viva cresce **1.8×** (152 → 276). Portadoras por dx de frente caem de 3.02 para 0.57.

**2. O corpo não pode estar congelado.** **37% das vivas têm `|v|` exatamente zero** — 45 pelo
pin de `ρ_b ≥ 0.8`, 56 pelo pin químico `c_n < 0.6`, que dispara dentro dos braços
(`c_n` p25 = 0.486).

**3. A força tem que ser gradiente de campo no interior.** A P3 provou que magnitude constante
cria cavaleiros independentes; o E5 tem `cs` plano dentro e só a interface é empurrada.

## O requisito REAL, que corrige o item 1

A pergunta "os braços serão preenchidos desde o início da formação?" tem resposta medida:

| | |
|---|---:|
| partículas recrutadas do zero | 783 |
| **que chegaram ao quórum** | **33 (4.2%)** |
| tempo recrutamento → quórum | p50 = **29.0 s** |
| tempo da frente passar por 1 dx | **0.73 s** |
| **razão** | **40×** |

**A frente passa em 0.73 s; a recruta leva 29 s para virar colônia** — e 96% nem chegam,
param em `ρ_b` ≈ 0.046. O braço seria preenchido ~29 s **depois** de formado, e a trajetória
continuaria mostrando desconexão.

**Logo `r_growth` não é a alavanca.** 5× reduz 29 s para ~6 s, contra os 0.73 s necessários;
fechar de vez exigiria ~40× (`r_growth` = 0.8), taxa de duplicação de 0.9 s — absurdo
biológico e explosão numérica.

> **O requisito correto: a portadora nova precisa NASCER no quórum, não crescer até ele.**
> Biologicamente é o que acontece — a filha nasce com a densidade da mãe, não como traço que
> amadurece.

## A célula não testada

A `BiomassColonization` com **alvo-doador** (lição #59, K3) já faz isso: relaxa para
`Σ V ρ² W / Σ V ρ W`, a densidade do doador. E o K3 mediu que **funciona perto do núcleo**
(p90 = 0.217 em r=0.6, acima do quórum) e **falha nos braços** (0.097 em r=1.0).

O motivo está na própria lição: o **filler foi excluído como doador**, porque 93% dos doadores
eram valores históricos congelados. Só que **nos braços o filler É o material** (`ρ_b` ~0.3), e
sem ele sobram as poucas vivas esparsas.

**Teste:** recrutamento com alvo-doador **incluindo o filler dos braços como doador**, para a
recruta nascer a `ρ_b` ~0.3 — colônia no instante em que nasce, acompanhando a frente em vez
de vir 29 s depois.

**Risco conhecido (lição #59):** recrutas sub-quórum afogam `cs`. Menos grave agora, porque
sem o teto o `cs` não satura por área — mas é o que precisa ser medido.

---

## P2 e P4 — recrutamento e ancoragem

> **2026-09-02**, `runs/P2_fillerdonor` e `runs/P4_attach`. Medidos sob a definicao fechada
> nesta data: colonia = `rho_b >= 0.1` ou `is_filler`; `0 < rho_b < 0.1` e **agar engolido**.

**P2 (`COL_FILLER_DONOR=1`)** — o filler passa a contar como doador na `BiomassColonization`,
que relaxa para a densidade do doador. Nos braços o filler **e** o material (`rho_b` ~0.3);
sem ele os doadores eram poucas vivas esparsas e a recruta nascia em ~0.1.

**P4 (`WAKE_ATTACH=1.05`)** — o wake so deposita se o material mais proximo estiver a
**menos** de 1.05 dx (antes so havia limite inferior, `WAKE_PROX=0.7`). Medido no E5: **69%
dos depositos nasciam a mais de 1.05 dx de qualquer material**, mediana 1.44, cauda 3.7.

| | E5 | **P2** | P4 |
|---|---:|---:|---:|
| **colonia** (particulas) | 3350 | **4744** (+42%) | 3898 |
| traco (`rho_b`<1e-6) | 1436 | **201** (−86%) | 152 |
| colonia / disco | 15.6% | **23.0%** | 20.7% |
| maior comp. @1.15dx | 30% | **58%** | 68% |
| **alcance** @1.15dx | 38% do raio | **83%** | 54% |
| `C5a` medio (t 35–50) | 0.24 ± 0.04 | **0.45 ± 0.05** | 0.35 ± 0.16 |
| AR / dedos | 7.61 / 24 | 8.16 / 25 | **8.79** / 25 |
| `a_press` med / picos | 2.58 / 7% | 2.59 / **0%** | 2.60 / 0% |
| `frac(cs>0.45)` | 52.4% ❌ | **42.3%** ✅ | 45.6% ✅ |
| wake nasce conectado | 31% | 43% | **61%** |

**P2 e o ganho estabelecido.** `C5a` 0.24 → 0.45 com separacao de 4× o ruido; converte agar
engolido em colonia (+42%, traco −86%); AR melhora; `frac(cs>0.45)` passa dos 50% pela
primeira vez; picos de pressao zeram.

**P4 engaja o mecanismo** (depositos nascendo conectados 43% → 61%) e da a melhor fracao na
escala estrita (41% vs 26%), **mas encurta o alcance** — 54% do raio contra 83% do P2.
A ancoragem rejeita justamente os depositos distantes, que eram os que estendiam a
componente pelos braços. **Erro de leitura a nao repetir:** julguei o P4 melhor pela
FRACAO, que favorece quem tem menos material (3938 vs 4792); o `C5a` (alcance) dizia o
contrario e estava certo.

**O que nenhum dos dois resolve:** os braços seguem desconexos em toda escala ate 1.5 dx.
