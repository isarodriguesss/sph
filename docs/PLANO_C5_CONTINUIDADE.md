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
