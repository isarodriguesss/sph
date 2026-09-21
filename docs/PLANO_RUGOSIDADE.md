# Plano — Pass L: superfície rugosa como rede de pilares no plano

> **Estabelecido 2026-09-14.** Pass L desbloqueado por decisão da usuária. Pelo registro, os
> critérios do §2.2 NÃO estão atingidos (sem tip-splitting; o P2 vira disco depois de t≈55,
> lição #90). Portanto a rugosidade entra como **extensão física a ser caracterizada**, não
> como validação da morfologia de referência, e toda comparação é **rugoso × liso com a mesma
> semente**, nunca rugoso × `reference.jpg`.

**Decisões tomadas (2026-09-14):** obstáculos no plano; rede periódica de pilares; bloqueio da
difusão decidido pela literatura (seção 1).

> **ATUALIZAÇÃO 2026-09-17 — baseline e geometria reancorados.** O baseline passou a ser o
> **P2R23** (P2R21 + filler conduzindo `c_s`), e ele tem braços **4 a 6× mais largos** que o P2 para
> o qual este plano foi dimensionado. Medido em t=50 (mesma régua nos dois):
>
> | em r=2.5 | P2 (base original do plano) | **P2R23 (baseline atual)** |
> |---|---:|---:|
> | braços | 19 | 16 |
> | largura do braço | 1.6 dx | **10.5 dx** |
> | passo entre braços | 15.4 dx | 18.2 dx |
> | fresta livre entre braços | 13.8 dx | **7.8 dx** |
>
> Com isso, o pilar de 5 dx da tabela da seção 2 deixaria de ser obstáculo (viraria metade da largura
> do braço) e Λ = 18 dx entraria em **ressonância** com o passo entre braços (18.2 dx) — exatamente o
> acoplamento que a escolha de rede triangular pretendia evitar. Geometria reancorada, preservando as
> três razões originais (pilar ≈ largura do braço; fresta ≫ largura; Λ fora de ressonância):
>
> | parâmetro | valor original | **novo** | derivação |
> |---|---:|---:|---|
> | diâmetro A | 5 dx | **10 dx** (R_p = 5 dx ≈ 0.27) | ≈ largura do braço (9-10.5 dx): o braço não engole o pilar nem é barrado por ele |
> | espaçamento Λ | 18 dx | **28 dx ≈ 1.51** | livre caminho médio `1/(n·(A+w))` = 35 dx com `n = 1/(Λ²·0.866)`, para ~1.5 contatos por braço no percurso de 52 dx até t=50 |
> | fresta Λ − A | 13 dx | **18 dx** | 1.8× a largura do braço — passa, mas com atrito |
> | fração de área | 7% | **~11.6%** | `π R_p²/(Λ²·0.866)` |
> | nº de pilares em [-7,7]² | ~240 | **~95** (~78 partículas cada, ~7400 no total) | — |
>
> A exclusão em r < 1.2 e a rede triangular seguem como estavam. O controle liso da série é o
> **P2R23** (`runs/swarm/P2R23_conduz`), mesma semente (20260806) e 10 threads. As rodadas desta
> série vão para **`runs/rugosidade/`**.
>
> **V0 executado em 2026-09-17: PASSOU.** Com `use_pilares = False`, as 62 arrays do HDF5 e as 45
> colunas do log saem bit-a-bit identicas ao P2R23 nas iterações 200 e 400 (10 threads), inclusive o
> `t` do dt adaptativo. Implementado: constantes `PILAR_*` e `SwarmApp._faz_pilares` em `main.py`
> (rede triangular; os próprios pontos da grade viram o array `pilar`; grava `pilares.json`), e os
> pilares como fonte de `SummationDensity`, `KernelSum` e `ViscousForce` em `src/scheme.py`.
> **Falta a força de contato** (seção 4), que exige a calibração offline de `K` antes de V1.

---

## 0. Pré-requisitos antes de escrever código

1. Registrar o desbloqueio no CLAUDE.md: §1 Objetivo 1, §2.2 (gatilho do Pass L), §10
   (proibições) e §12 (Pass L), com a ressalva acima.
2. Hook `guard_pass_l.py`: usar a escotilha `# pass-l-aprovado: <motivo>` que ele já prevê
   (uma vez por trecho novo com as palavras-chave), ou desligá-lo em
   `.claude/settings.json` — decisão da usuária.
3. Congelar a semente (`SEED = 20260806`, já fixa) e o número de threads (10) para toda a série.

---

## 1. Modelo físico

**O que um pilar representa.** Uma elevação da superfície mais alta que o filme da colônia
(alguns µm). Na vista de topo 2D, isso é uma região que a colônia não pode ocupar e que ela
precisa contornar. Rugosidade de amplitude MENOR que a espessura do filme não é representável
neste modelo 2D sem campo de altura — fica fora do escopo.

**Por que no plano e não nas paredes do domínio.** A colônia só alcança a parede em t≈75,
fora da janela morfológica do P2 (t ≲ 55). Paredes onduladas não teriam efeito mensurável.

**Difusão química: fluxo zero para `cs` e `c_n` nos pilares.**

| campo | física | no modelo | referência |
|---|---|---|---|
| `cs` (ramnolipídeo) | camada superficial própria, sob os tendrils, que avança junto com eles (r = 0.99) | bloqueado: um pilar mais alto que o filme interrompe a camada | Kasallis, Bru et al. 2023 |
| `cs` junto a quinas | biossurfactante gera escoamento ao longo de quinas CÔNCAVAS < ~60° (0.6-1.3 mm/h) | não se aplica: pilar cilíndrico sobre substrato plano faz 90° | Li, Sanfilippo, Kearns & Yang 2022 |
| `c_n` (nutriente) | fica no volume do gel e passa por BAIXO de relevo apoiado sobre o ágar; matriz sólida que atravessa a espessura é impermeável | bloqueado no plano, mas o suprimento vertical do ágar já é a fonte local `k_src(1-c_n)`, que o pilar não remove | Bhattacharjee et al. 2020 [T2 para `k_src`] |

Implementação do fluxo zero: deixar os pilares FORA das somas de Brookshaw da
`SurfactantEquation` e da `OxigenConsumption` — omitir as partículas de parede do laplaciano é o
contorno adiabático de Cleary & Monaghan (1999). Custo zero de código.

**Estimativa ex-ante do efeito químico: PEQUENO.** O raio do pilar (2.5 dx = 0.13) é muito menor
que o alcance difusivo do `cs` no ágar (`L = sqrt(D_ext/λ_ágar) = sqrt(0.08/0.075) ≈ 1.03`) e do
nutriente (`sqrt(D_n/k_src) = 0.41`). O campo contorna o pilar; a sombra química a jusante e o
acúmulo a montante devem ser de ordem `R_p/L ~ 0.1-0.3`. **O efeito dominante esperado é
MECÂNICO.** Verificação pós-rodada: `cs` e `c_n` num anel de 2 dx em volta dos pilares contra o
mesmo raio longe deles.

---

## 2. Geometria

| parâmetro | valor proposto | justificativa |
|---|---|---|
| rede | **triangular** (6 vizinhos) | a rede quadrada tem simetria de ordem 4, que casa com a semente `cos(8θ)` e travaria dedos nos eixos da rede por construção; a triangular (ordem 6) não é submúltiplo de 8 |
| espaçamento Λ | **18 dx ≈ 0.97** | espaçamento entre braços em r≈2.5 com ~16 braços: `2π·2.5/16 ≈ 1.0` |
| diâmetro A | **5 dx ≈ 0.27** (R_p = 2.5 dx) | ~2x a largura atual do braço (2.5-2.8 dx): o pilar é obstáculo real, mas a fresta Λ−A = 13 dx deixa o braço passar |
| exclusão | nenhum pilar com centro em r < **1.2** | não perturbar o inóculo nem a junção núcleo-braço (onde o P2 faz a ligação) |
| fase | rede gerada a partir da origem com rotação fixa, gravada no run | reprodutível; uma rodada com a rede girada 30° separa efeito de rede de efeito de orientação |

Fração de área ocupada: `π R_p² / (Λ² · 0.866) ≈ 7%`. Com o domínio [-7,7]² saem ~240 pilares
(~20 partículas cada, ~4800 no total).

**Construção sem sobreposição.** Os pilares são feitos com os PRÓPRIOS pontos da grade de fluido
que caem dentro de cada disco: esses pontos saem do array `fluid` e entram num array novo
`pilar`. Mesma rede, mesma massa `dx²`, então a densidade vista pelos vizinhos não muda.
Posições e raio dos pilares vão para `pilares.json` no diretório do run (pós-processamento).

**Chave de liga/desliga.** `use_pilares = False` não cria o array → código idêntico ao baseline,
verificável bit-a-bit. O array `solid` das paredes do domínio fica como está (inerte).

---

## 3. Numérica SPH — quem enxerga os pilares

| equação | inclui `pilar` como fonte? | motivo |
|---|---|---|
| `SummationDensity` | **sim** | sem isso o fluido junto ao pilar lê rarefação falsa (partição da unidade truncada, Violeau §3.6); paredes na soma de densidade como em Adami, Hu & Adams 2012 |
| `KernelSum` (`sigma_a`) | **sim** | idem; e sem isso a inserção C3.4 dispara no `sigma_a<0.85` junto ao pilar e põe filler DENTRO dele |
| `KernelGradientCorrection` | não | os gradientes somam só fluido; a KGC corrige exatamente o suporte truncado (Bonet & Lok 1999, teste-ouro no §12) |
| `MarangoniForce`, `FlagellarForce`, `BiomassGradient` | não | `cs` e `rho_b` não existem no pilar; KGC cobre a Marangoni |
| `SurfactantEquation`, `OxigenConsumption` | não | fluxo zero (seção 1) |
| `BiomassColonization`, `BiomassGrowth` | não | pilar não é doador nem recruta |
| `MomentumEquation` (pressão + visc. artificial) | não | a pressão da EOS é ~219x menor que a Marangoni (lição #85): parede por pressão não segura a colônia |
| `ViscousForce` | **sim** (`u=v=0` no pilar) | condição de não-deslizamento por atrito viscoso com a parede |
| **`ForcaContornoPilar`** (nova) | **sim** | contato mecânico — seção 4 |
| `ParticleShift` | **sim** | com fonte só-fluido, `-∇C` aponta PARA DENTRO do pilar e o shifting empurraria partículas para lá |

O integrador só avança o `fluid`: os pilares ficam estáticos por construção.

---

## 4. Força de contato

**Forma: força de contorno por partículas de Monaghan & Kajtar (2009)**, que usa o próprio kernel
e não tem a singularidade do Lennard-Jones de Monaghan (1994). Para cada par fluido *i* –
pilar *k*:

```
f_ik = (K / β) · (x_ik / r²) · W(r, h) · 2 m_k / (m_i + m_k)
```

com `β` = razão entre o espaçamento das partículas do pilar e o do fluido (= 1 aqui). A LJ
12-4 de Monaghan 1994 (Liu & Liu 2003 [T6], cap. 4, partículas virtuais) fica como alternativa.

**Por que a calibração de `K` é o ponto crítico.** O corpo da colônia opera comprimido
(`rho/rho0` p99 de 5-7 no P2, lição #88; espaçamento local ~0.5-0.7 dx). O ponto de vazamento é
o meio do vão entre duas partículas da superfície do pilar, a **0.5 dx** de cada uma. `K` precisa:

1. **Segurar:** força em r = 0.5 dx maior que a aceleração motriz máxima (`a_mar + a_flag`, ~15
   no máximo do log do P2).
2. **Não colapsar o dt:** o PySPH limita `dt_force = 0.25 · sqrt(h/|a|)`. Com |a| ≤ 60 no contato
   mais profundo, `dt_force ≈ 0.25 · sqrt(0.097/60) ≈ 0.010`, comparável ao dt atual (0.007-0.016).

Calibração feita OFFLINE antes de rodar: calcular a curva `|f(r)|` para o kernel cúbico com
h = 1.8 dx e escolher o `K` que satisfaz (1) e (2). Se não houver `K` que satisfaça os dois,
é sinal de que o contato precisa de r mínimo maior (pilar com 2 camadas de superfície
desencontradas, fechando o vão de 0.5 dx), não de `K` maior.

**Instrumentação da força:** arrays `ax_rep`/`ay_rep`, e o resíduo `a_press_max`/`a_press_med`
passa a descontá-la (senão a coluna de pressão vira a coluna de contato).

---

## 5. Armadilhas — tudo que hoje pode pôr material DENTRO de um pilar

| mecanismo | como falha | correção | estado |
|---|---|---|---|
| inserção C3.4 (`_insere_vacuo`) | candidatos testados contra a árvore de fluido; `sigma_a` baixo junto ao pilar dispara | `_livre_pilar` com o MESMO `prox` já aplicado ao fluido | ✅ 2026-09-18 |
| wake (`_deposita_rastro`) | o anel de até 7 a 0.75 dx e (com `WAKE_SEG`) os pontos do segmento podem cair no pilar | o mesmo teste dentro de `livre()` | ✅ 2026-09-18 |
| shifting | `-∇C` só-fluido aponta para dentro do pilar | `ParticleShift(sources=src_par)` | ✅ 2026-09-18 |
| `void_fraction` (C1 do §2.5) | área sólida conta como vazio da colônia | máscara no denominador | ✅ 2026-09-18 |
| figura e réguas de forma | pilar não tem ágar → lido como COLÔNIA | `campo()` perfura o sólido; `carrega()` acha `pilares.json` sozinha | ✅ 2026-09-18 |
| largura EDT (`ciclo_largura`) | `binary_fill_holes` preenche o pilar cercado → largura inflada | perfura DEPOIS do preenchimento | ✅ 2026-09-18 |
| `R99` | — | **nenhuma**: o pilar sai de `fluid` no `extract_particles` (array `pilar` próprio no HDF5) e `_registra` mede sobre `fluid` | ✅ já estava correto |
| `n_pen`, `n_pen_sup`, `n_contato`, `a_rep_*` no log | ainda não existem | seção 6 | ✅ 2026-09-18 |
| resíduo `a_press_*` contando a força de contato | a coluna de pressão viraria a de contato | subtrai `ax_rep`/`ay_rep` | ✅ 2026-09-21 |

**Magnitude do que foi removido, medida no V1 (t=20, 92 pilares de 5 dx):**

| | sem máscara | com máscara | baseline liso |
|---|---:|---:|---:|
| `void_15` | **8.576%** | **0.000%** | 0.02% |
| `void_10` | 11.074% | 0.338% | 0.25% |

Sem isso **toda dose reprovaria em C1** por área que é sólida. E a guarda de deposição não era
hipotética: o V1 rodou sem ela e já tinha **5 fillers dentro de pilares em t=20** (0 em t=12.4),
com 45 a menos de `prox` — crescendo, e o V1 só chegou a 40% do tempo previsto.

**Verificação (2026-09-18):** com `use_pilares=False` o código é **bit-a-bit idêntico** ao anterior
— `runs/swarm/_verif_mask_novo` × `_verif_mask_ctrl`, 10 threads, iterações 0/200/400/600:
`log.csv` idêntico caractere a caractere nas 47 colunas e `max|dif| = 0.000e+00` nas 79 arrays de
cada frame, com o `t` do dt adaptativo batendo até a última casa. Wake disparou em 100-600 e
insert em 200/400/600, então as guardas estiveram no caminho de execução.

**Segunda verificação (2026-09-21)**, cobrindo a instrumentação e o desconto do resíduo,
que foram acrescentados depois da primeira: `runs/swarm/_verif_instr_novo` × `_verif_instr_ctrl`
(HEAD), 10 threads, iterações 0/200/400 — **47 colunas comuns do `log.csv` idênticas**,
`max|dif| = 0.000e+00` nas 79 arrays comuns, e as **5 colunas novas todas em zero** com
`use_pilares=False`. `LOG_HEADER` e `writerow` alinhados em 52 (o modo de falha do §9,
*CSV column drift* do Pass F, está descartado).

Toda propriedade persistente nova tem que ser setada nos dicts de inserção (lição #67-I) — aqui
não há propriedade nova no `fluid`, só o array `pilar`.

---

## 6. Instrumentação e pós-processamento

- **Log:** `n_pen` (fluido com centro a < R_p − 0.5 dx de um centro de pilar — deve ser 0),
  `n_contato` (fluido a < 1 dx de uma partícula de pilar), `a_rep_max`, `a_rep_med` (nos em
  contato). Colunas novas no FIM do `LOG_HEADER`, como as de pressão.
- **`pilares.json`** no run: centros, raio, Λ, rotação, exclusão.
- **`plots/fig_tese.py`:** pilares desenhados em cinza nos 4 painéis e excluídos do campo da
  colônia e do fechamento de buracos de 3.5 dx.
- **Métricas novas** (script `tools/diag_pilares.py`):
  - eventos de contato por braço (braço que toca um pilar entre dois quadros) e o desfecho:
    **desvio** (o braço contorna), **divisão** (dois braços saem de um — tip-splitting
    induzido por obstáculo), **parada** (a ponta estaciona);
  - espectro azimutal: potência em m = 6 (assinatura da rede) contra o liso;
  - `cs` e `c_n` a montante/jusante dos pilares (verificação da seção 1).

---

## 7. Sequência de rodadas

**O estudo é uma COMPARAÇÃO: a rugosidade atrapalha o swarm?** (decisão da usuária, 2026-09-18).
Isso reordena o desenho: o padrão é o pilar, a alavanca é a **dose**, e a orientação (sulcos) sai
do caminho crítico — ela responde "como desvia", não "se atrapalha".

**Por que dose-resposta e não um par liso × rugoso.** Com `Λ = 28 dx` a fresta vale 18 dx = **3.3
larguras de braço** (braço = 5.4 dx no P2R23): o braço passa sem encostar. Um "não atrapalha"
nessa geometria mede que a colônia não tocou os obstáculos — geometria, não biofísica, e é a
lição #102 (`P2R18`) com outro rótulo. Um resultado negativo só vale se o experimento tinha poder
de detectar o efeito, e é a monotonicidade ao longo das doses que dá isso.

| rodada | Λ | fresta/braço | φ | pilares | **contatos/braço** |
|---|---:|---:|---:|---:|---:|
| **liso** | — | — | 0 | 0 | **0** — controle = **P2R23** |
| **V0** ✅ | — | — | 0 | 0 | bit-a-bit com `use_pilares=False` |
| **V1** ✅ | 28 dx | 3.3 | 11.6% | 92 | `n_pen`=0, dt inalterado, t≈20 |
| **R1** | 28 dx | 3.3 | 11.6% | 92 | **1.50** |
| **R2** | 20 dx | 1.9 | 22.7% | 180 | **2.38** |
| **R3** | 16 dx | 1.1 | 35.4% | 256 | **4.02** |
| **R4** | — | — | (casada) | 0 | **controle de atrito** (ver abaixo) |

Contatos por braço medidos por Monte Carlo sobre a rede real (pré-passo de 2026-09-18,
`docs/CRITERIOS_RUGOSIDADE.md` seção 8). **Nenhuma dose é transparente** — 100% dos braços
interceptam ao menos um pilar nas três. A fresta estática enganava: ela diz que o braço *pode*
passar entre dois pilares, não que passe ao percorrer 53 dx radiais numa rede 2D.

`φ = (π/2√3)(A/Λ)²`, com `A = 10 dx` FIXO em todas — 2× a largura do braço, que é a transposição
da regra da literatura (feição = 1-4× o tamanho do objeto que se move; Jayathilake usa 2× e 4×).
Uma alavanca por vez: só `Λ` varia, e a dose é reportada como `φ`.

**R4 — controle de atrito pareado, e ele não é opcional.** Os pilares entram como fonte da
`ViscousForce`, então parte da desaceleração de R2/R3 é **atrito adicional, não bloqueio
geométrico**. R4 roda um campo contínuo `γ_eff(x) = γ·[1 + κ·φ_r(x)]` com `κ` calibrado para
igualar o sumidouro de momento dos pilares em `φ` casado, **sem nenhum sólido**. A diferença
R3 − R4 isola a contribuição geométrica. Sem ele, "a colônia ficou mais lenta" não vira "a
colônia foi bloqueada". Detalhe da arquitetura e por que a Opção 1 ficou como controle e não
como mecanismo: análise preditiva de 2026-09-18 (Opção 1 × 2 × 3).

**Sulcos (radiais e concêntricos)** ficam para depois de fechada a curva de dose: mesma
infraestrutura, trocando só o recorte, com a orientação como alavanca única.

Controle liso = **P2R23**, mesma semente e **mesmo número de threads** (lição #88: thread
diferente já é outra realização). Como a morfologia diverge caoticamente assim que algo muda
(§2.2), toda comparação usa janela **t ∈ [35, 50]** com média ± desvio, contra o ruído entre as
duas realizações do P2 como piso.

**Confundidor de normalização:** a rodada rugosa nasce com menos fluido (−11.6% em R1 a −35.4%
em R3), porque `_faz_pilares` extrai as partículas. Nenhuma comparação pode usar contagem ou
massa absoluta — só grandezas **por área** ou já mascaradas.

---

## 8. Predições e critérios

> **SUPERSEDIDA em 2026-09-18** pela série de dose-resposta. Os critérios pré-registrados da
> série R1-R4 estão em **[docs/CRITERIOS_RUGOSIDADE.md](CRITERIOS_RUGOSIDADE.md)** — em `docs/`,
> e não em `runs/<R>/CRITERIOS.md`, porque `tools/archive_run.sh` faz `rm -rf` no destino e já
> apagou o CRITERIOS.md pré-registrado do P2R22. Copiar para o run **depois** de arquivar.
> O que está abaixo é o registro do desenho de rodada única, mantido como histórico.

### 8-antigo. Predições e critérios para R1 (rodada única, superseded)

**Predições:**
- **Contatos:** livre caminho médio de um braço até um pilar `≈ 1/(n·(2R_p + w)) ≈ 35 dx ≈ 1.9`,
  com `n = 1/(Λ²·0.866)` e w ≈ 3 dx. Um braço anda de r=1.2 a r≈4 até t=50 (~52 dx), então
  **~1.5 contatos por braço, ~20 eventos na colônia** — estatística suficiente para classificar
  desfechos.
- **Efeito químico pequeno** (seção 1): diferença de `cs` montante/jusante < 30%, de `c_n` < 10%.
- **R99(t=50):** −5 a −15% contra o liso (7% da área bloqueada + desvios).
- **Numérica:** `n_pen = 0`; dt dentro de 0.5x do liso.

**Critérios:**
1. Contato impenetrável: `n_pen = 0` em todos os quadros.
2. Numérica: sem colapso de dt; massa sem runaway; `a_press_med` (sem a força de contato) ≤ 2x o liso.
3. A colônia segue na classe (b) Fingering do Trinschek (§11) — a rugosidade modula, não destrói.
4. **O resultado científico é o desfecho dos contatos**, não um alvo numérico: qualquer das
   três respostas (desvio, divisão, parada) é informativa; registrar a distribuição.

---

## 9. Fora do escopo desta etapa

- Campo aleatório de pilares e espectro de rugosidade (R3+).
- Escoamento de canto por biossurfactante (Li et al. 2022) — só relevante com quinas < 60°.
- `c_n` atravessando o pilar (só se a verificação da seção 1 mostrar diferença > 10%).
- Paredes do domínio com relevo e interação colônia-parede.
- Rugosidade menor que a espessura do filme (exigiria campo de altura).

---

## Referências

- Adami, Hu & Adams 2012, *J. Comput. Phys.* 231:7057 — condição de parede generalizada (paredes na soma de densidade).
- Bhattacharjee, Amchin, Ott, Kratz & Datta 2020, *bioRxiv* 10.1101/2020.08.10.244731 — matriz sólida impermeável × hidrogel permeável a nutriente.
- Bonet & Lok 1999, *CMAME* 180:97 — correção do gradiente sob suporte truncado.
- Cleary & Monaghan 1999, *J. Comput. Phys.* 148:227 — condução em SPH e contornos adiabáticos.
- Kasallis, Bru, Chang, Zhuo & Siryaporn 2023, *Curr. Opin. Solid State Mater. Sci.* — camada de surfactante e tendrils avançando juntos.
- Li, Sanfilippo, Kearns & Yang 2022, *Microbiology Spectrum* — escoamento de canto induzido por biossurfactante.
- Liu & Liu 2003 [T6], cap. 4 — partículas virtuais de contorno (seção exata a confirmar).
- Monaghan 1994, *J. Comput. Phys.* 110:399 — força de contorno repulsiva (Lennard-Jones).
- Monaghan & Kajtar 2009, *Comput. Phys. Commun.* 180:1811 — força de contorno por partículas baseada no kernel.
- Violeau 2012 [T7], §3.6 — partição da unidade.
