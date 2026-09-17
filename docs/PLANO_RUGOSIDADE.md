# Plano — Pass L: superfície rugosa como rede de pilares no plano

> **Estabelecido 2026-09-14.** Pass L desbloqueado por decisão da usuária. Pelo registro, os
> critérios do §2.2 NÃO estão atingidos (sem tip-splitting; o P2 vira disco depois de t≈55,
> lição #90). Portanto a rugosidade entra como **extensão física a ser caracterizada**, não
> como validação da morfologia de referência, e toda comparação é **rugoso × liso com a mesma
> semente**, nunca rugoso × `reference.jpg`.

**Decisões tomadas (2026-09-14):** obstáculos no plano; rede periódica de pilares; bloqueio da
difusão decidido pela literatura (seção 1). **Baseline da rugosidade:** o vencedor entre o P2
e o P2S (`WAKE_SEG`, em rodada) — a decisão sai do veredito do P2S.

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

| mecanismo | como falha | correção |
|---|---|---|
| inserção C3.4 (`_insere_vacuo`) | candidatos testados contra a árvore de fluido; `sigma_a` baixo junto ao pilar dispara | pilar em `KernelSum` + teste geométrico: rejeitar ponto a < R_p + 0.5 dx de um centro |
| wake (`_deposita_rastro`) | o anel de até 7 a 0.75 dx e (com `WAKE_SEG`) os pontos do segmento podem cair no pilar | o mesmo teste geométrico dentro de `livre()` |
| shifting | `-∇C` só-fluido aponta para dentro do pilar | pilar como fonte do `ParticleShift` |
| métricas (`void_fraction`, R99, borda da figura) | pilar conta como vazio ou fura a colônia | máscara geométrica dos pilares em todas |

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

| rodada | o quê | duração | critério |
|---|---|---|---|
| **V0** | `use_pilares=False` contra o código atual | 600 passos | HDF5 e log bit-a-bit (10 threads) |
| **V1** | rede completa, só até o primeiro contato | t ≈ 20 (~5 min) | `n_pen = 0`; dt ≥ 0.5x o do liso no mesmo t; sem partícula ejetada |
| **R1** | rede completa, t = 50 | ~30 min | seção 8 |
| R2 | R1 com a rede girada 30° | ~30 min | separa efeito da rede de efeito de orientação |
| R3+ | uma alavanca por vez: Λ (12 e 26 dx), A (3 e 8 dx), depois campo aleatório | — | — |

Controle liso = o baseline escolhido (P2 ou P2S), mesma semente e threads. Como a morfologia
diverge caoticamente assim que algo muda (§2.2), toda comparação usa janela t ∈ [35, 50] e o
ruído entre as duas realizações do P2 como piso.

---

## 8. Predições e critérios para R1 (pré-registrar em `runs/<R1>/CRITERIOS.md`)

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
