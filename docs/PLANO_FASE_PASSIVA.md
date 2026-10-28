# Plano — física da fase passiva (filler): Marangoni, molhabilidade e influxo osmótico

> Escrito em 2026-10-28. **Nenhuma rodada deste plano roda antes de terminar a série de rugosidade**
> (R1″, R2, R3, R4, padrão sorteado — `docs/CRITERIOS_RUGOSIDADE.md`), porque todas as alavancas
> abaixo mexem no modelo liso e quebrariam o pareamento com o P2R23.
> Origem: a análise do modelo do professor (nota Obsidian *Modelo do professor (1 partícula = 1 célula)
> × nosso modelo*) e a pergunta sobre molhabilidade no lugar do pin químico.

## 0. O problema

O filler é 94% da colônia (#57) e é a crítica mais provável de um revisor. Ele tem três defeitos
diferentes, e cada parte deste plano ataca um:

| defeito | parte do plano |
|---|---|
| é **inerte**: não sente Marangoni nem flagelo, só é empurrado por pressão, viscosidade e shifting | série F, passo F1 |
| quando ganha motor, **espalha** pela lateral e a colônia vira blob (M4, #67-M) | série F, passo F2 (molhabilidade) |
| é **fabricado**: nasce por regra geométrica (rastro, wake, inserção), não por física | série O (influxo osmótico) |

Na teoria de duas fases de [T2] Srinivasan, o filler ocupa o lugar da **fase passiva** (fluido): ela
entra do ágar por influxo osmótico, é puxada pela Marangoni ([T1]) e espalha conforme a molhabilidade
([T17]). Este plano dá ao filler essas três propriedades, uma de cada vez.

## 1. Medidas que fundamentam o plano (pós-processamento, sem rodada nova)

### 1.1 O pin químico não segura baías (P2R23, t = 25–50)

Segura **55–64 células vivas, todas no núcleo** (mediana r/R99 = 0.07), em volta do anel pinado por
`rho_b` ≥ 0.8. As baías são ágar e ficam paradas por inércia (pressão zero, nenhum motor; #58).
**Conclusão: a molhabilidade não substitui o pin químico** — agiriam em lugares diferentes.

### 1.2 A lateral do dedo hoje está parada por construção

Borda = partícula de colônia com ≥ 30% de ágar no suporte 2h; força normal líquida = γ_eff·v_n
(regime superamortecido).

| borda, t ≥ 35 | composição | força normal líquida p10 / p50 / p90 |
|---|---|---:|
| ponta (14–17 líderes) | vivas | 0.0 / +4.8 / +7.0 |
| lateral (430–1017) | **95% filler** | −0.004 / +0.003 / +0.024 |

### 1.3 Quando o filler ganha motor, a lateral espalha (M4, base C4)

| borda, t ≥ 35 | C4 (filler passivo) | M4 (filler sente Marangoni) |
|---|---:|---:|
| lateral: força normal líquida p50 / p90 | −0.002 / +0.04 | **+0.62 / +1.0** |
| lateral: fração avançando | 41% (deriva ≈ 0) | **87%**, a 0.17 dx por unidade de t |
| ponta: força normal líquida p10 / p50 | (pontas pinadas por fome no C4) | +1.65 / +6.8 |

### 1.4 Janela da molhabilidade no baseline de hoje

Marangoni que a lateral de filler sentiria no P2R23 (reconstruída com a fórmula do solver, KGC, filler
como fonte do gradiente porque hoje ele conduz `c_s`), t = 35–50:

| | p10 | p50 | p90 |
|---|---:|---:|---:|
| lateral de filler: componente normal para fora | −0.17 | +0.05 | +0.37 |
| ponta: força motriz total (Marangoni + flagelo) | 3.75 | 5.9 | 9.2 |

| retenção F_w | laterais retidas | pontas acima (só Marangoni) |
|---:|---:|---:|
| 0.2 | 78% | 88% |
| **0.5** | **95%** | 88% |
| **1.0** | **100%** | 85% |
| 2.0 | 100% | 59% |

**Janela: F_w ∈ [0.5, 1.0]**, uma ordem de grandeza entre o p90 da lateral (0.37) e o p10 da ponta
(3.75). Os 12% de pontas com Marangoni normal zero andam pelo flagelo e ficam acima pela força total.
A lateral seria empurrada **12× menos que no M4** (p50 0.05 contra 0.62).

### 1.5 De onde vem o filler e onde o influxo entraria (P2R23)

| t | rastro (ágar convertido) | wake | inserção | massa |
|---:|---:|---:|---:|---:|
| 29.0 | 3367 | 551 | 13 | 11.4 |
| 50.0 | **8250 (82%)** | 1509 (15%) | 297 (3%) | 29.2 |

Produção em t ∈ [25, 50]: **~290 partículas (~0.8 de massa) por unidade de t**, nascendo na frente.

φ = fração de vivas entre as partículas de colônia no suporte 2h (t = 50):

| r/R99 | 0–0.1 | 0.1–0.3 | 0.3–0.5 | 0.5–0.7 | 0.7–0.85 | líderes |
|---|---:|---:|---:|---:|---:|---:|
| φ p50 | 0.34 | 0.15 | 0.02 | 0.00 | 0.00 | 0.06 |

**O influxo de [T2] entraria no núcleo; o filler de hoje nasce na frente.** Para trocar um pelo outro,
a fase passiva teria de levar ~290 partículas por unidade de t por 1.5–3 unidades de comprimento, na
velocidade da frente (`dR/dt` = 0.076). Hoje o filler anda a |v| p50 = 4.5×10⁻⁵ (p90 5.6×10⁻⁴):
**1700× (p50) a 135× (p90) mais devagar.** Mesmo com Marangoni (F1, |a| p50 ≈ 0.4), com o arrasto
atual (γ_eff ≈ 67 em `rho_b` = 0.4) daria v ≈ 0.006 — ainda **~12× curto**.

## 2. Série F — o filler como fase passiva móvel

### F1 — filler sente Marangoni (alavanca única sobre o P2R23)

- **Mudança:** o filler entra como destino da `MarangoniForce` e como fonte do gradiente de `c_s`
  ([T1]: a Marangoni puxa o fluido; [T2]: a fase passiva é advectada). Hoje ele está fora pela lição
  #39, cuja razão era o `c_s` congelado do filler; no P2R23 o filler conduz e decai
  (`FILLER_CS_CONDUZ` = 1), então a razão não vale mais — **é o primeiro ponto a vigiar.**
- **Predição:**

  | grandeza | P2R23 | F1 |
  |---|---:|---:|
  | motor nas vivas (`a_mar_bio_p95`) | 0.73 | sobe (M4: +7% sobre o C4) |
  | `mean_v` | — | sobe (M4: 2.6×) |
  | largura do dedo a 0.6 R | 5.0 dx (EDT) | +0.4 dx (p50) a +3 dx (p90) até t = 50 |
  | classe morfológica | fingering | fingering (sem blob) |

- **Falsificadores:** blob (a reconstrução estática errou a escala); motor nas vivas cai (a #39
  voltou por outro caminho).

### F2 — retenção da linha de contato (molhabilidade) sobre o F1

- **Forma:** ancoragem da linha de contato (histerese do ângulo de contato; de Gennes 1985,
  *Rev. Mod. Phys.* 57, 827). Na borda, a componente normal para fora da força motriz vira
  max(F_n − F_w, 0); **nunca recuo**. Rejeitadas: (i) molhamento parcial como força de recuo — contrai
  a colônia, porque a pressão da EOS (~10⁻³) não equilibra (#66-C, P3: R99 4.62 → 1.93); (ii) forças
  de par fluido-substrato (Tartakovsky & Meakin 2005, *PRE* 72, 026227) — o substrato não é
  representado por partículas (o ágar é espaço não colonizado).
- **Borda no solver:** o gate de |∇rho_b| da Marangoni, ou a fração de ágar no suporte (Liu §6.5).
- **Valor:** F_w = 0.7 (meio da janela de 1.4).
- **Predição:** largura volta à do P2R23 (±0.5 dx); velocidade da ponta −10 a −20%; ganho de motor do
  F1 mantido; AR e número de dedos dentro do ruído do P2R23.
- **Falsificadores:** a lateral ainda alarga (a janela estática não vale na dinâmica); pontas param
  (F_w acima da força total delas).

## 3. Série O — retomar o influxo osmótico

### 3.1 O que já falhou e por quê

| tentativa | o que fez | por que falhou | o que muda agora |
|---|---|---|---|
| M-A.1 (2026-05-05) | campo de osmolito `c_o`, massa somada à partícula ∝ \|∇c_o\| | massa sem dt (×11), `c_o` saturou, influxo isotrópico no rim (#16, #17) | a fonte passa a ser a fórmula de [T2] em φ, não \|∇c_o\| |
| X5 (análise, #54) | influxo como fonte de **biomassa** | o influxo traz solvente, não células: dilui `rho_b` | aqui o influxo cria **fase passiva**, que é exatamente solvente — a objeção não se aplica |
| OSM1 (#87) | ágar vizinho de viva vira colônia à taxa ∝ densidade de vivas, sem deposição | (a) taxa sem saturação: viva parada absorve tudo à volta (seletividade por taxa ≠ por tempo); (b) sem wake/inserção, o fosso em volta do núcleo é vazio, sem partícula para converter | (a) a fórmula de [T2] satura: ao diluir as células, φ cai a φ₀ e o influxo para; (b) a deposição fica até a troca de orçamento (O3) |
| P13 (#79) | pressão no ágar | o meio pressurizado resiste à frente | não envolvido |

### 3.2 O bloqueio central (§1.5): transporte

A água entra onde há células (núcleo); o filler é necessário na frente. Sem uma fase passiva que
**escoe**, o influxo incha o núcleo (o "nucleus maturation" da #18 com outro nome). Por isso a série O
começa pelo arrasto da fase passiva, não pela fonte.

### O0 — diagnóstico offline (sem rodada)

1. Mapa de V₀(x) = Q₀·(φ/(1−φ) − φ₀/(1−φ₀))₊ no P2R23 para φ₀ ∈ {0.05, 0.1, 0.2}; φ limitado a
   φ_max = 0.8 (no núcleo φ chega a 0.79 e a fórmula diverge em φ → 1).
2. **Q₀ calibrado para repor a produção atual:** ∫V₀ dA = ~0.8 de massa por unidade de t (§1.5).
3. Tempo de trânsito núcleo → frente em função do arrasto da fase passiva γ_f, com a força medida
   no F1. **Critério de avanço:** existe γ_f que leve a água à frente em ≤ 1/`dR/dt` sem tornar o dt
   inviável (limite do arrasto explícito, dt < 2/γ, e CFL).

### O1 — arrasto próprio da fase passiva (sobre o F2)

- **Mudança:** `LinearDrag` com γ_f para o filler, separado do γ das células ([T2]: arrasto entre
  fases distinto do atrito célula-substrato).
- **Estimativa a confirmar no O0:** para v ≈ 0.075 com |a| ≈ 0.4, γ_f ≈ 5 (contra ~67 hoje).
- **Risco:** lateral mais móvel → a janela de F_w (medida em força, não em velocidade) continua válida,
  mas tem de ser re-medida; fosso em volta do núcleo pinado (#73) pode abrir.

### O2 — influxo como fonte de partículas da fase passiva, em paralelo com o rastro

- **Mudança:** cada partícula de colônia acumula volume a taxa V₀(φ); ao chegar a dx², nasce uma
  partícula de fase passiva adjacente (`is_filler` = 1, `is_wake` = 1 — **não pinada**, senão fica
  presa no núcleo; inserção na rede dx com guarda de proximidade, como a inserção C3, para evitar o
  empilhamento da #27).
- **O rastro continua.** O2 só mede: fluxo de influxo, para onde a água vai, se chega à frente.
- **Nutriente:** a partícula nova não consome `c_n` (é água). A #39/S5 mostrou que tirar o consumo do
  filler desestabilizou o crescimento; aqui o filler antigo continua consumindo — **vigiar `c_n` e massa.**
- **Falsificadores:** o núcleo incha (a água não escoa); a frente não recebe nada.

### O3 — troca de orçamento

`RASTRO_W` reduzido na proporção do que o influxo entrega à frente no O2; produção total de fase passiva
mantida em ~290 partículas por unidade de t. Critérios: as réguas de forma (#95, #99), C5 (§2.2),
fração do filler nascido do influxo.

### O4 — sem rastro

Só se O3 mantiver a classe morfológica. É o objetivo: o filler deixa de ser fabricado.

## 4. Ordem e regras

1. Terminar a série de rugosidade.
2. F1 → F2 → O0 (offline) → O1 → O2 → O3 → O4, **uma alavanca por rodada** (§2.3 do CLAUDE.md).
3. Antes de cada rodada: teste de efeito curto (#102) e pré-registro neste documento.
4. Comparação sempre contra a rodada anterior da série, com o mesmo modo de execução.

## 5. Ligação com o modelo do professor

O O2 usa a ideia de nascimento de partículas do modelo dele (Cesgrav 2026, §3.2), mas para a fase
**passiva** (água), não para células — o que evita os problemas da mitose (#81, #85, #86). E a lei de
mobilidade da fase passiva perto dos pilares (O1) pode vir dele: a mobilidade de Kozeny–Carman
(§2.5 dele) medida na escala do poro e homogeneizada para a nossa escala.

## 6. Como não repetir o que já falhou — mapa de precedentes

Cada passo do plano tem um parente que já reprovou. Para cada um: o que reprovou, o que é diferente
agora (medido, não argumentado), o **portão offline** que tem de passar antes de rodar, e o **sinal
precoce** que aborta a rodada.

| passo | precedente que falhou | mecanismo da falha | o que é diferente agora | portão offline (antes de rodar) | sinal precoce de aborto |
|---|---|---|---|---|---|
| F1 | **M4** (#67-I/M) | filler com motor espalha pela lateral → blob | lateral empurrada 12× menos (p50 +0.05 contra +0.62, §1.4), porque o filler já conduz `c_s` | reconstruir no P2R23 a força normal da lateral **e** quanto ágar de baía fica a < 2h do filler deslocado (ver linha "laço #90") | largura a 0.6 R crescendo > 1 dx por 5 unidades de t; amplitude < 0.8× P2R23 |
| F1 | **#39** (filler no gradiente matou o motor) | `c_s` congelado do filler como fonte achatava ∇`c_s` | o filler conduz e decai (`FILLER_CS_CONDUZ` = 1) | reconstruir `a_mar` nos **líderes** com o filler como fonte vs. sem: diferença < 10% | `a_mar_bio_p95` < 0.6 (P2R23: 0.73) |
| F1, O1 | **laço #90** (disco tardio) | filler é doador da colonização (`COL_FILLER_DONOR` = 1, `rho_b` = 0.4): onde ele chega, recruta o ágar vizinho; filler que entra na baía a fecha | — **ainda não medido** | fração do ágar de baía (0.5–0.7 R99) que passa a ter filler no suporte se a lateral andar +0.4 / +3 dx | agar na faixa 0.5–0.7 R99 e limbo novo por quadro (os sinais da #104) |
| F2 | **#86-G** (frear a frente apagou o relevo) | frente mais lenta que a proliferação perde a selecao | a retenção age na lateral (força ~0), não no motor; na ponta freia 10–20%, contra os 50% do B1-motor | velocidade da ponta prevista com a força total medida, por percentil | amplitude < 0.8× F1 |
| F2 | **#66-C / P3** (coesão em região sub-densa contraiu a colônia) | força para dentro sem pressão que equilibre | retenção nunca empurra para dentro: só corta a parte para fora | — (por construção) | `R99` caindo entre quadros |
| O1 | **B2** (#73: ágar com arrasto baixo) | meio móvel é arrastado junto em vez de sair da frente | o arrasto baixo é só da fase passiva da **colônia**, não do ágar | dt explícito com γ_f e a velocidade prevista; pares < 0.05 dx no P2R23 por classe | pares < 0.05 dx (`tools/pares_dt.py`) acima do P2R23 no mesmo t |
| O1 | **cascata de pareamento** (§88: 98% dos pares são filler) | filler denso sem repulsão efetiva pareia e derruba o dt | — **risco aumenta** com filler mais móvel | idem | iter/t > 1.5× o P2R23 |
| O2 | **inserção C3.2** (#34: runaway) | partícula inserida que cresce e produz `c_s` realimenta a inserção | a partícula nova é água: não cresce, não produz `c_s` | — (por construção) | massa: razão 2ª/1ª metade > 1.3 |
| O2 | **#57 / #61** (`rho_b` com três papéis) | se a água nascer com `rho_b` > 0 ela vira doadora e "colônia" nas métricas; se nascer com 0, tem pressão zero e não escoa (#53) | **pré-requisito novo:** separar a densidade mecânica da biológica (`phi_s` da #61) antes do O2 | a água precisa de pressão **sem** ser doadora: verificar no código que nenhum gate biológico a enxerga | — |
| O2 | **OSM1** (#87) | taxa sem saturação; fosso sem partícula | saturação em φ₀ pela fórmula de [T2]; deposição mantida | O0: com o Q₀ calibrado, a água criada no núcleo por unidade de t não pode exceder o que O1 transporta | φ no núcleo caindo abaixo de φ₀ (núcleo diluindo) ou densidade no núcleo subindo (inchando) |
| O2 | **S5** (#39: filler sem consumo de nutriente desestabilizou) | o consumo do filler freava o crescimento | o filler antigo continua consumindo; só a água nova não consome | fração da massa de fase passiva que é água nova ao fim do O2 | `c_n` nos braços > 0.75 com massa acelerando |
| O3/O4 | **I3/I4** (#85: cortar deposição) e **OSM0** | sem deposição, a viva cava túneis vazios; `R99` cai | a deposição só sai na proporção do que o influxo entrega à frente (medido no O2) | O2 tem de mostrar água chegando a > 0.7 R99 | `R99` < 0.9× da rodada anterior |
| todos | **#76** (métrica interna premiou o blob) | continuidade e uniformidade são máximas num disco | — | toda métrica de qualidade interna reportada **junto** de dedos, AR, `R99` | — |
| todos | **#88** (ruído entre realizações) | uma rodada só, diferença dentro do ruído | — | piso de ruído: P2R23 em outra realização | — |

### Regras de processo que saem do mapa

1. **Portão offline primeiro.** Nenhuma rodada sem a reconstrução da coluna "portão" sobre quadros
   existentes. Foi isso que, nesta sessão, mostrou que a molhabilidade não substitui o pin e que o
   influxo entra no lugar errado — duas rodadas economizadas.
2. **Teste de efeito curto** (#102) antes de pré-registrar: a alavanca move algum bit em ~600 iterações?
3. **Sinal precoce com critério numérico** em cada pré-registro, checado a cada quadro, para matar a
   rodada cedo em vez de esperar t = 50.
4. **Uma alavanca por rodada**, comparada à rodada imediatamente anterior da série, no mesmo modo de
   execução (1 thread, como o R1).
5. **Pré-requisito estrutural antes do O2:** separar densidade mecânica e biológica (#61). Sem isso, a
   série O cai na #57 qualquer que seja a fonte.

## 7. Portões offline EXECUTADOS (2026-10-28, sobre o P2R23, t ≥ 35; nenhuma rodada)

### 7.1 Portão da #39 — o filler como fonte do gradiente de `c_s` (F1)

Marangoni reconstruída com a fórmula do solver (gate de |∇`rho_b`|, KGC), nas duas versões da fonte.
Validação: a versão atual reproduz a `au_mar` gravada (correlação 0.9999, razão 1.0001).

| | sem filler na fonte (hoje) | com filler na fonte (F1) | |
|---|---:|---:|---|
| **líderes** (n = 92): \|a\| p50 / p95 | 4.58 / 9.47 | **3.51 / 9.45** | razão p50 **0.86**; direção igual (cos 1.000) |
| **vivas do corpo** (n = 2802): \|a\| p50 / p95 | 0.081 / 0.51 | **0.291 / 1.59** | **3.3×** |
| vivas do corpo: componente para fora p50 | +0.068 | **+0.262** | para fora em 92% (94% hoje) |

**Veredito:** o critério pré-registrado (diferença < 10% nos líderes) **reprova por pouco: −14%**.
O p95 e a direção não mudam. O efeito grande é outro e é favorável: **o corpo da colônia ganha motor
(força para fora ×3.9)**, que é o que a lição #86 identificou como ausente ("o meio anda devagar por
falta de motor"). Predição do F1 revista: motor nos líderes −10 a −15% (não "sobe"); motor no corpo
×3–4. **O critério de aborto do F1 passa a ser `a_mar_bio_p95` < 0.55**, não 0.6.

### 7.2 Portão do laço da #90 — exposição do ágar de baía ao filler

| t | ágar de baía (0.5–0.7 R99) | já a < 2h do filler | novo se a lateral andar +0.4 dx / +3 dx | `c_s` do ágar novo: > 0.30 / > 0.25 / p99 |
|---:|---:|---:|---:|---|
| 35.8 | 298 | 95% | 12 / 16 | 0.000 / 0.125 / 0.277 |
| 41.7 | 761 | 86% | 37 / 107 | 0.000 / 0.047 / 0.273 |
| 50.0 | 1810 | 75% | 85 / 406 | 0.000 / 0.020 / 0.267 |

**Veredito: a exposição não é o limitante** — 75–95% do ágar de baía já está ao alcance do filler.
O laço está fechado pelo **gate de `c_s`** da colonização: nenhum ágar de baía passa de 0.30, mas o
p99 fica em 0.27–0.29 (a margem zero do §2 do `CRITERIOS_RUGOSIDADE.md`). **O sinal precoce do F1 é o
`c_s` máximo do ágar a < 2h do corpo**, não a geometria: se o filler móvel carregar `c_s` para a
baía e esse máximo passar de 0.30, o laço do disco reabre.

### 7.3 O0 — influxo de van't Hoff, calibração e transporte

V₀ = Q₀·(φ/(1−φ) − φ₀/(1−φ₀))₊, φ = fração de vivas entre as partículas de colônia no suporte 2h,
limitado a 0.8. Q₀ calibrado para ∫V₀ dA = produção atual de fase passiva (290 partículas = 0.84 de
área por unidade de t).

| t | φ₀ | partículas-fonte | Q₀ | influxo em r < 0.3 R99 | 0.3–0.5 R99 | > 0.7 R99 |
|---:|---:|---:|---:|---:|---:|---:|
| 41.7 | 0.05 | 2575 | 0.45 | 66% | 32% | 1.9% |
| 41.7 | 0.10 | 2041 | 0.55 | 69% | 29% | 2.1% |
| 41.7 | 0.20 | 1135 | 0.85 | 75% | 22% | 2.7% |
| 50.0 | 0.05 | 3228 | 0.53 | 79% | 21% | 0.1% |
| 50.0 | 0.10 | 2374 | 0.70 | 84% | 16% | 0.0% |
| 50.0 | 0.20 | 879 | 1.20 | 93% | 7% | 0.0% |

**O influxo cai quase todo no núcleo** (66–93% dentro de 0.3 R99; < 3% além de 0.7 R99), qualquer
que seja φ₀. A distância média ponderada da fonte até 0.85 R99 é **2.1–2.75 (40–51 dx)**: mesmo
andando na velocidade da frente, a água levaria **28–36 unidades de t** para chegar — a janela inteira.

**O que o transporte exige, por continuidade.** Para entregar 0.84 de área por unidade de t através
do anel r ≈ 1.5, por onde a fase passiva só passa pelos ~14 dedos de ~6 dx (~4.5 de perímetro), a
velocidade de escoamento tem de ser **v ≈ 0.84/4.5 ≈ 0.19** — **2.5× a velocidade da frente**. Com a
força que o F1 daria ao corpo (|a| p50 ≈ 0.29), isso exige **γ_f ≈ 1.5**, contra ~70 hoje no filler
(γ·(1 + 1.5·0.4²)): **~50× menos arrasto**.

**Veredito do O0:** numericamente viável (o limite explícito do arrasto, dt < 2/γ_f ≈ 1.3, fica muito
acima do dt ≈ 0.016), mas **é uma mudança de regime**: a fase passiva passaria a escoar quase sem
atrito, 2.5× mais rápido que a frente. Consequências a vigiar no O1: espalhamento lateral (a janela
F_w da §1.4 foi medida em força e continua válida, mas a lateral móvel tem de ser re-medida),
pareamento (a cascata da §88 é de filler) e o fosso em volta do núcleo pinado (#73).

### 7.4 O que os portões mudam no plano

1. **F1 continua de pé**, com predição corrigida (líderes −10 a −15%, corpo ×3–4) e dois sinais
   precoces: `a_mar_bio_p95` < 0.55 e `c_s` máximo do ágar vizinho > 0.30.
2. **O O1 é o passo decisivo da série O**: sem fase passiva ~50× mais móvel, o influxo de [T2] incha o
   núcleo. Antes do O1, medir no F1 (quando rodar) a velocidade real do filler do corpo, para ancorar
   γ_f em dado dinâmico e não na estimativa estática acima.
