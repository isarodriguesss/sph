# CLAUDE.md

Este arquivo governa todas as interações do Claude Code com este repositório. As regras aqui descritas sao **obrigatorias** e tem precedencia sobre qualquer comportamento padrao.

---

## 1. Projeto e Referencial Teorico

**Tese de Doutorado:** *"Simulacao do Swarming de Pseudomonas aeruginosa em superficies rugosas utilizando Hidrodinamica de Particulas Suavizadas (SPH)"* (Rodrigues, 2025).

Simulacao 2D em SPH (via [PySPH](https://pysph.readthedocs.io/)) da dinamica de biomassa bacteriana: crescimento logistico, forcas de Marangoni induzidas por gradientes de surfactante (ramnolipideo), fluxo viscoso, arrasto diferenciado e coesao via EOS. O objetivo final e reproduzir o padrao de *swarming* dendritico de *P. aeruginosa* (fingering de Mullins-Sekerka) observado experimentalmente (Michiels et al.).

### Objetivos de longo prazo (guiam decisoes arquiteturais)

1. **Superficies rugosas** *(em aberto — bloqueado ate atingir morfologia de `reference.jpg`)* — Partículas de contorno estaticas com geometria irregular que interagem mecanicamente com o fluido e alteram os campos de difusao/escoamento. Toda decisao de refatoracao deve preservar a capacidade de substituir as paredes planas atuais por topografias arbitrarias. **PRE-REQUISITO OBRIGATORIO (§2.2):** a IA so pode sugerir implementar rugosidade (Pass L) apos a simulacao reproduzir a morfologia dendritica de `reference.jpg` via mecanismos hidrodinamicos puros (Marangoni + Flagelar + EOS). Rugosidade e **refinamento fisico**, nao muleta para compensar motor insuficiente.
2. **Pressao osmotica (van't Hoff)** *(em aberto)* — Reativar `OsmoticForce` com formulacao termodinamicamente consistente (`Pi = iCRT`) acoplada ao campo de biomassa, substituindo a abordagem atual via ramo atrativo da EOS quando estável.
3. **Motilidade flagelar orientada por gradiente** *(✅ implementado — Pass K)* — Forca propulsiva de **magnitude constante** (`f0 * gate`) alinhada a `-∇cs` (quimiotaxia para agar fresco), gateada a swarmers de borda (`rho_b ∈[0.1, 0.6]`). Ver `FlagellarForce` em [src/equations.py](src/equations.py) e §3.2 Frente 5.

**Regra de Alinhamento Teorico (Rodrigues, 2025):** Antes de sugerir qualquer mudanca no codigo, a IA **deve** verificar se a alteracao e consistente com as equacoes governantes do projeto (SPH, Navier-Stokes incompressivel com termos-fonte biologicos) e com os tres objetivos de longo prazo acima. Mudancas que comprometam a extensibilidade para superficies rugosas, pressao osmotica ou motilidade flagelar devem ser sinalizadas e justificadas.

---

## 2. Regra de Validacao Baseada em Dados

A IA **nunca** deve sugerir mudancas cegas em parametros fisicos. Antes de calibrar qualquer variavel (`beta`, `sigma`, `gamma`, `D`, `lambda_`, `r_growth`, `c0`, `alpha_mon`, `f0`, etc.), o seguinte protocolo e **obrigatorio**:

### 2.1 Analise quantitativa do `log.csv`

- Verificar tendencias temporais de `mean_v`, `n_fast`, `a_marangoni`, `a_flag`, `mean_cs`, `max_cs`, `contrast_cs`.
- Identificar colapsos (motor de Marangoni morrendo), saturacoes (`cs` uniforme), e outliers (`n_fast = 1` indica particula rogue, nao expansao real).
- Comparar valores atuais com o **orcamento de aceleracoes** (Secao 8) antes de propor ajustes.

### 2.2 Validacao visual da morfologia

- Solicitar ou levar em consideracao a analise dos frames gerados em `main_output/movie/` para confirmar se a morfologia e fisicamente valida.
- **Referencia visual obrigatoria 1 (experimental):** [reference.jpg](reference.jpg) (Michiels et al.) — colonia com ~15-20 dendritos radiais longos, `AR >= 1:5`, separados por agar limpo, nucleo central coeso com picos de surfactante nas pontas (painel C da referencia). Foto biologica de *P. aeruginosa* PA14 em agar swarming.
- **Referencia visual obrigatoria 2 (numerica):** [reference_result.png](reference_result.png) painel **(b) Fingering** (Trinschek, John, Thiele 2018, Soft Matter — T1 §3.0). Resultado de simulacao thin-film 2D com `W = 0.1` (alta wettability) e `Γ_max = 0.5` (alta producao de surfactante): colonia com **7-9 dedos finos** radiando de um nucleo compacto, baias estaticas entre dedos, campo de surfactante (Γ) extendendo-se alem da biomassa em um halo radial. **Esta e a meta morfologica primaria do projeto** porque e gerada por equacoes governantes da mesma classe que as nossas (Marangoni + wettability + producao bioativa de surfactante) — se Trinschek conseguiu com thin-film, devemos conseguir com SPH. O painel (a) Modulated mostra o estado intermediario (rugosidade radial sem dedos definidos), (c) Circular o estado de falha (motor fraco), (d) Arrested o de bloqueio.
- **Criterio de sucesso (Pass I-K):** aproximacao monotonica da morfologia de `reference.jpg` (PA14) e do painel (b) de `reference_result.png` (Trinschek) — dendritos cada vez mais finos, longos e separados, com halo de surfactante visivel alem da biomassa.
- **Criterio de falha:** expansao circular uniforme, anel oco, colapso para disco compacto, ou blob ameboide com muitos bumps curtos (AR ~1:2).
- Nao declarar um Pass como "bem-sucedido" baseando-se apenas em metricas escalares — a forma da colonia **comparada a `reference.jpg`** e o validador final.
- **Gatilho para Pass L (rugosidade):** so propor apos frames mostrarem dendritos de AR >= 1:5 separados por agar limpo. Ate la, refinar mecanismos hidrodinamicos (Pass I-K).

### 2.3 Predicao antes da acao

- Para cada mudanca de parametro, a IA deve apresentar uma **predicao quantitativa** do efeito esperado (ex: "reduzir `D` de 4e-3 para 2e-3 deve aumentar `|nabla cs|` em ~40% e reduzir `L_D` de 1.5h para 1.1h").
- Apos a simulacao, comparar predicao vs. resultado. Discrepancias > 2x devem ser investigadas antes do proximo passo.

### 2.4 CRITERIOS DUPLOS — Core Pinning Mecanico + Saturacao Quimica (revisado 2026-04-29, pos-K.19)

**Historico:** o "core pinning failure" diagnosticado em K.15 (2026-04-22) era um sintoma de **dois bloqueios fisicamente independentes**, originalmente confundidos como um. A sequencia K.16-K.19 demonstrou que resolver o bloqueio mecanico nao e suficiente — a morfologia dendritica de `reference.jpg` exige resolver **ambos**.

**Contexto biologico:** na literatura (Michiels et al., Tremblay et al., Kearns — P. aeruginosa PA14 em agar swarming), a colonia exibe **nucleo denso imovel** (matriz EPS madura) com **dendritos finos sendo a unica frente de expansao**. Interior e baias entre dendritos sao **estaticos**, e o campo de surfactante e **localizado nas pontas** (nao saturado uniformemente).

#### A) Bloqueio Mecanico (RESOLVIDO por K.17, 2026-04-28)

Nucleo deve ser estaticamente imovel para que pontas tenham trajetoria propria.

**Falha original (K.15):** colonia inteira inflava radialmente, dendritos eram bumps arrastados sem trajetoria propria.

**Solucao implementada:** hard pinning (`u = v = 0` em `rho_b >= 0.8`) em [src/scheme.py:44-49](src/scheme.py#L44-L49). Custom Euler step zera velocidade do nucleo a cada timestep, congelando matriz EPS madura.

**Criterios de aceitacao (todos satisfeitos por K.17):**
1. ✅ `a_pressure` rotina < 30 em t > 20s (era 15-30 com picos > 50 em K.15).
2. ✅ `mean_v` < 0.005 (era 0.008-0.011 em K.15).
3. ✅ `n_fast` em queda apos bootstrap (era 900-1300 sustentado em K.15).
4. ✅ Razao efetiva `v_core / v_tip` ≈ 0 — nucleo literalmente congelado.

**Lever permitidas (apenas refinamento se necessario):** `gamma_mature` (core drag adicional), `c0` (EOS stiffness), `edge_fade` geometria. Nao sao mais bloqueio principal.

#### B) Bloqueio Quimico (NOVO — descoberto K.18-K.19, 2026-04-28→29)

Mesmo com nucleo mecanicamente congelado, o campo `cs` **satura uniformemente** porque a producao na borda excede o sumidouro. Toda a interface ativa (rho_b ∈ [0.1, 0.8]) produz cs amplificado por `motile_boost · tip_boost`, gerando inflacao quimica que mata a selecao Mullins-Sekerka.

**Falha demonstrada (K.19, t=0-100s):** `mean_cs` cresceu monotonicamente 0.13 → 2.20 antes do motor morrer (n_fast → 0 em t = 80s). Frames mostram colonia esfericamente expandindo com bumps uniformes (~15-20 protuberancias AR ~1:1), sem dendritos longos.

**Fluxo de massa de surfactante (steady-state interno):**
```
Producao_interna ≈ σ · qs(rho_b) · (1.2-rho_b) · tip_boost · motile_boost ≈ 0.3-0.5/s (zona ativa)
Sumidouro_interno = (λ_int + k_consume·rho_b) · cs ≈ 0.65·cs em rho_b=1 (com k_consume=0.5)
Drenagem borda  = -D_ext · ∇cs|_fronteira (proporcional ao salto cs_int - cs_ext)
```

Se Producao > Sumidouro + Drenagem, `mean_cs` cresce sem limite → `contrast_cs` colapsa → motor morre.

**Metricas diagnosticas (saturacao quimica):**
- `mean_cs` crescendo monotonicamente alem de 0.5 indica afogamento iminente.
- `contrast_cs` em queda monotonica (lição §11): mesmo com `mean_v` aparentemente saudavel, motor real esta morrendo.
- `n_fast` subindo enquanto `max_v` degrada — ativacao bulk espuria (lição §9 #8).
- `max_cs / mean_cs < 4` — cs distribuido uniformemente, nao concentrado em pontas.

**Criterios de aceitacao adicionais (saturacao quimica):**
5. `mean_cs` plateau em [0.15, 0.45] sustentado em t > 30s (sem crescimento monotonico).
6. `contrast_cs` plateau > 100 sustentado em t > 30s (sem queda monotonica).
7. `max_cs / mean_cs > 4` — cs concentrado nas pontas, nao distribuido uniformemente.
8. Motor vivo em t > 100s (`n_fast > 50`, `mean_v > 0.001`).

**Alavancas para resolver bloqueio quimico (DESBLOQUEADAS pos-K.19):**
- `k_consume` em `SurfactantEquation.post_loop` ([equations.py:147](src/equations.py#L147)) — atualmente 0.5, considerar 1.5-2.5.
- Coeficiente do `motile_boost` (atualmente `1+50·|v|`) — reduzir para 10-20× localiza producao.
- Coeficiente do `tip_boost` (atualmente `1+3·|∇ρ_b|`) — reduzir para diminuir amplificacao na zona de transicao.
- `σ` em conjunto com mudancas de sumidouro (nunca isolado — perde atribuicao).

#### C) Bloqueio Geometrico/Mullins-Sekerka (PROVAVELMENTE NECESSARIO apos B)

Mesmo com producao quimica controlada, **focalizacao geometrica** entre pontas requer que o campo cs no agar tenha alcance lateral comparavel a distancia tip-tip:
```
L_D_ext = √(D_ext / λ_ext)
τ_decay_ext = 1/λ_ext
τ_propagate_tip-tip ≈ d_tip-tip² / D_ext
```
Para competicao Mullins-Sekerka: `τ_decay_ext ≥ τ_propagate_tip-tip`, equivalente a `L_D_ext ≥ d_tip-tip`.

**Lições K.18-K.19 (críticas para próximas calibracoes):**
- `D_ext` e `λ_ext` agem em direcoes opostas para drenagem na borda. **Mexer simultaneamente perde atribuicao e quebra Pass J.**
- Reduzir `λ_ext` abaixo de `λ_int` enquanto motile_boost ativo destroi drenagem na fronteira agar-biofilme. **K.19 fez isso e mean_cs explodiu para 2.2 (vs 0.71 em K.18 com λ_ext = 5·λ_int).**
- A janela funcional `D_ext ∈ [0.008, 0.012]` da lição #5 era valida sob K.13 (motor afogando). **Pos-K.17 com bloqueio quimico resolvido, janela precisa ser re-medida.**

#### D) BLOQUEIOS REMANESCENTES

- **Pass L (rugosidade) permanece BLOQUEADO** ate frames mostrarem dendritos AR ≥ 1:5 com baias estacionarias e tip-splitting visivel.
- **Combinacoes simultaneas de parametros de surfactante proibidas** (ex: mexer em `σ` e `k_consume` na mesma rodada). Uma alavanca por vez para preservar atribuicao.
- **Reverter Pass J inadvertidamente** (reducir `λ_ext` para igualar `λ_int` sem entender que isso destroi drenagem na borda) e **proibido**. Lição K.19.

#### E) Caminho prescrito (Passes K.20+)

1. **K.20 — aumentar `k_consume`:** 0.5 → 2.0 (4×). Mantem Pass J intacto (`λ_ext_ratio = 5.0`). Predicao: `mean_cs` plateau em ~0.17, `contrast_cs` > 100 sustentado, motor vivo em t > 100s. **Validar criterios #5, #6, #7, #8.**

2. **Se K.20 estabilizar mas dendritos curtos:** atacar bloqueio C com cuidado de uma alavanca por vez:
   - K.21: `D_ext` 0.04 → 0.08 (mantem `λ_ext` alto). `L_D_ext` = √(0.08/0.75) = 0.327. Drenagem preservada, alcance lateral +41%.
   - Criterio: contrast_cs > 150, baias visiveis nos frames.

3. **Se K.20 + K.21 produzir baias mas pontas curtas:** ajustar β/f0 para regime Mullins-Sekerka (Etapa 2 — capilaridade vs flux).

4. **Se K.20 nao estabilizar mean_cs em < 0.45:** alavanca alternativa — reduzir `motile_boost` factor 50 → 20 ou desativar amplificacao por velocidade.

---

## 3. Skill: SPH Bacterial Swarming Specialist

### 3.0 Fundacao teorica — literatura referencial

Cinco artigos formam a base teorica do projeto. A leitura cruzada destes papers foi consolidada em 2026-05-02 e revela tensoes importantes com nossas hipoteses anteriores.

**[T1] Trinschek, John, Thiele 2018 — *Soft Matter* 14, 4464.** Modelo thin-film 2D com surfactante insoluvel + wettability + crescimento bioativo. Reproduz 4 morfologias (arrested, circular, modulated, fingering) controladas pelos parametros `W` (wettability) e `Γmax` (concentracao maxima de surfactante). Mecanismo de tendrils: gradiente forte de Γ nos tips + Γ saturado nas baias → baias arrested por wettability + Marangoni nos tips. Inspiracao original do projeto. **Painel (b) Fingering** (`W=0.1, Γmax=0.5`) — ver [reference_result.png](reference_result.png) — e a **meta morfologica primaria** do projeto: 7-9 dedos finos com `AR >= 1:5`, baias estaticas entre dedos, halo de Γ extendendo alem da biomassa. Suas equacoes governantes (Navier-Stokes thin-film + Marangoni + reacao-difusao de surfactante) sao da mesma classe que as nossas SPH; se Trinschek atingiu (b) com thin-film, devemos conseguir reproduzir com SPH. **Limitacoes:** modelo passivo, sem QS, sem flagelo, sem nutriente explicito.

**[T2] Srinivasan, Kaplan, Mahadevan 2019 — *eLife* 8, e42697.** Teoria multifase generalizada que unifica swarms e biofilmes via duas fases (ativa + passiva fluida) e equacoes de balanco de massa/momento. **Insight crucial:** swarming e biofilm sao regimes distintos:
- **Swarming (nutrient-rich):** `c ≈ c0` constante, capilaridade dominada, steady-state com velocidade `V = C₁·g₀·H·Ca^(1/3)`. Mecanismo: osmolytes secretados pelas bacterias → pressao osmotica → influxo de fluido do agar via van't Hoff `V₀(x) = Q₀·(φ/(1-φ) - φ₀/(1-φ₀))`.
- **Biofilm (nutrient-limited):** `Γ/c₀ ~ O(1)`, transient, dirigido por stress osmotico de polimero EPS via Flory-Huggins.

Validado experimentalmente em B. subtilis (este trabalho), E. coli (Wu & Berg, Ping et al.), V. cholerae (Yan et al.). **Implicacao critica para nosso projeto:** o modelo Mimura-Murray de consumo de nutriente (`dc/dt = -k·ρ_b·c`) descreve regime de **biofilme**, NAO swarming. Para swarming de P. aeruginosa em CA (casamino acids), nutrientes sao abundantes e nao-limitantes.

**[T3] Giverso, Verani, Ciarletta 2016 — *Biomech. Model. Mechanobiol.* 15, 643.** Modelo continuum 2D com sharp interface comparando expansao volumétrica (`Γ = K_γ·ρ·n`) versus chemotactic (`m = χρ∇n`). Linear stability analysis mostra dispersion curves com wavenumber caracteristico. Volumetrico → instabilidade k=1 (assimetria, translacao do centro de massa). Chemotactic → padroes mais simetricos com multiplos dendritos. Fingers crescem com `t^0.45 ≈ √t` (diffusion-limited). Validacao numerica via finite element no FreeFem++.

**[T4] Potomkin, Tournus, Berlyand, Aranson 2017 — *J. R. Soc. Interface* 14, 20161031.** Modelo individual de microswimmer com flagelo flexivel em fluxo de cisalhamento. Resultados-chave: (a) **flagella bending reduz viscosidade efetiva** em suspensoes diluidas SEM tumbling (vs Haines et al. que requeria tumbling); (b) flagela buckling assiste **escape de paredes**. Complexidade dependente da rigidez `K_b` do flagelo. Valida nossa abordagem de Frente 5 (motilidade flagelar) mas sugere extensao futura para incluir flexibilidade do flagelo.

**[T6] Liu & Liu 2003 — *Smoothed Particle Hydrodynamics: A Meshfree Particle Method* (World Scientific).** Texto-referencia obrigatorio para qualquer modificacao em consistencia SPH, instabilidade numerica, ou refinamento. Capitulos criticos para o projeto:
- **§3.3 — Particle Approximation Consistency:** condicoes de consistencia de ordem zero e um; requisitos minimos de vizinhos (~20 para gradiente confiavel, ~35 para Laplaciano).
- **§3.4 — Conservation Properties:** condicoes para conservacao exata de massa, momento linear e momento angular sob formulacao SPH antisimetrica.
- **§4.1-4.3 — SPH for Navier-Stokes:** formulacao da equacao de momentum, derivacao da pressao via EOS, papel da viscosidade artificial Monaghan.
- **§6.4 — Tensile Instability:** mecanismo de instabilidade sob pressao negativa (p<0 → kernel atrativo); criterios de Monaghan para diagnostico; mitigacoes (artificial stress, kernel gradient correction).
- **§6.5 — Treatment of Free Surface:** sub-amostragem do kernel em bordas e particulas isoladas; tradeoff acuracia vs estabilidade.

**[T7] Violeau 2012 — *Fluid Mechanics and the SPH Method: Theory and Applications* (Oxford U. Press).** Complemento ao Liu para fundamentacao teorica rigorosa. Capitulos criticos:
- **§3.4-3.6 — Discrete Consistency / Partition of Unity:** definicao formal de `σ_a = Σ_j V_j W_aj` e seu papel como medida direta do erro do operador SPH. Trigger v2.4 do Pass N baseado nesta secao.
- **§5.3 — Conservation in SPH:** prova rigorosa de que momento linear e exatamente conservado com formulacao antisimetrica; momento angular conservado apenas aproximadamente — refinamento parcial agrava drift.
- **§7.4-7.5 — Particle Refinement / Coarsening:** baseado em Vacondio 2013 e Feldman 2006; deriva razao otima `ε/α = 1`, padrao hexagonal 2D, requisitos de simetria para preservacao do CM e momento angular.

**Aplicabilidade obrigatoria (2026-05-28):** qualquer mudanca em (a) refinamento adaptativo (Pass N e variantes), (b) integradores customizados (CustomEulerStep, pinning), (c) novas equacoes SPH (forcas, difusao, EOS), (d) trigger ou criterio de qualidade do kernel, **deve citar explicitamente o capitulo/secao de Liu ou Violeau que fundamenta a alteracao**. Solucoes "tentativa e erro" sem ancoragem teorica nestas referencias sao **rejeitadas em revisao**. Ver §10 Proibicoes.

### 3.0.1 — Literatura numerica SPH para o problema do vacuo (adicionada 2026-06-02)

Tres familias de tecnicas atacam a perda de suporte de kernel (`σ_a < 1`, Violeau §3.6 / Liu §3.3) nos braços sub-resolvidos — o problema do "vacuo". Levantadas na analise comparativa de 2026-06-02 (ver §12 "Rotas para o vacuo"). **Pass N (splitting) e apenas UMA delas.**

**[T8] Particle Shifting Technique (PST) — Xu, Stansby & Laurence 2009 (*JCP* 228, 6703); Lind, Xu, Stansby & Rogers 2012 (*JCP* 231, 1499); Adami, Hu & Adams 2013 (*JCP* 241, 292 — transport-velocity).** Preenche vacuos REDISTRIBUINDO particulas existentes via deslocamento Fickiano `δr = -D∇C` (C = concentracao), sem criar particulas → **custo de dt ZERO**. Lind estendeu para superficie livre / kernel truncado — exatamente o caso Liu §6.5. Cura clumping/tensile instability (Liu §6.4). **Implementado como Rota A (2026-06-02), ver §12.**

**[T9] Timesteps individuais/locais — Springel 2005 (*MNRAS* 364, 1105, GADGET-2); Saitoh & Makino 2009 (*ApJL* 697, L99).** A resposta canonica ao problema "dt = min-h global" (lição #29): particulas de h pequeno integram no SEU passo curto (blocos hierarquicos), o bulk no passo grande. Saitoh & Makino provam que e preciso um **time-step limiter** (vizinhos nao podem ter dt muito diferentes) para nao corromper a fisica. Viabilizaria Pass N com α pequeno SEM os 85× de penalidade — mas exige reescrever o integrador (CustomEulerStep + pin). Rota B (futura).

**[T10] Kernel Gradient Correction (KGC) / CSPM — Bonet & Lok 1999 (*CMAME* 180, 97); Chen & Beraun / Liu & Liu (CSPM, dentro do [T6]).** Restaura consistencia de 1a ordem do operador de gradiente (reproduz ∇ de campo linear exatamente) MESMO com vizinhanca incompleta — corrige `∇cs` espurio nos braços a **custo de dt ZERO**, sem preencher o vacuo fisico. Ataca diretamente a justificativa cientifica (validade do gradiente de Marangoni). Aplicavel so aos operadores de gradiente (`MarangoniForce`, `FlagellarForce`, `BiomassGradient`). Rota C (futura).

**Distincao OBJETIVO ↔ MECANISMO (2026-06-02):** o objetivo e `σ_a ≈ 1` / operadores consistentes na zona ativa (obrigatorio — Liu §3.3). O mecanismo (splitting / shifting / KGC / timesteps individuais) e escolha de engenharia. Conflundir os dois foi o que custou 85× de dt na v2.4. **O vacuo e primariamente problema de DISTRIBUICAO (PST) e de OPERADOR (KGC), e so secundariamente de CONTAGEM (splitting).**

**[T5] Bru, Kasallis, Zhuo, Høyland-Kroghsbo, Siryaporn 2023 — *Biophys. Rev.* 4, 031305.** Review especifico de swarming em P. aeruginosa. Pontos cruciais:
- **Marangoni nao e dominante:** experimento de Yang et al. — adicionar surfactante (Triton X-100) DEVERIA reduzir gradiente de tensao superficial e enfraquecer Marangoni, mas EXPERIMENTALMENTE aumentou o swarming. Conclusao: pressao-osmotica (van't Hoff) e o motor dominante, NAO Marangoni.
- **Modelo multilayer:** bacteria + camada de surfactante + agar. Ramnolipidos produzem camada distinta da camada bacteriana (confirmado por IRIS imaging).
- **Henrichsen 1972 vs realidade P. aeruginosa:** Henrichsen definiu swarming via aggregates flagela-dependentes; sliding via expansao por crescimento sem flagelo. P. aeruginosa NAO forma aggregates organizados, sugerindo que swarming P. aeruginosa e combinacao **sliding + swarming**, nao swarming puro.
- **Diferenca PA14 vs PAO1 (MPAO1):** PA14 forma camada de surfactante + tendrils. MPAO1 nao produz surfactante na superficie → nao forma tendrils. Surfactante e necessario.
- **Papel do flagelo:** propor que flagelo aumenta producao de osmolytes (LPS, EPS) no rim, drenando fluido do agar. Hyperflagellation (Deforet et al.) causa hyperswarming via aumento da producao de osmolytes via metabolic turnover.
- **Knowledge gap:** mecanismo exato pelo qual flagelo causa influxo de fluido permanece nao caracterizado em P. aeruginosa.

### 3.1 Balanco critico de forcas

O swarming dendritico emerge do balanco entre:

| Papel | Forca | Mecanismo |
|-------|-------|-----------|
| **Desestabilizadora** | Marangoni | Gradiente de surfactante (ramnolipideo) gera tensao superficial diferencial que puxa a interface para fora |
| **Desestabilizadora** | Motilidade flagelar | Propulsao ativa de magnitude constante alinhada a `-∇cs` — sobrevive a saturacao do reservatorio |
| **Estabilizadora** | Viscosidade SPH | Difusao de momento suaviza perturbacoes de curto comprimento de onda |
| **Estabilizadora** | Arrasto diferenciado | Nucleo maduro (EPS) e imóvel; borda (swarmers) e motil |
| **Estabilizadora** | EOS coesiva | Pressao negativa em rarefacao mantem a colonia coesa |
| **Desestabilizadora** | Ruido estocastico | Perturbacoes na producao de surfactante semeiam instabilidades |

A condicao de fingering exige que as forcas desestabilizadoras superem as estabilizadoras **apenas nas pontas** (interface), enquanto o interior permanece coeso. Se elas dominam em todo lugar → explosao. Se nunca dominam → disco compacto.

### 3.2 Cinco frentes de modelagem fisica

**Frente 1 — Crescimento bacteriano logistico**
- `drho_b/dt = r * rho_b * (1 - rho_b/rho_max)` ([equations.py:16](src/equations.py#L16))
- Controla a velocidade da frente de saturacao. `r_growth` baixo (0.4, Pass K) estende o tempo de residencia na zona ativa (`rho_b ~ 0.5`), permitindo acumulo de surfactante e sustentando swarmers de borda (gate `rho_b ∈[0.1, 0.6]`).

**Frente 2 — Dinamica de reacao-difusao de surfactantes (Quorum Sensing / Hill)**
- Producao: `sigma * qs(rho_b) * (1.2 - rho_b) * noise`, com `qs = rho_b^2/(rho_b^2 + K^2)` ([equations.py:74-78](src/equations.py#L74-L78))
- Difusao bi-escala: Laplaciano de Brookshaw com `D_eff` variavel — `D_int` no biofilme, `D_ext` no agar, smoothstep em `rho_b ∈ [0.1, 0.5]` (Pass I.6/I.7)
- Decaimento espacialmente dependente (Pass J): `lambda_eff = lambda` no interior, `lambda_ext = 3*lambda` no agar
- O modelo Hill (`K^2 = 0.01`) implementa quorum sensing: producao so ativa para `rho_b > 0.1`. O fator `(1.2 - rho_b)` (Pass H) cria frente movel de producao com piso de 17% em `rho_b=1.0`.
- **Acoplamento c_n (Pass M-B.2):** `production *= c_n_factor` onde `c_n_factor = c_n/(c_n+0.1)` — gate de disponibilidade metabolica.
- **TENSAO COM XAVIER et al. 2011 (2026-05-15):** biologicamente, o operon rhlAB (responsavel pelo ramnolipideo em P. aeruginosa) e **inversamente** proporcional a taxa de crescimento via CCR — bacteria com c_n ALTO (crescimento rapido) suprime rhlAB; com c_n BAIXO (fase estacionaria) expressa rhlAB. A implementacao correta seria `c_n_factor_xavier = (1 - c_n)/(1 - c_n + 0.1)`. **POR QUE NAO INVERTEMOS:** o c_n_factor invertido concentraria producao de cs no interior profundo (c_n=0, rho_b=1) onde: (a) D_int=1.5e-3 confina cs (L_D_int=0.1 << raio do braco), (b) gate Marangoni e flagellar sao zero (rho_b≥0.8 → pinning K.17, grad_rho_b≈0). O motor morreria por starvation no rim. O c_n_factor atual (direto) produz cs na **zona de transicao do braco** (c_n intermediario + rho_b intermediario) — onde as forcas atuam — gerando o perfil cs correto para Marangoni outward, como confirmado pelo run sem motile_boost (2026-05-15). **Para artigo:** descrever c_n_factor como gate de "disponibilidade metabolica de substrato para biossintese de surfactante", nao como regulacao direta de rhlAB. Citar Xavier para a relacao inversa CCR como contexto biologico que o modelo aproxima via localizacao geometrica (zona de transicao, nao interior profundo).

**Frente 3 — Fluxo induzido por Marangoni**
- `F = -beta * nabla_cs * gate(|nabla_rho_b|)` ([equations.py:137-138](src/equations.py#L137-L138))
- Gradiente simetrico SPH: `cs_ij = cs_j - cs_i` (Pass E.1); negacao aplicada via `self.beta = -beta` no construtor
- Gate smoothstep em `grad_rho_b_mag in [0.05, 0.6]` restringe a forca a interface da colonia.
- Requer campo `cs` com camada de transicao resolvida em >= 1.5*h.

**Frente 4 — Pressao Osmotica (Equacao de van't Hoff)** *(a implementar — possivelmente unificada com Frente 6 pos-2026-05-02)*
- `Pi = iCRT` — pressao proporcional a concentracao de solutos intra/extracelulares.
- `OsmoticForce` existe em [equations.py:318-363](src/equations.py#L318-L363) mas esta desabilitada. A EOS coesiva atual (ramo atrativo) funciona como substituto simplificado.
- **Reformulacao 2026-05-02 (apos T2/T5):** o paper Srinivasan-Kaplan-Mahadevan (eLife 2019) e o review Bru et al. 2023 mostram que a pressao osmotica gerada por **osmolytes secretados pelas bacterias** (LPS, EPS, surfactantes) e o **motor dominante** do swarming, nao Marangoni. Mecanismo: bacterias em alta densidade secretam osmolytes → diferencial osmotico colonia-agar → influxo de fluido do agar via van't Hoff `V₀ = Q₀·(φ/(1-φ) - φ₀/(1-φ₀))` → expansao volumetrica steady-state. Isso unifica Frente 4 (osmotica) com Frente 6 (substrato/osmolyte field) numa unica frente.
- Meta atualizada: implementar campo `c_o` (osmolyte produzido pelas bacterias) com acoplamento ao influxo de massa SPH para substituir/complementar a EOS atrativa atual.

**Frente 5 — Motilidade Flagelar orientada por gradiente** *(✅ implementado — Pass K)*
- `F_flag = -f0 * gate(rho_b) * n̂(∇cs)` — magnitude **constante** (`f0 * gate`), direcao oposta a `∇cs`.
- **Escolha-chave:** magnitude nao depende de `|∇cs|` — por isso sobrevive a saturacao do reservatorio (quando `|∇cs|` cai mas direcao persiste). Diferente da Marangoni, que colapsa quando `cs` uniformiza.
- `motility_gate`: smoothstep em `rho_b ∈ [0.1, 0.6]` com pico em `rho_b = 0.35` — so swarmers de borda; nucleo maduro (`rho_b ≈ 1`) e agar livre (`rho_b ≈ 0`) sao imunes.
- Implementacao usa padrao loop/post_loop (como `BiomassGradient`): `loop` acumula `grad_cs_x/y` pelo gradiente SPH simetrico, `post_loop` normaliza o vetor total e aplica a forca com sinal `-f0` (mesma direcao da Marangoni — para agar fresco).
- `f0 = 3.0` (K.27); orcamento teorico `f0 ~ 0.5*gamma*v_term = 3`.
- **Invariante critico:** `FlagellarForce` precisa executar **apos** `SurfactantEquation` para usar `cs` atualizado. Ver ordem em [scheme.py](src/scheme.py).

**Frente 6 — Osmolyte/substrato (campo escalar)** *(a implementar — Pass M, REQUISITO ANTES DE PASS L)*

**REFORMULACAO 2026-05-02 (apos analise de T2/T5):** A proposta original de Pass M (substrato consumivel `c_n` a la Mimura-Murray) foi refinada apos leitura do Srinivasan 2019 e Bru 2023. Esses papers mostram que swarming de P. aeruginosa e regime **nutrient-rich** (`c ≈ c0` constante), nao nutrient-limited. O modelo Mimura-Murray descreve **biofilm**, nao swarming. A correcao do "halo de baia" em K.27 nao requer consumo de nutriente, mas sim representacao explicita dos osmolytes que dirigem o influxo de fluido do agar.

**Duas formulacoes possiveis para Pass M (a decidir empiricamente):**

**Pass M-A (osmolyte produzido):**
```
dc_o/dt = D_o · ∇²c_o + k_o · ρ_b · (1 - c_o/c_o_max) - λ_o · c_o
                       ↑ producao biomassa-dependente
                       ↑ saturacao para evitar runaway
```
Acoplamento: `c_o` modula a EOS atrativa (representando influxo de fluido)
ou adiciona source term na continuidade SPH para representar massa entrando do agar.

**Pass M-B (substrato consumivel — original):**
```
dc_n/dt = D_n · ∇²c_n - k_n · ρ_b · c_n
production_cs *= c_n/(c_n + K_n)   [Michaelis-Menten]
```
Mantem proposta Mimura-Murray. Pode ser apropriado se nosso regime experimental estiver na transicao swarm→biofilm (que ocorre em ageing colonies).

- Diagnosticado em K.27 (2026-05-01) como omissao estrutural do modelo. Sem este mecanismo, o modelo nao consegue suprimir as baias entre dendritos — particulas no rim recebem push radial uniforme, gerando halo isotropico sobreposto ao padrao dendritico.
- **Mecanismo biologico (Pass M-A, preferido pos-2026-05-02):** bacterias secretam LPS/EPS/surfactantes que funcionam como osmolytes. Bays cercadas por colonia que ja produziu osmolyte → gradiente local saturado, baixo influxo de fluido. Tips em agar virgem → gradiente forte, alto influxo de fluido por van't Hoff. Resultado: tips continuam expandindo via influxo, bays estagnam. Mass flow funnels do bulk para os tips via influxo localizado.
- **Mecanismo biologico (Pass M-B, alternativa):** bacterias consomem nutriente do agar. Bacterias na baia esgotaram o nutriente local; bacterias nos tips alcancaram agar virgem com nutriente fresco. Resultado: bays param de crescer, tips continuam. Apropriado se regime virar nutrient-limited.
- **Equacao a implementar:**
  ```
  dc_n/dt = D_n * Laplacian(c_n) - k_n * rho_b * c_n
  ```
  com `c_n(t=0) = 1.0` em todo o dominio (agar virgem uniforme), `D_n` pequeno (nutriente difunde devagar), `k_n` controla taxa de consumo.
- **Acoplamento com producao de cs:** modificar termo de producao em `SurfactantEquation.post_loop` para depender de `c_n`:
  ```python
  production = sigma * qs(rho_b) * (1.2 - rho_b) * c_n / (c_n + K_n) * tip_boost * motile_boost
  ```
  Onde `K_n` e constante de Michaelis-Menten — quando c_n esta esgotado (baias), producao cai a zero; quando c_n esta fresco (tips), producao saturada.
- **Acoplamento com crescimento (CRITICO — nao opcional):** `r_growth_eff = r_growth * c_n / (c_n + K_n)`. Bacterias na baia tambem param de crescer biomassa. **Diagnosticado em M-B.3 (2026-05-07): sem este acoplamento, swarmers das baias (c_n esgotado) continuam crescendo rho_b ate atingir 0.8 e entrar para o nucleo, expandindo-o. O c_n acoplado apenas a surfactante e insuficiente — o bloqueio e na biomassa, nao no surfactante.** Ver licao #18.
- **Resultado esperado:** baias suprimidas → halo radial desaparece → mass flui para tips → dendritos finos com agar limpo entre eles, matching reference.jpg.

### 3.3 Fundamentacao das direcoes de forca quimiotatica (FUNDAMENTAL — NAO VIOLAR)

Esta secao codifica uma confusao recorrente sobre a direcao das forcas `MarangoniForce` e `FlagellarForce`. Toda mudanca em qualquer fator de producao/sumidouro de cs (`c_n_factor`, `growth_headroom`, `sigma`, `tip_boost`, `motile_boost`, `qs`, `k_consume`, `lambda`) DEVE passar pela analise prescrita abaixo antes de ser proposta.

#### 3.3.1 Convencao de sinal — derivada diretamente do codigo

**`MarangoniForce.loop`** ([src/equations.py:217](src/equations.py#L217), [:254](src/equations.py#L254)):
```python
self.beta = -beta            # input beta=+1.0 → self.beta = -1.0
cs_ij = s_cs - d_cs          # cs_j - cs_d
acc_x = gate * self.beta * vol_j * cs_ij * DWIJ[0]
# = -gate × Σ_j vol_j × (cs_j - cs_d) × DWIJ
# = -gate × ∇cs|_d (forma SPH symmetric gather)
```
Portanto: **F_marangoni / m = −gate · ∇cs**

**`FlagellarForce.post_loop`** ([src/equations.py:662](src/equations.py#L662)):
```python
acc_x = -self.f0 * gate * gx / mag   # gx = ∇cs_x (acumulado no loop)
```
Portanto: **F_flag / m = −f0 · gate · n̂(∇cs)**

**Ambas as forcas apontam na direcao `−∇cs` — de regioes de ALTO cs para regioes de BAIXO cs.**

#### 3.3.2 Condicao para push outward (expansao da colonia)

Para que a borda da colonia seja empurrada para fora (em direcao ao agar), `−∇cs` deve apontar para fora nas particulas com biomassa. Equivalentemente: o **gradiente cs deve apontar para dentro**, ou seja, **cs deve DECRESCER monotonicamente para fora** ao longo de uma regiao com biomassa contigua.

```
Perfil cs radial OUTWARD-DRIVING (correto):     Perfil PEAK-AT-EDGE (patologico):

cs ┤██                                          cs ┤        ██
   │  ██                                            │      ██  ██
   │    ██  biomassa aqui                           │    ██      ██
   │      ██  recebe push outward                   │  ██          ██
   │        ██     ✓                                │██  push inward  ██  push outward
   │          ██___                                 │  em toda a       └────── (mas e agar
   └─────────────────── r                           │  biomassa antes           sem biomassa)
   interior   rim   agar                            └────────────────── r
   alto      baixo  ~0                              interior  pico   agar
                                                    (baixo) (alto)  (zero)
```

**O push outward acontece APENAS em particulas POSTERIORES (radialmente) ao pico de cs.** Particulas ANTERIORES ao pico (interior do pico) recebem push INWARD (compactacao).

#### 3.3.3 Onde o pico de cs cai em funcao do c_n_factor

A escolha do `c_n_factor` determina onde no perfil radial o cs e produzido e portanto onde o pico se forma:

| `c_n_factor` | Producao maxima em | Pico cs em | Particulas com push outward |
|---|---|---|---|
| `c_n/(c_n+0.1)` (direto) | c_n alto = borda/agar adjacente | mid-arm ou outer rim (depende de qs·(1.2-ρ)) | depende exato do pico — ver tabela 3.3.4 |
| `(1−c_n)/((1−c_n)+0.1)` (Xavier/CCR) | c_n baixo = interior depletado | interior profundo (rho_b≈1) | toda a colonia exceto o centro — **expansao uniforme isotropica** |

#### 3.3.4 Posicao do pico × morfologia esperada

| Pico de cs em zona | Biomassa com push outward | Biomassa com push inward | Morfologia |
|---|---|---|---|
| `rho_b ≈ 1.0` (interior profundo) | rho_b < 1 (~99% da colonia) | Apenas centro geometrico | Expansao isotropica uniforme, **sem seletividade dendrítica** |
| `rho_b ≈ 0.5` (mid-arm, swarmer zone) | rho_b < 0.5 (rim externo) | rho_b > 0.5 (compacta nucleo) | **Dendrítica seletiva** quando localizado azimutalmente por motile_boost |
| `rho_b ≈ 0.2-0.3` (outer rim) | Apenas particulas em rho_b < 0.2 ou agar | Maior parte da colonia | Predominio de push inward — colonia comprime |
| `rho_b ≈ 0` (agar virgem, externa) | Nenhuma (so agar sem biomassa) | TODA a colonia | **Patologica** — colonia colapsa |

**Caso M-B.10 (melhor morfologia ate o momento):** `c_n_factor` direto × `(1.2 − rho_b)` × motile_boost × qs(rho_b) → pico em mid-arm (rho_b ≈ 0.5) com localizacao azimutal pelos swarmers ativos. Outer rim (rho_b < 0.5) recebe push outward, nucleo (rho_b > 0.5) recebe push inward (compactacao do nucleo + extensao dos dendritos).

**Caso Xavier puro (atual, K.27→N=6):** `c_n_factor` invertido × `growth_headroom=1.0` × qs → pico em interior profundo (rho_b ≈ 1). Toda a colonia recebe push outward de magnitude similar — expansao isotropica sem seletividade (~25 bumps uniformes, `contrast_cs` em queda).

#### 3.3.5 Erro de raciocinio recorrente (CASO DOCUMENTADO PARA NAO REPETIR)

**Falacia:** "Os red spots de cs nas pontas observados em frames significam que cs concentrado na borda externa puxa a colonia para fora."

**Por que e falsa:** A forca calculada pelo codigo e `−∇cs`. Um pico de cs no rim externo empurra apenas o que esta IMEDIATAMENTE APOS o pico no sentido radial. Se este "apos" e agar virgem (sem biomassa), a forca atua em particulas SPH sem biomassa e nao impulsiona o swarming. Por outro lado, particulas ANTES do pico (a maior parte da colonia) recebem push INWARD.

**Como evitar:** sempre calcular cs_∞ (produção / decaimento total) em 4-5 zonas representativas e identificar onde o maximo cai ANTES de propor qualquer alteracao em fator de producao ou sumidouro de cs.

#### 3.3.6 Protocolo obrigatorio — analise pre-mudanca de fator de producao/sumidouro de cs

Antes de propor mudanca em `c_n_factor`, `growth_headroom`, `sigma`, `tip_boost`, `motile_boost`, `qs`, `k_consume` ou `lambda`:

1. **Calcular `cs_∞` em 4 zonas representativas:**
   - Interior profundo (rho_b=1, c_n=0)
   - Mid-arm (rho_b=0.5, c_n=0.4)
   - Outer rim (rho_b=0.2, c_n=0.8)
   - Agar (rho_b=0, c_n=1)

   `cs_∞ = (sigma × qs × growth_headroom × c_n_factor) / (lambda_eff + k_consume × rho_b)`

2. **Identificar a zona com maior `cs_∞`** (= localizacao do pico de cs no perfil radial).

3. **Verificar se ha biomassa contigua (rho_b > 0.1)** entre o pico e a fronteira com o agar.

4. **Validar o sinal do push esperado:**
   - Se a biomassa "apos o pico" e desprezivel (caso patologico) → REJEITAR a mudanca.
   - Se o pico cai no interior profundo e cs satura uniformemente → push outward uniforme, sem seletividade dendritica → mudanca pode ajudar o motor mas nao a morfologia.
   - Se o pico cai em mid-arm/swarmer zone E e localizavel azimutalmente (por motile_boost ou equivalente) → caminho para seletividade dendritica.

5. **Documentar a predicao**: incluir tabela `cs_∞` por zona antes e depois da mudanca, e marcar explicitamente em qual zona o pico se moveu.

**Esta analise e obrigatoria mesmo (especialmente) para mudancas "biologicamente motivadas".** A coerencia biologica nao implica push outward correto — a direcao do push e definida exclusivamente pelo perfil espacial de cs, nao pelo significado biologico do termo de producao.

### 3.4 Meta de superficies rugosas

Partículas de contorno (`solid` em [particles.py](src/particles.py)) atualmente formam paredes planas de 2 camadas. O objetivo e substituí-las por:
- Topografias irregulares (rugosidade controlada por amplitude e comprimento de onda).
- Interacao mecanica bidirecional: a rugosidade canaliza o fluxo e altera os caminhos de difusao.
- Condicoes de contorno de nao-deslizamento (no-slip) ou deslizamento parcial conforme o substrato.

Toda refatoracao da geometria de contorno ou do NNPS deve ser compativel com esta meta.

---

## 4. Comandos

```bash
make run          # Limpa output e executa simulacao (python main.py)
make view         # Visualiza output com PySPH viewer
make run_view     # Executa e visualiza
make format       # Lint e formatacao com Ruff
make paraview     # Exporta HDF5 para VTK (ParaView)
```

A simulacao gera arquivos HDF5 em `main_output/` a cada `print_freq` iteracoes.

---

## 5. Setup do Ambiente

```bash
conda create -n pysph_env python=3.10 numpy scipy matplotlib -c conda-forge
conda activate pysph_env
conda install pysph cython mako -c conda-forge
pip install PySPH
conda install mayavi -c conda-forge
conda install mpi4py -c conda-forge
```

---

## 6. Arquitetura

**Entry point:** [main.py](main.py) — define `SwarmApp(Application)`, parametros fisicos como globais, cria particulas, scheme e solver. O hook `post_step` imprime diagnosticos por iteracao e opcionalmente trata divisao celular (`use_splitting = False`).

**[src/particles.py](src/particles.py)** — `create_initial_state()` constroi a grade 2D: particulas fluidas com campo Gaussiano `rho_b_grown` (biomassa) + perturbacoes, e particulas-fantasma solidas (2 camadas) formando paredes.

**[src/scheme.py](src/scheme.py)** — `MyBiomassScheme(Scheme)` conecta todas as equacoes em dois `Group`s PySPH:
1. Pre-step (non-real): `SummationDensity` -> `BiomassEOS`
2. Main step: `MomentumEquation` (Monaghan `alpha`) -> `BiomassGrowth` -> `BiomassGradient` -> `MarangoniForce` -> `ViscousForce` -> `LinearDrag` -> `SurfactantEquation`

**Invariantes criticos:**
- `BiomassGradient` **deve** executar antes de `MarangoniForce` — Marangoni usa `grad_rho_b_mag` como gate de interface.
- `SurfactantEquation` **deve** executar antes de `FlagellarForce` — motilidade flagelar usa o campo `cs` atualizado via `grad_cs_x/y`.

`CustomEulerStep` estende `EulerStep` para integrar `rho_b_grown` e `cs`, com clamping: `rho_b_grown` em [0, 1], `cs` >= 1e-9.

**[src/equations.py](src/equations.py)** — Equacoes SPH customizadas:

| Equacao | Descricao | Referencia |
|---------|-----------|------------|
| `BiomassGrowth` | Crescimento logistico gateado por `rho_b < 0.8` e `c_n` (M-B.3): `rate = r * (1 - rho_b/rho_max) * c_n/(c_n+0.1)` | [equations.py:4-22](src/equations.py#L4-L22) |
| `BiomassGradient` | Gradiente SPH simetrico de biomassa; magnitude usada pelo gate Marangoni |[equations.py:21-49](src/equations.py#L21-L49) |
| `SurfactantEquation` | Reacao-difusao: Hill QS + frente movel `(1 - rho_b)` + Brookshaw + decaimento | [equations.py:52-83](src/equations.py#L52-L83) |
| `MarangoniForce` | `F = -beta * nabla_cs * gate(grad_rho_b)`, gradiente simetrico |[equations.py:86-150](src/equations.py#L86-L150) |
| `LinearDrag` | `F = -gamma_eff * v`, com `gamma_eff = gamma_base + gamma_mature * rho_b^2` |[equations.py:153-198](src/equations.py#L153-L198) |
| `ViscousForce` | Viscosidade SPH padrao com `mu_eff = mu * min(rho_avg, 1)` | [equations.py:215-259](src/equations.py#L215-L259) |
| `BiomassEOS` | "Soft Interior, Cohesive Edge" — repulsao quadratica + atracao leve + edge_fade |[equations.py:262-315](src/equations.py#L262-L315) |
| `OsmoticForce` | **Desabilitada.** Substituida pelo ramo atrativo da EOS. Mantida para referencia | [equations.py:318-363](src/equations.py#L318-L363) |
| `FlagellarForce` | **Pass K.** Motilidade quimiotactica de magnitude constante: acumula `grad_cs_x/y` em `loop`, aplica `f0 * gate(rho_b) * n̂(∇cs)` em `post_loop` | [equations.py:408-483] (src/equations.py#L408-L483) |

---

## 7. Parametros Calibrados (branch `ram-8827-v1`)

Calibracao pos-Passes A-I.7. **MARCO I.7:** Transicao blob→dendritico confirmada. `a_mar` pico **148**, selecao competitiva de dedos, tendril dominante formado. Difusao bi-escala (D_int/D_ext) e o mecanismo-chave para Mullins-Sekerka.

| Parametro | Variavel | Valor atual | Notas |
|-----------|----------|:-----------:|-------|
| Viscosidade | `mu` | **0.020** | K.23: I.2 revert parcial — coesão viscosa (I.2=0.012, pre-I.2=0.025) |
| Arrasto (base) | `gamma` | 60.0 | `a_drag` = `gamma*v_term` = 6 em v=0.1 |
| Arrasto (nucleo) | `gamma_mature` (scheme) | `1.5*gamma` | Pass I.4: razao core/edge 2.5x (era 0.3*gamma) |
| Coef. Marangoni | `beta` | **5.0** | Pass T2d: T1 10→5 — reduz tracao amplificada na fronteira pela metade; razao `|F_mar|/B_tension` ~3000× |
| Producao surfactante | `sigma` | **10.0** | Pass T2g: T2d 20→10 — desacelera saturacao cs; reduz mean_cs e tracao no rim. Baseline estavel ate t=100s sem fragmentacao (contrast_cs plateau ~12, a_pressure 3.0) |
| Saturacao cs (T1) | `cs_max` | **0.5** | Pass T1: substitui growth_headroom+tip_boost+motile_boost+c_n_factor+k_consume por `(1-cs/cs_max)` |
| Decaimento (agar) | `lambda_eff` (T1) | `0.5*lambda` = 0.075 | Pass T1: decay LENTO no agar (gera halo); `2*lambda`=0.30 no biofilme |
| Forca flagelar | `f0` | **3.0** | K.22: 0.5→3.0 — alvo teorico §8 (f0~γ·v_term/2=3); K.21 diagnosticou regime subcritico |
| Difusao (biofilme) | `D` (`D_int`) | 1.5e-3 | Pass I.3: gradiente afiado na interface (`L_D_int`=0.93h) |
| Difusao (agar) | `D_ext` | **0.08** | K.18 0.01→0.04; M-B.2 0.04→0.08; `L_D_ext=0.327`; habilita focalizacao Mullins-Sekerka |
| Decaimento | `lambda_` | 0.15 | Manter — confina `cs` mas permite penetracao de `L_D_ext` |
| Sumidouro enzimatico | `k_consume` | **1.0** | K.20 0.5→2.0; revisado M-B para 1.0; cs_∞(interior)≈0.10 |
| Taxa crescimento | `r_growth` | **0.02** | M-B.7: K.21 0.15 → 0.02 — elimina tip pumping (bloqueio E) |
| Modelo producao | `sigma*qs*(1.2-rho_b)*noise*tip_boost*motile_boost*c_n_factor` | — | K.6 adicionou `tip_boost`, K.7 adicionou `motile_boost`, M-B.2 adicionou `c_n_factor` |
| Visc. artificial Monaghan | `alpha_mon` | **0.12** | K.23: I.2 revert parcial — previne instabilidade de tração SPH (I.2=0.06, pre-I.2=0.15) |
| Vel. som (EOS) | `c0` | 0.8 | B ~ 0.09; tensao `tension_ratio=0.08` (M-B.9b) |
| Tensao superficial EOS | `tension_ratio` | **0.30** | Pass T2b: 0.20 → 0.30 — otimo local da alavanca (T2c=0.40 regrediu spikes +53/75% sem ganho morfologico); plateau de retornos atingido |
| Nutriente (consumo) | `k_n` | 0.5 | M-B.2: taxa consumo bacteriano de c_n |
| Nutriente (D agar) | `D_n` | 0.02 | M-B.2: difusao no agar livre |
| Nutriente (D biofilme) | `D_n_int` | **1e-4** | M-B.4: bi-escala — EPS bloqueia nutriente; `L_D_n_int=0.018` |
| Gate BiomassGrowth | `c_n` | smoothstep [0.4, 0.8] | M-B.6/M-B.8: cresce 0 em c_n<0.4, 1 em c_n>0.8 |
| Pinning quimico | `c_n` | hard cutoff < 0.4 | M-B.8 [scheme.py:58]; M-B.9a confirmou — pinning rigido protege coesao (lição §22) |
| Smoothing kernel | `h_factor` | 1.8*dx | ~35 vizinhos por particula |
| Timestep | `dt` | 5e-5 | Adaptivo, CFL=0.4 |
| Grade | `x_dim, y_dim` | **187x187** | M-B.10: 150→187 para preservar dx≈0.054 em dominio expandido |
| Dominio | `x/y_min/max` | **[-5, 5]^2** | M-B.10: 8x8 → 10x10 — evita colisao de dendritos com fronteira |
| Perturbacao inicial (rho_b) | `azimuthal_perturb` | `0.8*cos(8*theta)` | Pass K.12: N unificado em 8 (K.4 era N=5, K.11 N=10) |
| Perturbacao inicial (noise) | `noise` | `1.0 + 0.6*sin(8*theta) + 0.01*rand()` | Pass K.12: N=8 coerente com rho_b |

---

## 8. Orcamento de Aceleracoes (v_term = 0.1)

| Termo | Meta | Formula |
|-------|:----:|---------|
| `a_marangoni` (liquida) | ~6 | `beta * |nabla_cs| * gate` (so interface/pontas) |
| `a_flag` (Pass K) | ~2 | `gamma * v_term = 60*0.1` |
| `a_drag` | ~6 | `gamma * v_term = 60*0.1` |
| `a_pressao` EOS | <3 | `B * excess^2 * edge_fade` |
| `a_viscosa` | ~1.5 | `mu * v / h^2` |
| **Total |a|** | **5-15** | Equilibrio: Marangoni ~ Drag nas pontas |

---

## 9. Historico de Passes e Failure Modes

### Resolvidos

- **"Thermal death" (t~30s):** Producao de `cs` proporcional a `rho_b` colapsava em disco de baixa densidade. **Fix:** modelo Hill QS (satura producao para `rho_b > 0.3`).
- **"Diffusive death" (t~12s):** `D = 3e-4` muito pequeno. **Fix:** `D = 1.5e-2`, remocao do bias `(0.3 + grad_rho_b)`.
- **"Kernel asymmetry collapse" (t~11s, Pass E.1):** Gradiente SPH nao-simetrico capturava apenas deficit transiente de kernel. **Fix:** forma simetrica `cs_ij = cs_j - cs_i`. Resultado: `a_mar = 14.6` sustentado por ~3s.
- **"Diffusive flatness" (t~15s, Pass E.2):** `L_D = 2.9*h` excessivo. **Fix:** `D: 1.5e-2 -> 4e-3` (`L_D = 1.5*h`). `a_mar` subiu de 14.6 -> 16.28.
- **"Production-saturation lock" (t~15s, Pass F):** Hill saturava uniformemente, criando `cs` plano com `nabla_cs -> 0`. **Fix:** fator `(1 - rho_b)` criou frente movel de producao.
- **CSV column drift (Pass F):** Writer em[main.py:188-205](main.py#L188-L205) corrigido para 13 colunas. Revelou: outlier de particula unica (`n_fast=1`) poluia `max_v` e `a_drag`.
- **"Production-front extinction" (t~17s, Pass G):** Motor da frente movel muito fraco. **Fix:** `sigma: 0.4 -> 1.2` (x3) + `r_growth: 2.5 -> 0.8`. Resultado: `a_mar` pico 49, `n_fast` pico 657, morfologia dendritica confirmada.
- **Morfologia dendritica (Pass G):** Frame 023 (t~30s) mostra colonia r~1.7 com perimetro irregular e ~30 protuberancias — fingering de Mullins-Sekerka no regime dominado por ruido. Fase ativa 3-37s.
- **"Tail extinction" (t~37s, Pass H → resolvido por Pass J):** Colonia saturava globalmente, fator `(1 - rho_b)` zerava producao, `cs` drenava exponencialmente. Fix inicial (Pass H): `growth_headroom = 1.2 - rho_b` garantiu 17% de producao em `rho_b=1.0`. Insuficiente sozinho — drenagem interna uniformizava `cs` (`mean_cs` subia, `contrast_cs` caía de 146 → 13). Fix definitivo (Pass J): decaimento `lambda_ext = 3*lambda` no agar removeu drenagem interna. Resultado: `contrast_cs` estavel em 54, `mean_v` nao decai ate t=24.5s.
- **Pass J - Motor sustentado indefinido:** Decaimento espacialmente dependente em `SurfactantEquation.post_loop` (`lambda=0.15` interior, `lambda_ext=0.45` exterior, smoothstep em [0.1, 0.5]). `log.csv` confirma: `mean_v` estavel em 0.0075, `n_fast` 150-300, `max_cs` cresce monotonicamente ate 8.2. `a_mar` pico 182 (+23% vs I.7 isolado) — Pass K tambem contribui.

### Em andamento

- **Pass K.8 — "One-tip tyranny" (confirmada visualmente):** Apos K.5-K.8 (ver §12), frame 034 (t~28s) mostrou **um unico dendrito** radial com AR ≥ 1:10 e nucleo compacto preservado — morfologia individual **identica a `reference.jpg`**. Porem referencia pede **10-15 dendritos**, nao um. Causa: winner steals — um vencedor consome todo o gradiente de cs via difusao lateral, extinguindo competidores (`n_fast` cai de 670 → 130 entre t=3s e t=13s).
- **Pass K.9 — few-tip tyranny:** `sigma: 0.6→0.8` relaxou one-tip para **3-4 dendritos** (frame 039, t=30s). `a_mar` pico 300, `contrast_cs` 95-155, sustentado. Melhor que K.8 mas ainda sub-alvo.
- **Pass K.10 — reducao D_ext falhou:** `D_ext: 0.01→0.005` tentou localizar cs para prevenir winner-steal. Resultado: **pior morfologia** (5-6 bumps curtos, AR ~1:2). Lecao: reduzir D_ext mata o campo long-range necessario para focalizacao geometrica Mullins-Sekerka. **D_ext=0.01 e otimo**, nao reduzir mais.
- **Pass K.11 — multi-seeding N=10 (rho_b) + N=8 (noise):** `D_ext: 0.005→0.01`, `cos(10θ)` em particles. `a_mar` pico 368 (vs K.9 300, +23%), `n_fast` pico 560. Porem `contrast_cs` em queda monotonica **188→83** em 16s e `max_cs` declinou apos pico (51.5→33.4). **Regressao morfologica**: N=10 semeou muitos tips que se sobrepuseram via difusao lateral, homogeneizando o campo. Mismatch N=10/N=8 (rho_b vs noise) violou licao I.1 de coerencia de modos.
- **Pass K.12 — unificar em N=8:** `cos(8θ)` tambem em particles.py. `a_mar` pico identico (368), `contrast_cs` cai igualmente (204→81). **Unificacao morfologicamente neutra**. Dado novo: `n_fast` apresentou **pulsacao forte** (627→165→567 em 12s) nao vista em K.11. Frames revelaram blob lobular com ~14 bumps uniformes AR 1:2 — pre-semeadura nao resistiu a homogeneizacao difusiva dos plumes.
- **Pass K.13 — avanco parcial + falha terminal diagnosticada tardiamente (2026-04-22):** `lambda_ext_ratio=5.0` em [src/scheme.py:131-138](src/scheme.py#L131-L138). `lambda_ext=0.75`, `L_D_ext=0.115`. Frames 090-130 (t=28-40s) mostraram dendritos visualmente discerniveis (avanco genuino vs K.12 bumps AR 1:2), porem **morfologia incorreta vs P. aeruginosa**: lobulos bulbosos, SEM tip-splitting, SEM estrutura fractal. **Falha terminal confirmada via analise temporal extendida (t=0→53s):** `contrast_cs` **colapsa monotonicamente 131→32** entre t=30s e t=53s, `mean_cs` **cresce sem limite 0.26→1.17** (saturacao global do interior), `n_fast` sobe 430→1787 (4x) com `max_v` degradando — **ativacao bulk disfarcada de motor seletivo**. Causa raiz: **ausencia de sumidouro fisico** para cs. Com `(1.2-rho_b)=0.2` + `motile_boost` + `lambda_int=0.15`, estado estacionario interior cs_∞ = σ·0.17/λ ≈ 0.91 — afogamento inevitavel. Falha auto-diagnosticada apos frame 056 ser declarado sucesso prematuramente. **Pass L permanece BLOQUEADO.**
- **Pass K.14 — sumidouro biomassa-dependente (Etapa 1, 2026-04-22):** adicionado termo `-k_consume * rho_b * cs` em [src/equations.py:147](src/equations.py#L147) `SurfactantEquation.post_loop`. `k_consume=0.5` — decaimento efetivo interior = `lambda + k_consume*rho_b` = 0.65 em rho_b=1 (4.3x maior que K.13). Exterior (rho_b≈0) preserva mecanismo inalterado. **Fundamentacao biologica**: ramnolipideos sao adsorvidos nas membranas bacterianas e degradados por enzimas rhlE/rhlB reguladoras intrinsecas em alta densidade (QS homeostase). **Resultado**: motor prevenido de afogar mas sub-amplitude — `a_mar` pico 172 (vs K.13 368). Morfologia lobular sem dendritos. Evoluiu para K.15.
- **Pass K.15 — sigma 0.8→1.2 para restaurar amplitude (2026-04-22):** mantido `k_consume=0.5`, aumentado `sigma: 0.8 → 1.2` em[main.py:45](main.py#L45). Resultado (t=0→33.8s): **primeira morfologia dendritica real** (~15-18 dendritos, AR 1:5 a 1:7 em frames 097/110/132), motor sustentado (`a_mar` pico 329, pulsos 100-260), `contrast_cs` plateau 55-95 apos transitorio (inicial 204→52 em t=26s, recuperou para 94 em t=29s), `mean_cs` plateau 0.47-0.55 (afogamento estabilizado), `max_cs` 30-45. **Motor resolvido (Etapa 1 ✓).** Porem: **core pinning failure diagnosticado** (§2.4) — `a_pressure` 15-30 rotina com picos >50 (orcamento <3), `n_fast` 900-1300 (ativacao em massa, nao localizada nas pontas), razao drag core/borda apenas 2.2× (`gamma_mature=1.5·gamma`). Conclusao: colonia inteira infla radialmente, dendritos sao bumps arrastados. Tip-splitting ausente. **Pass K.16 abordara core pinning (§2.4).**

- **Pass K.16a-e — varredura de alavancas mecanicas falhada (2026-04-28):** testou todas as 3 alavancas §2.4 originais isoladamente. **K.16a** (`gamma_mature: 1.5x → 5x`): drag combate velocidade ja estabelecida, nao a fonte de pressao. **K.16c** (edge_fade invertido com fade_repulsao=fade_atracao=0 em `rho_b≥0.8`): SPH pair-wise pressure persiste; sem coesao no core, a casca de transicao 0.5-0.8 absorve toda carga de pressao. **K.16d** (`c0: 0.8 → 0.35`): SPH auto-regula via `excess²` quando B cai. **K.16e** (edge_fade assimetrico — `fade_repulsao=0`, `fade_atracao=1` em core): `B_tension = 0.02·B` e negligível (50× mais fraca que repulsao); manter atracao no core nao mudou dinamica. **Conclusao: alavancas mecanicas isoladas insuficientes para core pinning.**

- **Pass K.16f — `gamma_mature = 8x` parcialmente eficaz (2026-04-28):** `a_pressure` 8-22 rotina (vs 15-35 K.16e), `mean_v` 0.005-0.008 (-30%), `n_fast` 300-450 (-30%), `contrast_cs` plateau 122-180. **Primeiro avanco mensuravel da sequencia K.16.** Frames 060-072 mostraram 4-5 dendritos distribuidos angularmente. Porem extensao para t=40s: `mean_cs` saturou 0.43→0.65, `contrast_cs` colapsou 138→81, `n_fast` explodiu 300→972 — **recidiva K.13 syndrome**. Conclusao: drag pesado contem expansao mecanica mas **nao previne saturacao quimica** (motile_boost ainda ativo na borda).

- **Pass K.17 — hard pinning rigido (2026-04-28):** implementado `u = v = 0` para `rho_b ≥ 0.8` em [src/scheme.py:44-49](src/scheme.py#L44-L49) `CustomEulerStep.stage1`. Reverteu `gamma_mature` para baseline 1.5x. **Bloqueio mecanico (§2.4-A) RESOLVIDO**: `mean_v` 0.001-0.005, `n_fast` 200-340, nucleo congelado literalmente. Critério #4 `v_core/v_tip < 0.2` passa diretamente. Porem: **morfologia continuou anel uniforme** (~15-20 bumps AR ~1:1) em vez de dendritos. **Diagnostico**: pinning resolveu sintoma errado — bloqueio nao era cinematico (nucleo movendo) mas **quimico** (cs saturando uniformemente apesar de nucleo estatico).

- **Pass K.18 — `D_ext: 0.01 → 0.04` (2026-04-28):** hipotese: `L_D_ext` muito pequeno impede focalizacao Mullins-Sekerka. Calculo: `L_D_ext = √(D_ext/λ_ext) = √(0.04/0.75) = 0.231` (assumindo `λ_ext_ratio=5`); razao `L_D_ext/R = 0.15`. Predicoes: `contrast_cs` subir, `max_cs` subir, `mean_cs` cair. Resultado real: `mean_cs` 0.13→0.71 (cresceu 5×), `contrast_cs` continuou colapso. **Predicoes falharam.** Diagnostico: `D_ext` maior aumenta drenagem na borda mas nao toca producao interna; insuficiente para deter saturacao quimica. Frame 078 (t=52s): expansao mais contida (R=1.5 vs K.16f R=2.0) mas sem dendritos.

- **Pass K.19 — `λ_ext_ratio: 5.0 → 1.0` (2026-04-29) — FALHA CATASTROFICA:** hipotese: aumentar `τ_decay_ext` de 1.3s para 6.7s permite cs no agar alcancar pontas vizinhas (`L_D_ext = 0.516`, razao 0.34). Resultado completo (t=0-100s): `mean_cs` cresceu 0.13 → **2.20** (saturacao 17×, vs 0.71 em K.18). Motor morreu em t=80s (`n_fast → 0`, `mean_v → 0.0001`). Frame 345 (t=99s): colônia parou de crescer em R=2.5, painel cs cinza uniforme (cs decaído quase a zero), |a_flag| zero. **Causa raiz**: reduzir `λ_ext` para igualar `λ_int` permite cs acumular no agar → gradiente na fronteira agar-biofilme cai → fluxo difusivo que drenava cs do bulk **diminuiu drasticamente** → `mean_cs` interno explode. **Reverteu Pass J inadvertidamente.** Lição: a premissa "pontas roubam gradiente das baias" requer **bulk com cs bounded primeiro**; quando `mean_cs` satura para 2.0+, gradiente lateral tip-tip e irrelevante (eclipsado pelo salto vertical biofilme-agar). **Pos-K.19, §2.4 revisada para reconhecer bloqueio quimico como independente do mecanico.**

- **Pass K.20 — `k_consume: 0.5 → 2.0` (2026-04-29, REVERTEU K.19 + aplicou k_consume):** Bloqueio quimico §2.4-B resolvido. Steady-state interior cs_∞ = σ·0.17/2.15 ≈ 0.10 (vs 0.85 K.13). `mean_cs` plateau 0.08-0.18 ✅, contrast_cs > 80 sustentado. Porem motor sub-amplitude (mean_v 0.001-0.002, n_fast 0-46). Evoluiu para K.21.

- **Pass K.21 — `r_growth: 0.4 → 0.15` (2026-04-30):** Eliminado K.17 Trap (saturacao global de rho_b que congelava toda a colonia em t~68s com r_growth=0.4). Tempo de residencia swarmer estendido de ~60s para ~150s. Diagnostico pos-K.21: Blockers A e B ✅ confirmados. **"Motor Starvation" identificado** — pontas tentam avançar mas estagnam. Fronteira "ferve" sem elongacao. AR ~1:1 nas poucas protuberancias que sobrevivem. Ver lição #14 e §2.4-E para analise completa.

- **Pass K.22 — `f0: 0.5 → 3.0` (2026-04-30, FALHA — motor ativo mas material frágil):** f0=3.0 aumentou a_flag de 0.5→3.0 e restaurou ignição flagellar (a_flag=3.0 constante no log). Porem morfologia ficou **mais estagnada que K.21**. Causa raiz diagnosticada: **instabilidade de tração SPH** ("brittle neck"). Hard pinning (K.17) ancora o nucleo; f0=3.0 puxa a ponta com força 3.0; o "pescoço" do dendrito tem coesão EOS de B_tension = rho0·c0²/7·0.02 = 0.00035 — quase zero. F_tração/F_coesão ≈ 8500 → fratura antes de elongar. Adicionalmente, motile_boost=3.5× com f0=3.0 elevou produção de cs, mean_cs subiu 0.08→0.31 em 100s (Bloqueio B re-emergindo lentamente). **Diagnostico confirmado pelo usuario (2026-05-01):** Pass I.2 deixou o fluido excessivamente "fino" (mu=0.012, alpha_mon=0.06), criando material quebradiço em vez de viscoelástico. Solução: revert parcial I.2 para fortalecer coesão.

- **Pass K.23 — revert parcial I.2: `mu: 0.012→0.020`, `alpha_mon: 0.06→0.12` (2026-05-01):** Fortaleceu coesao SPH contra fratura. Resultado: 2 dendritos coerentes formaram mas estagnaram em ~25s por K.17 trap. Confirmou que `alpha_mon` e o lever direto contra instabilidade de tracao.

- **Pass K.24 — `r_growth: 0.15 → 0.05` (2026-05-01):** Estendeu swarmer pool ~3×. Diagnosticado snapping/anel voador devido a f0=3.0 violento + r=0.05. Frames mostraram 2 dendritos formando + ring de particulas voadoras descoladas do corpo.

- **Pass K.25a — `threshold: 0.8→0.95 + f0: 3.0→1.5` (2026-05-01, FALHA):** Borbulhamento sem estrutura. a_pressure caiu 3-5× (1.5-3.5 vs 8-22 K.23). Ancora destruida + motor sub-critico = colonia gel uniforme sem dendritos. **Confirmou: hard pinning rigido em 0.8 e estrutural — nao negociavel.**

- **Pass K.26 — `threshold: 0.95→0.8 mantendo f0=1.5` (2026-05-01, SUCESSO MORFOLOGICO):** **PRIMEIRA MORFOLOGIA DENDRITICA MULTIPLA PERSISTENTE** desde inicio do projeto. Frames 050-079 mostram **8 dendritos formados** (matching cos(8θ) seed), AR ~1:4-1:6, dendritos persistentes (sem swallowing), nucleo ancorado (sem ballooning), sem snapping. Limitacao: dendritos bulbosos (largura 4-5 particulas vs 1-2 da reference.jpg) + halo radial de baia visivel.

- **Pass K.27 — `f0: 1.5 → 3.0` (2026-05-01, DIAGNOSTICO ESTRUTURAL):** Manteve sucesso K.26 e amplificou ambos os fenomenos simultaneamente — os 8 dendritos cresceram (n_fast pico 113, +50%) E o halo radial das baias tambem cresceu. **Diagnostico fundamental:** o motile_boost atual (1+50·v) gera diferencial tip/bay de apenas 2.3× — insuficiente para criar gradiente azimutal de cs capaz de funilar massa das baias para os tips. Bays expandem em paralelo, nao param. **CONCLUSAO: calibracao parametrica esgotada para suprimir baias. Mass flow para os tips requer mecanismo de SUPRESSAO ATIVA das baias — nutriente consumivel `c_n` (Frente 6, Pass M).**

- **Pass M-A.1 — campo de osmolito `c_o` com influxo de massa via van't Hoff (2026-05-05, FALHA — borbulhamento isotropico):** primeira tentativa de Pass M-A (Frente 6 osmolito) implementada conforme `PASS_M_HANDOFF.md`. Adicionados: campo `c_o`, equacao `OsmolyteProduction` (difusao Brookshaw + producao biomassa-dependente Hill `k_o · qs · (1-c_o)` + decaimento `λ_o · c_o` + influxo `dm += Q0 · gate(rho_b) · |∇c_o|`). Parametros: `D_o=1e-3, k_o=0.5, λ_o=0.05, Q0=5e-4`. Inserida entre `SurfactantEquation` e `FlagellarForce` em [scheme.py:159-169](src/scheme.py#L159-L169), integrada em `CustomEulerStep` com clamp em [0,1]. Diagnosticos no log: 6 colunas novas (`min_c_o`, `max_c_o`, `mean_c_o`, `contrast_c_o`, `mass_total`, `dm_influx`).

  **Resultado morfologico (frames 000-062, t=0→30s):** ~30 bumps curtos e uniformes (AR ~1:2) borbulhando radialmente, **sem selecao competitiva**. Pior que K.26-K.27 (que tinham 8 dendritos persistentes). Um filamento fino atravessa NE em frames 050-062 — **artefato de ejecao numerica** apos pico de velocidade `max_v=1.34` em t=18.5s (vs ~0.5 K.26), nao dendrito real.

  **Tres falhas independentes diagnosticadas:**
  1. **Massa runaway (BUG DE INTEGRACAO):** `d_m += Q0 · gate · |∇c_o|` no `post_loop` **NAO foi multiplicado por dt**. Acumula a cada chamada do equation evaluator (1/dt = 20 000×/s com dt=5e-5). `mass_total` cresceu de 36.5 → 406.6 em 30s (**11×**), `dm_influx` = 370. Taxa observada ≈ 12.3 unidade-massa/s, consistente com `Q0·gate·|∇c_o|/dt`. Cada particula do rim incha → SPH pressure radial → expansao tipo balao.
  2. **Saturacao instantanea de c_o:** `max_c_o = 1.0` ja em t=2.86s (iteracao 400). Steady-state interior sob k_o/λ_o = 10:1 e `c_o_∞ = k_o·qs/(k_o·qs+λ_o) ≈ 0.91`. Apos saturacao, `(1-c_o) → 0` zera producao mas decaimento e lento — campo congela proximo de 1 dentro da colonia. `mean_c_o` cresce so de 0.006 → 0.015 em 30s (penetracao difusiva insignificante: `L_D_o = √(D_o/λ_o) = 0.14 ≈ 1.3·h`).
  3. **Sem assimetria tip-baia (FALHA CONCEITUAL):** mecanismo `dm ∝ |∇c_o|` aplicado por particula do rim NAO discrimina pontas de baias. Todas as particulas com `rho_b ∈ [0.1, 0.6]` recebem influxo de magnitude similar porque `|∇c_o|` na fronteira biomassa-agar e dominado pelo **salto vertical** c_o_int(~1) → c_o_ext(~0), nao pela curvatura azimutal do campo. Para criar a assimetria tip-baia esperada (Srinivasan 2019, Bru 2023), seria necessario que `c_o` **acumulasse no agar das baias** (entre dois dendritos, agar confinado) e **ficasse baixo no agar das pontas** (virgem). Isso exigiria `D_o` muito maior (`L_D_o ~ d_tip-tip ≈ 0.5-1.0`) e mecanismo de confinamento entre dendritos. Implementacao atual nao tem isso.

  **Metricas confirmando falha (log.csv linhas 2-71):**
  - `mass_total` 36.5 → 406.6 (runaway)
  - `dm_influx` 0 → 370 (sem cap)
  - `max_c_o` saturado em 1.0 desde t=2.86s
  - `mean_c_o` plateau 0.012-0.015 (campo nao penetra agar)
  - `contrast_c_o` plateau 65-95 enganoso — alto so porque `mean_c_o` e diminuto
  - `max_v` pico 1.34 em t=18.5s, `n_fast` 280 (alta atividade dirigida pelo influxo, nao pela quimiotaxia)
  - `a_marangoni` 100-150 sustentado (motor cs OK, mas dominado pelo influxo)

  **Pass L permanece BLOQUEADO.** Caminhos para Pass M-A.2 (a discutir):
  - **Fix obrigatorio:** integrar massa com dt-scaling — `d_m += dt * Q0 * gate * grad_mag` ou via `d_am` (analogo a `BiomassGrowth.am`). Sem isso, qualquer ajuste de Q0 e dominado pelo bug de integracao.
  - **Lever 1 — Inverter producao:** osmolitos NAO produzidos pelas bacterias mas **consumidos**, criando acumulo no agar virgem (mais analogo a Pass M-B nutriente consumivel).
  - **Lever 2 — Difusao alta + confinamento bi-escala:** `D_o` 50-100× maior + `λ_o_ext` muito baixo no agar exterior (analogo a Pass J), permitindo `c_o` acumular nas baias entre dendritos onde fluxo difusivo lateral converge.
  - **Lever 3 — Acoplamento via EOS atrativa em vez de massa:** `c_o` modula `B_tension` da EOS. Tips em `c_o` baixo → `B_tension` alta (agar puxa colonia para fora); baias em `c_o` alto → atrocao zerada. Evita o bug de integracao de massa e mantem balanco de massa fechado.
  - **Reconsiderar Pass M-B (Mimura-Murray):** apesar de §3.0 T2/T5 indicarem que swarming P. aeruginosa e nutrient-rich, a logica "consumo cria gradiente que distingue baia de ponta" e geometricamente mais robusta que producao. M-B nao tem o problema de saturacao instantanea e cria assimetria natural via consumo localizado pelas bacterias na trajetoria do rim.

- **Pass M-B.1 — primeiro run com `NutrientConsumption` (2026-05-06, FALHA por anti-difusao — BUG DE SINAL):** `NutrientConsumption` implementada com `D_n=1e-3, k_n=1.0` (corrigido depois para `D_n=0.02, k_n=0.5`). Alem dos parametros incorretos, havia um bug critico de sinal no Laplaciano de Brookshaw: `cn_ij = s_c_n - d_c_n` em vez do correto `cn_ij = d_c_n - s_c_n` (convencao de `SurfactantEquation`). O sinal oposto cria **anti-difusao**: gradientes de `c_n` sao amplificados em vez de suavizados. `max_c_n` explodiu de 1.0 → 10^17 em t=2.4s → `inf` em t=37s. Quando `c_n = inf`, `c_n_factor = inf/(inf+0.1) = nan` em Python → `production = sigma * ... * nan = nan` → `cs = max(1e-9, nan) = 1e-9`. Motor completamente morto em t=40s. Adicionalmente, `c_n_factor` havia sido adicionado incorretamente em `FlagellarForce` (nao planejado), matando tambem a propulsao quando c_n esgotava. Cinco correcoes aplicadas antes do segundo run.

- **Pass M-B.2 — run corrigido (2026-05-06) — AVANCO MORFOLOGICO CRITICO:** Correcoes: (1) `cn_ij = d_c_n - s_c_n` (sinal correto); (2) clamp superior `min(c_n, 1.0)` no integrador; (3) `c_n_factor` removido de `FlagellarForce`; (4) `D_n=0.02, k_n=0.5`; (5) `total_sim_time=100`. Parametros atuais: `sigma=1.5, r_growth=0.05, k_consume=1.0, D_ext=0.08`.

  **Resultado c_n (correto):** `min_c_n` cai de 0.97 → ~0 em t=15s (core depletado), `max_c_n` permanece 1.0 (agar exterior fresco), `mean_c_n` decresce de 1.0 → 0.76 em t=100s, `contrast_c_n` cresce monotonicamente 0.03 → 1.26 (gradiente crescente). Campo funcionando exatamente como esperado — core depletado, exterior fresco.

  **Resultado cs motor:** `a_marangoni` 5-102 durante t=0-100s (**motor vivo todo o run** ✓), `mean_cs` 0.04-0.17 (abaixo do alvo [0.15, 0.45] em parte da simulacao), `contrast_cs` 57-460 (oscilante, frequentemente >100 ✓), `max_cs/mean_cs` ~100-140 (>>4 ✓).

  **Resultado morfologico (frames 050-161, t=40-100s):** **MELHOR MORFOLOGIA JA VISTA NO PROJETO.** Frames 090-161 mostram **8-10 dendritos distintos** com AR ~1:4-1:6, aumentando progressivamente ao longo do tempo. Secondary branching visivel nos frames mais tardios. Nucleo compacto ancorado (hard pinning funcionando). cs concentrado nos tips (halo visivel alem da biomassa nos frames de cs). Flagellar force ativa nas pontas.

  **Comparacao com referencias:**
  - vs Trinschek (b) Fingering: **aproximando da meta** — 8-10 dedos (vs alvo 7-9), AR crescendo, cs halo visivel. Baias parcialmente suprimidas mas nao completamente limpas.
  - vs reference.jpg (PA14): ainda aquem — PA14 tem 15-20 dendritos com baias completamente limpas. Porem a morfologia geral (bracos radiais a partir de nucleo coeso) ja reconhecivel.

  **Problema remanescente — expansao radial do nucleo:** O nucleo (regiao branca rho_b ≥ 0.8) expande ao mesmo ritmo que os bracos crescem. Causas: `BiomassGrowth` cresce `rho_b` E `m` mesmo em particulas ja pinadas (rho_b >= 0.8) — swarmers na borda interna maturaram para biofilme, alargando o nucleo. Como resultado, AR dos bracos plateau em ~1:2.5-3 (nucleo e bracos crescem igualmente). `mass_total` cresce 36.5 → 547.6 (15×) em 100s — excessivo.

  **Diagnostico pos-M-B.2 (2026-05-07) — "Nucleus Maturation" — causa raiz da expansao do nucleo:**

  O gate `rho_b < 0.8` ja estava implementado em `BiomassGrowth` (impede que nucleo ja formado cresça mais). Porem a expansao persistia por **tres mecanismos independentes**, todos decorrentes da ausencia de c_n no `BiomassGrowth`:

  1. **Swarmers das baias maturam para rho_b = 0.8 sem restricao:** quando c_n → 0 nas baias (acontece antes de t=10s, confirmado por `min_c_n → 0`), as bacterias das baias param de produzir surfactante (via `c_n_factor` em `SurfactantEquation`) mas **continuam crescendo rho_b** a taxa completa `r_growth = 0.05`. Um swarmer em rho_b=0.5 atinge 0.8 em ~9s e e pinado — mas estava fisicamente dentro da baia. A fronteira rho_b=0.8 avanca para dentro das baias, preenchendo-as.
  2. **Massa cresce 15x sem controle:** `d_am ∝ rate * m` era aplicado a todos os swarmers (tips E baias). `mass_total` cresceu 36.5 → 547 em 100s — inflacao mecanica via pressao EOS uniforme.
  3. **Interiores dos bracos engrossam:** swarmers atras das pontas (c_n esgotado, mas rho_b < 0.8) continuavam crescendo rho_b → bracos espessavam → AR nao melhorava apesar das pontas avancando.

  **Fix implementado (Pass M-B.3, 2026-05-07):** acoplamento de c_n ao `BiomassGrowth` via fator Michaelis-Menten — ver [equations.py:16-21](src/equations.py#L16-L21):
  ```python
  c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + 0.1)
  rate = self.r_growth * (1.0 - d_rho_b_grown[d_idx] / self.rho_max) * c_n_factor
  ```
  Pontas (c_n ≈ 1.0) crescem normalmente. Baias (c_n ≈ 0) param de crescer biomassa e massa.

  **Resultado M-B.3 (2026-05-11) — FALHA: gating ineficaz por reabastecimento difusivo:**

  Frames e log confirmaram: expansao radial do nucleo e halo entre bracos **persistiram identicos a M-B.2**. `mass_total` cresceu 36.5 → ~550 em 100s (15×, igual a M-B.2). Causa raiz: o Michaelis-Menten funciona apenas se `c_n ~ 0` no local gateado. Com `D_n = 0.02` uniforme (sem bi-escala), o comprimento de penetracao do nutriente e:
  ```
  L_D_n = sqrt(D_n / (k_n · rho_b)) = sqrt(0.02 / 0.3) = 0.258
  ```
  Comparable ao raio da colonia (~1.5). O agar exterior age como reservatorio `c_n = 1.0` que reabastece livremente:
  - **Fronteira do nucleo** (rho_b ≈ 0.6): steady-state `c_n ≈ 0.57` → `c_n_factor = 0.85` → nucleo avanca a 85% da taxa irrestrita.
  - **Baias** (geometricamente abertas ao agar): c_n ≈ 0.3-0.5 → `c_n_factor_baia ≈ 0.80` vs `c_n_factor_tip ≈ 0.90` → seletividade 1.12× — insuficiente para suprimir baias.

  Ambos os problemas tem a **mesma causa raiz**: `NutrientConsumption.loop` usa `D_n = 0.02` uniforme, tratando o biofilme EPS como agar livre para transporte de nutriente.

- **Pass M-B.4 — difusao bi-escala para D_n (2026-05-11, IMPLEMENTADO — aguardando validacao):**

  Analogia direta com Pass I.6 (que resolveu o mesmo problema para cs): a matriz EPS bloqueia o transporte de nutriente dentro do biofilme, assim como bloqueia a difusao de ramnolipideo. Fix implementado em [equations.py](src/equations.py) `NutrientConsumption.loop`:
  ```python
  rho_b_avg = 0.5 * (d_rho_b_grown[d_idx] + s_rho_b_grown[s_idx])
  if rho_b_avg < 0.1:
      D_eff = D_n        # agar livre: 0.02
  elif rho_b_avg < 0.5:
      gate = smoothstep(...)
      D_eff = D_n + (D_n_int - D_n) * gate
  else:
      D_eff = D_n_int    # biofilme EPS: 1e-4
  ```
  Novo parametro: `D_n_int = 1e-4` em [main.py](main.py) e [scheme.py](src/scheme.py).

  **Efeito calculado:**
  ```
  L_D_n_int = sqrt(1e-4 / 0.3) = 0.018  (vs 0.258 antes)
  ```
  Com `L_D_n_int << h` (kernel radius), o nutriente fica confinado ao agar exterior — o biofilme e uma barreira real para o campo `c_n`. A transicao entre nucleo e baia e depletada em `tau ~ h²/D_n_int = 0.072²/1e-4 ≈ 52s` (mas consumo ativo acelera isso).

  **Predicoes M-B.4:**
  - `mass_total` em t=100s: cai para ≈ 50-80 (apenas pontas crescem, baias paradas)
  - Baias: c_n ≈ 0 apos ~5-10s (barreira EPS impede reabastecimento) → suprimidas
  - Nucleo: fronteira permanece fixada no hard pinning (c_n_factor ≈ 0 no interior) → sem expansao radial
  - Motor vivo: tips em agar virgem recebem D_n_ext = 0.02 → c_n ≈ 1.0 → taxa plena
  - AR esperado: ≥ 1:5 (bracos nao engrossam, nucleo nao expande)

  **Resultado M-B.4 (2026-05-11, t=0→79s) — AVANCO MORFOLOGICO REAL + halo residual:**

  Frames 090-110 (t=40-50s) mostram **a melhor morfologia ja vista no projeto**: ~10-12 dendritos longos e finos com bifurcacoes secundarias emergentes (estilo PA14 conforme reference.jpg). `mass_total` cresceu 4.36× em 79s (vs 15× em M-B.3 — melhoria de 3.4×). `mean_c_n` plateau em 0.86 — agar exterior preservado virgem ✅. Motor vivo (`a_marangoni` 16-45, `n_fast` pico 196 em t=50s).

  **Problemas remanescentes diagnosticados (sugestao do usuario validada):**
  1. **Anel interno entre nucleo e dendritos:** baias com `rho_b ∈ [0.2, 0.4]` na smoothstep `[0.1, 0.5]` ainda tem `D_eff ≈ 0.014` (apenas 31% reducao vs agar livre). Reabastecimento difusivo de c_n nas baias persiste.
  2. **Michaelis-Menten suave permite crescimento residual:** `c_n_factor = c_n/(c_n+0.1)` da 50% growth em c_n=0.1, 33% em c_n=0.05. Baias **desaceleram** mas nao param.
  3. **Inflacao mecanica de massa (sintoma):** mass cresceu 4.36× → area 4× → expansao radial mecanica via SPH. Sintoma das causas 1 e 2.

- **Pass M-B.5 — barreira EPS antecipada + cutoff estrito de crescimento (2026-05-11, IMPLEMENTADO — aguardando validacao):**

  Duas correcoes ortogonais aplicadas em conjunto:

  **Fix 1 — antecipar barreira EPS** em [equations.py](src/equations.py) `NutrientConsumption.loop`:
  ```python
  if rho_b_avg < 0.05:   D_eff = D_n          # agar livre
  elif rho_b_avg < 0.25: D_eff = smoothstep(...)
  else:                  D_eff = D_n_int      # ja e biofilme/baia
  ```
  Em `rho_b=0.25` (baia tipica), `D_eff = 1e-4` (eliminacao total do reabastecimento).

  **Fix 2 — cutoff estrito de crescimento** em [equations.py](src/equations.py) `BiomassGrowth.loop`:
  ```python
  if c_n < 0.2:    c_n_factor = 0.0          # baia esgotada → crescimento ZERO
  elif c_n < 0.5:  c_n_factor = smoothstep(...)
  else:            c_n_factor = 1.0          # ponta em agar fresco → taxa plena
  ```
  Substitui Michaelis-Menten suave. Garante que baias param de crescer rho_b E ganhar massa.

  **Predicoes M-B.5:**
  - `mass_total` em t=80s: cai de 159 (M-B.4) para ≈ 60-80 (so pontas crescem)
  - Anel interno entre nucleo e dendritos: desaparece (baias param de maturar rho_b)
  - Dendritos preservam morfologia M-B.4 (frames 090-110) — pontas em c_n>0.5 inalteradas
  - AR esperado: ≥ 1:6 (bracos finos, nucleo fixo, baias vazias)

  **Resultado M-B.5 (2026-05-11, t=0→98s) — FALHA POR PERMISSIVIDADE INVERTIDA:**

  `mass_total` cresceu 14.1× (36.5 → 516 em 98s) — **PIOR que M-B.4 (4.36×)**. Frames 150-229 mostram dendritos visiveis mas anel interno engordou massivamente, colonia preencheu quase todo o dominio.

  **Causa raiz:** smoothstep `[0.2, 0.5]` substituiu Michaelis-Menten `c_n/(c_n+0.1)` na intencao de ser mais estrito, mas a faixa curta saturou em 1.0 para c_n ≥ 0.5 — **mais permissivo** que MM em todo regime c_n ∈ [0.5, 1.0]:

  | c_n | M-B.4 MM | M-B.5 [0.2,0.5] |
  |-----|----------|------------------|
  | 0.9 | 0.90 | **1.00** (+11%) |
  | 0.7 | 0.875 | **1.00** (+14%) |
  | 0.5 | 0.833 | **1.00** (+20%) |

  Como `mean_c_n` plateau em 0.89, a maioria dos swarmers do rim esta em c_n > 0.5 e cresce a **100%**. Inflacao mecanica acelerada.

- **Pass M-B.6 — rampa longa de crescimento [0.2, 0.8] + EPS [0.05, 0.3] (2026-05-11, IMPLEMENTADO):**

  **Fix 1** ([equations.py](src/equations.py) `BiomassGrowth.loop`) — alargar faixa de transicao:
  ```python
  if c_n < 0.2:    c_n_factor = 0.0
  elif c_n > 0.8:  c_n_factor = 1.0
  else:            t = (c_n - 0.2)/0.6; c_n_factor = smoothstep(t)
  ```
  Agora c_n=0.5 cresce a 50% (vs 100% M-B.5, 83% M-B.4) — finalmente mais estrito que MM em todo regime.

  **Fix 2** ([equations.py](src/equations.py) `NutrientConsumption.loop`) — EPS [0.05, 0.3]:
  Transicao ligeiramente mais larga que M-B.5, full block em rho_b=0.3.

  **Predicoes M-B.6:**
  - `mass_total` em t=100s: ≈ 50-80 (faixa de transicao longa freia rim significativamente)
  - Dendritos M-B.4 preservados, anel interno desaparece
  - mean_v previsto: 0.001-0.003 sustentado

  **Resultado M-B.6 (2026-05-11, t=0→99.5s) — MELHOR MORFOLOGIA + tip pumping remanescente:**

  **Frames 090-150 (t=40-60s): a melhor morfologia ja produzida no projeto.** ~15 dendritos longos e finos com bifurcacoes secundarias visiveis, nucleo compacto, baias livres. Padrao reconhecivel vs reference.jpg (PA14).

  Quantitativo: `mass_total` 36.5 → 354 em 99.5s (9.7× — **31% menos** que M-B.5 com 14×). Motor saudavel ate t=72s, depois degrada. `n_fast` 101 (t=48s) → 1 (t=99s). `mean_v` 0.002 → 0.001. `a_pressure` plateau 2.6-3.1 ✓ no orcamento §8 (exceto pico transitorio 20.4 em t=30s — ajuste inicial).

  **Frame 272 (t=99s) — DEGRADACAO:** nucleo expandiu radialmente, dendritos absorvidos pela massa central inflada, morfologia final e blob com extensoes marginais.

  **Causa raiz remanescente — "Tip Pumping":**
  - Pontas em agar virgem: `c_n ≈ 1.0` → `c_n_factor = 1.0` → growth a taxa plena `r_growth = 0.05`
  - `mean_v → 0.001` mas `mass_total` continua crescendo 236 → 354 nos ultimos 16s
  - Massa vem de `BiomassGrowth.d_am`, nao de movimento — pontas em c_n alto ganham massa mesmo sem motor flagellar ativo
  - Pontas eventualmente saturam (rho_b → 0.8) → K.17 pinning ativa → ponta para mas swarmers adjacentes maduram → nucleo radial expande

  **Pinning quimico (proposta intermediaria) NAO ataca tip pumping:** baias ja sao quasi-estaticas em M-B.6 (mean_v = 0.001 sem pinning quimico) porque growth zerou e Marangoni morre sem gradiente cs nas baias. Pinning de baias e redundante. Mass cresce nas pontas, nao nas baias.

- **Pass M-B.7 — reduzir `r_growth: 0.05 → 0.02` (2026-05-11, IMPLEMENTADO — aguardando validacao):**

  Ataque direto ao tip pumping. Alavanca unica em [main.py:55](main.py#L55).

  **Predicoes M-B.7:**
  - Taxa de bombeamento de massa nas pontas: 2.5× menor
  - `mass_total` em t=100s: ≈ 140 (vs 354 em M-B.6) — extrapolacao linear
  - Vida util dos swarmers em rho_b intermediario: estende 2.5× → mais tempo na zona ativa de growth+motilidade
  - Risco: motor de producao cs pode enfraquecer (mais swarmers em rho_b baixo, qs Hill abaixo de saturacao)
  - Morfologia esperada: frames M-B.6 (090-150) preservados ate t=100s

  **Resultado M-B.7 (2026-05-11, t=0→99.8s) — RECORDE HISTORICO DE CONTROLE DE MASSA:**

  `mass_total` cresceu 36.5 → 61.4 em 99.8s (**1.68× — minimo absoluto do projeto**, vs 9.7× em M-B.6, 15× em M-B.3). Motor **acelerando** no final do run: `n_fast` 37 (t=78s) → 189 (t=99.8s), `mean_v` 0.0016 → 0.0036, `a_marangoni` 17-35 sustentado. `mean_cs` 0.267 (alvo [0.15, 0.45] ✅), `contrast_cs` 62.5, `mean_c_n` 0.893 (agar exterior preservado ✅), `min_c_n` = 0 (core depletado ✅), `a_pressure` 2.9-5.0 (orcamento <10 ✅).

  **Bloqueio E (tip pumping) RESOLVIDO ✅:** predicao era `mass_total ≈ 140 em t=100s` — resultado real foi 61.4 (2.3× melhor que predicao). A combinacao `r_growth=0.02 + bi-scale D_n + smoothstep [0.2, 0.8]` conteve o bombeamento de massa nas pontas com eficacia superior ao esperado.

  **Dois problemas residuais identificados pelo usuario (2026-05-11):**

  1. **Halo radial (rim fantasma):** anel de particulas com ρ_b≈0.15 na borda externa da colonia nunca atinge c_n < 0.2 (threshold atual de starvation). Causa: consumo lento (`k_n × 0.15 × c_n ≈ 0.07/s`) + D_n_ext=0.02 reabastece mais rapido que o consumo. O rim continua movendo (Marangoni/flagelar ativo), criando halo isotropico sobre o padrao dendritico.

  2. **Dendritos extrapolando o dominio:** dominio [-3,3]² muito pequeno. Dendritos com v≈0.2-0.3 alcancam a fronteira em t≈40-50s. Paredes solidas sao invisiveis para as equacoes fluidas (`sources=["fluid"]` apenas) — dendritos atravessam sem resistencia, gerando artefatos pos-t=50s.

- **Pass M-B.8 — threshold starvation 0.2→0.4 + chemical pinning c_n<0.4 + total_sim_time=50 (2026-05-11, IMPLEMENTADO — aguardando validacao):**

  **Fix 1 — BiomassGrowth gate** ([src/equations.py](src/equations.py) linha ~25): threshold `c_n < 0.2 → 0.0` alterado para `c_n < 0.4 → 0.0`, faixa de transicao `[0.4, 0.8]` (era `[0.2, 0.8]`). Particulas com ρ_b≈0.15 no rim que tinham c_n≈0.3-0.4 agora tem c_n_factor=0 → crescimento zero.

  **Fix 2 — Pinning quimico** ([src/scheme.py](src/scheme.py) linha 53): `if d_rho_b_grown[d_idx] >= 0.8:` → `if d_rho_b_grown[d_idx] >= 0.8 or d_c_n[d_idx] < 0.4:`. Mecanismo biologico: motor flagelar requer ATP; sem nutriente (c_n < 0.4), flagelo para literalmente. Complementa o pinning mecanico (rho_b >= 0.8). Particulas do halo (c_n esgotado) ficam imobilizadas — continuam existindo no SPH mas nao se movem.

  **Fix 3 — total_sim_time=50** ([main.py](main.py)): captura apenas a fase limpa antes de colisao com fronteira. Evita artefatos de paredes invisiveis. Ativar paredes solidas (Pass L) requer nova fisica de reflexao e contorno — pertence a Pass L, nao agora.

  **Predicoes M-B.8:**
  - Halo radial desaparece — particulas do rim com c_n < 0.4 ficam imobilizadas
  - Motor nas pontas preservado — tips em agar virgem tem c_n ≈ 0.9-1.0 >> 0.4
  - `mass_total` em t=50s: ≈ 45-50 (crescimento confinado as pontas em c_n > 0.4)
  - Morfologia esperada: dendritos M-B.7 sem o anel externo fantasma; AR mais pronunciado
  - Risco: se threshold 0.4 for alto demais, alguns swarmers de borda legitimos (c_n=0.5) em transicao ficam congelados; monitorar n_fast em t=10-20s

  **Resultado M-B.8 (2026-05-11, t=0→50s) — MELHOR MORFOLOGIA DO PROJETO:**

  Frames 015-016 (t≈48-52s) mostram **8-10 bracos dendriticos distintos com ramificacao secundaria emergindo**, nucleo compacto ancorado, cs altamente concentrado em spots nos tips (nao halo uniforme). Esta e a morfologia mais proxima de reference.jpg produzida ate agora.

  **Aviso de diagnostico enganoso:** metricas de log (n_fast=0, mean_cs declining) pareceram indicar "motor morto" mas eram enganosas:
  - `max_v=0.294` a t=35.9s → partículas rapidas existem
  - `a_flag=2.9993` → flagellar force na amplitude maxima, constantemente
  - `a_marangoni=10.5` → Marangoni ainda ativo
  - `mean_cs` recuperou apos t=27s: 0.024 → 0.032 → 0.040 (crescendo)
  - n_fast=0 = sem particulas acima de 0.1, mas 4 tips com v≈0.05-0.3 — CORRETO e DESEJAVEL

  **Mecanismo novo identificado (Tip-Only Swarming):** pinning quimico `c_n<0.4` produziu o efeito biologicamente correto: colonia dormindo + poucos tips ativos. Em M-B.7 havia 100-200 particulas em movimento difuso; em M-B.8 ha 4 tips em movimento concentrado — mais proximo da biologia P. aeruginosa (nucleo estatico, so dendrite tips avancam).

  **Licao #21:** n_fast=0 e mean_v baixo NAO implicam motor morto quando: (a) a_flag e a_marangoni permanecem acima de 5, (b) max_v ainda tem picos 0.1-0.3, (c) frames mostram morfologia saudavel. A interpretacao de "motor colapso" baseada so em n_fast e INCORRETA quando a simulacao esta num regime de tip-only activation. SEMPRE validar morfologia por frames antes de diagnosticar falha.

  **Estado atual (M-B.8 validado):**
  - Bloqueio A (mecanico) ✅ K.17 hard pinning
  - Bloqueio B (quimico) ✅ k_consume=1.0
  - Bloqueio C (geometrico) ✓ parcial — D_ext=0.08, L_D_ext=0.327
  - Frente 6 M-B ✅ rampa longa + EPS [0.05, 0.3]
  - Bloqueio E (tip pumping) ✅ M-B.7 r_growth=0.02
  - Bloqueio F (halo radial) ✓ parcialmente resolvido — halo reduzido vs M-B.7, mas residual visivel
  - **Proximo objetivo: estender total_sim_time=100 para validar sustentabilidade pos-t=50s**

- **Pass M-B.8 validacao t=100s (2026-05-13) — degradacao pos-t=50s diagnosticada:** estendido `total_sim_time=100`. Frame t=50s reconfirmou melhor morfologia do projeto (8-10 dendritos finos com bifurcacao secundaria). Porem frame t=100s exibiu **fragmentacao interna dos bracos**: dendritos com buracos/fraturas no meio, picos anomalos de `a_drag=56` e `a_pressure=8.7` (rogue particles, lição §9 #6), `mass_total` taxa cresceu 2.4× na segunda metade (+6.6% em t=50-100s vs +2.7% em t=0-50s), `mean_cs` subiu monotonicamente 0.038 → 0.061. Tres causas hipoteticadas: (1) pinning quimico `c_n<0.4` cria descontinuidade dentro de bracos longos → brittle neck (lição §15); (2) tip pumping residual; (3) rogue particles. Causa #1 identificada como principal.

- **Pass M-B.9a — soft chemical pin via drag implicito (2026-05-13, FALHA MORFOLOGICA — fragmentacao distribuida):** removido cutoff `c_n < 0.4` do integrador, substituido por starvation drag suave em `LinearDrag.loop`. Tentativa 1 (drag explicito) explodiu numericamente — `a_drag` atingiu 138 000 com `max_v=25.5` em t=0.5s, pq adaptive dt do PySPH escala para ~2e-3 quando `v_max` baixa, violando `dt·γ_max < 2` com `γ_eff` ate 15000. Tentativa 2 (drag IMPLICITO no integrador — `v_new = v_inertial / (1 + dt·γ_eff)`, incondicionalmente estavel) completou os 100s em 307s wall-time.

  **Metricas Tentativa 2:** `contrast_cs=119` (+24% vs M-B.8), `mean_cs=0.047` (-22% vs M-B.8, melhor localizacao), `a_pressure=21` (alto), `a_marangoni=15.9` (+30%), `mass_total=71.18` (igual M-B.8). **Sem spikes anomalos de a_drag** (regime numerico mais limpo).

  **Frame t≈98s (Protocolo §11):** dendritos visiveis mas **piores que M-B.8** — halo grande de **particulas individuais soltas** orbitando a colonia, "trilhas" de particulas detachadas seguindo os bracos. Fragmentacao **distribuida** ao inves de localizada (M-B.8 tinha fraturas internas, M-B.9a tem particulas soltas pelo dominio inteiro).

  **Causa raiz (lição #22):** hard pin rigido mascarava fragilidade da coesao SPH. Substituir `u=v=0` por drag forte (mesmo com γ_eff=5000-15000) permite mobilidade marginal v_term ~ 0.001. Vizinhos com c_n diferente tem v_term diferentes → abre gap → coesao `B_tension=0.00069` nao segura → particula se solta. Pin cinematico era nao-fisico mas estruturalmente protetor. **Reverter para M-B.8 baseline e atacar coesao diretamente (Hipotese 3 obrigatoria).**

- **Pass M-B.9b — revert pin quimico + `tension_ratio: 0.02 → 0.08` (2026-05-13, SUCESSO PARCIAL — melhor frame absoluto do projeto):** mantida toda a config M-B.8; unica alavanca `BiomassEOS(tension_ratio=0.08)` em [src/scheme.py](src/scheme.py). `B_tension` sobe de 0.00069 → 0.00275 (4×).

  **Metricas (t≈97s):**
  - `contrast_cs = 283.3` — **3× MELHOR que M-B.8** (96.4), 2.4× melhor que M-B.9a (119)
  - `max_cs = 13.7` — 2.4× M-B.8 (5.84), cs altamente localizado nas pontas
  - `a_pressure = 2.99` — **dentro do orcamento §8** (vs picos 8.7 M-B.8, 21 M-B.9a)
  - `mass_total = 71.0` — identico a M-B.8
  - `mean_cs = 0.0485` — controle quimico preservado
  - `mean_c_n = 0.944` — agar exterior preservado

  **Frame 010 (t≈50s) — MELHOR FRAME DO PROJETO:** ~15 dendritos finos coerentes com bifurcacao secundaria incipiente, nucleo compacto, cs em pontos discretos nas pontas (red spots nitidos), **sem fragmentacao visivel**. Mais limpo que M-B.8 t=50s.

  **Frame 023 (t≈97s):** ainda apresenta dispersao de particulas dispersas em halo radial — atribuida a **colisao com fronteira do dominio** [-4,4]² (dendritos com `v_term~0.1` alcancam fronteira em t≈40-50s, paredes invisiveis ao fluido `sources=["fluid"]`). Fragmentacao interna dos bracos (M-B.8) e fragmentacao distribuida (M-B.9a) **resolvidas**.

  **Hipotese 3 (coesao reforcada) validada com qualificacoes.** Reforco do `tension_ratio` reduz fraturas internas E melhora dramaticamente localizacao do cs (mecanismo plausivel: braços mais coesos mantem swarmers no rim por mais tempo, sustentando producao localizada). Predicao de risco (re-aparecer halo nas baias) **invalidada** — pinning quimico hard ainda controla baias.

- **Pass M-B.10 — expansao do dominio `[-4,4]² → [-5,5]²` (2026-05-13, SUCESSO — primeira morfologia sustentavel ate t=100s):** unica alavanca em [main.py:34-42](main.py#L34-L42). `x_dim, y_dim: 150→187` para preservar `dx≈0.054`. Wall time 582s (vs 387s M-B.9b) — proporcional ao aumento de particulas (~35k vs 22.5k).

  **Metricas (t≈97s):**
  - `mass_total = 107.47` (vs 71.0 M-B.9b — diferenca de baseline, dominio maior tem mais particulas iniciais)
  - `contrast_cs = 148.4` (vs 283 M-B.9b — reducao proporcional ao espalhamento)
  - `max_cs = 6.31` (vs 13.7 M-B.9b — cs distribuido em mais area)
  - `a_pressure = 2.99` — orcamento §8 ✓
  - `a_marangoni = 11.86` — identico M-B.9b
  - `mean_c_n = 0.965` (vs 0.944 — mais agar virgem disponivel, esperado em dominio maior)
  - `n_fast = 52` — sustentado (vs 31 M-B.9b)

  **Frame 010 (t≈50s):** 8-9 dendritos limpos radiando do nucleo compacto, cs concentrado em red spots nos tips, **sem halo radial de particulas dispersas** — primeira vez no projeto.

  **Frame 023 (t≈97s):** ~10-12 bracos com **AR ~1:6-1:8** (vs ~1:4 truncado em M-B.9b), ramificacao secundaria visivel, baias limpas. **Halo de particulas soltas dramaticamente reduzido vs M-B.9b**. Padrao starfish radial reconhecivel vs reference.jpg.

  **Diagnostico:** confirmado que a fragmentacao pos-t=50s em M-B.8/M-B.9b era **dominada por colisao com fronteira do dominio**, NAO por brittle neck residual. Paredes invisiveis ao fluido (sources=["fluid"]) deixavam dendritos atravessarem a fronteira [-4,4]² em t≈40-50s, gerando halo isotropico de particulas dispersas como artefato numerico. M-B.10 e a primeira morfologia sustentavelmente dendritica ate t=100s do projeto.

  **Estado atual pos-M-B.10:**
  - Bloqueio A (mecanico) ✅ K.17 hard pinning
  - Bloqueio B (quimico) ✅ k_consume=1.0
  - Bloqueio C (geometrico Mullins-Sekerka) ✓ avancado — `D_ext=0.08`, `L_D_ext=0.327`, contraste ~150 sustentado
  - Bloqueio E (tip pumping) ✅ r_growth=0.02
  - Bloqueio F (halo radial baias) ✅ M-B.8 pinning quimico
  - Bloqueio G (brittle neck) ✅ M-B.9b `tension_ratio=0.08`
  - Bloqueio H (colisao fronteira) ✅ M-B.10 dominio [-5,5]²

  **Gaps vs reference.jpg:** dendritos ainda nao tao finos (largura ~3-4 particulas vs 1-2 PA14), tip-splitting incipiente mas nao fractal. **Pass L (rugosidade) agora desbloqueado em principio** — morfologia base esta validada por 100s sem artefatos. Porem considerar primeiro:
  - **Pass M-B.11 (opcional):** `tension_ratio: 0.08 → 0.12` para ganhar marginal contra largura excessiva (lição §23 caminho seguro).
  - **Pass M-B.12 (opcional):** explorar tip-splitting via β/f0 ratio ou ruido estocastico.

- **Pass T1 — Trinschek-like saturation (2026-05-18, SUCESSO TRANSITORIO + FALHA POS-25s):** substitui combinacao `growth_headroom · tip_boost · motile_boost · c_n_factor` + `k_consume` por mecanismo unico tipo Trinschek 2018 ([T1] §3.0): `production = σ · qs · (1 - cs/cs_max) · noise` com `cs_max=0.5`, e decaimento bi-modal `λ_eff = 0.5·λ` no agar (gera halo lento) / `2·λ` no biofilme (mantem contraste). Parametros: `σ: 1.5 → 5.0`, `β: 1.0 → 10.0` (compensa cs_max=0.5 que e ~10× menor que pico antigo).

  **Predicoes ex-ante (analise §3.3.6):** cs deve uniformizar em `cs_max` em todo `rho_b > 0.1` (saturacao local), gradiente apenas radial na fronteira biomassa-agar. Risco identificado: sem `motile_boost`, perda de assimetria azimutal de cs → expansao isotropica em vez de fingering.

  **Resultado (t=0→49s, 13 frames):**

  **Frame 005 (t≈21s) — MELHOR MORFOLOGIA TRANSITORIA JA PRODUZIDA:** 12-15 dendritos finos perfeitos radiando de nucleo coeso, halo cs nitido extendendo-se alem da biomassa. **Visualmente quase identica ao painel (b) Fingering de Trinschek 2018 ([reference_result.png](reference_result.png)).** Primeira vez no projeto que essa morfologia foi atingida.

  **Frames 010-012 (t≈40-45s) — DEGRADACAO CATASTROFICA:** dendritos perdem coesao e viram **linhas de particulas individuais** (single-particle radial streaks), nucleo COLAPSA (centro essencialmente vazio em t=45s), particulas ejetadas radialmente como halo de fragmentacao distribuida (sintoma lição §22).

  **Metricas log.csv (t=0→49.4s):**
  - `max_cs` cravado em `cs_max = 0.4904` desde t=3.8s — saturacao funcionando ✓
  - `mean_cs` 0.002 → 0.076 monotonico (transitorio de saturacao ainda em curso aos 49s)
  - `contrast_cs` colapso monotonico **348 → 6.4** (54×) — sintoma classico §2.4 / lição §8
  - `mean_v` baixissimo (pico 0.0013) — **2.8× menor que M-B.10 (0.0036)**. Colonia praticamente parada apos ejecao radial inicial.
  - `n_fast` pulsando 0-50 — particulas individuais escapando, NAO swarmers organizados nas pontas
  - `a_marangoni` 4.6 → 13.7 (orcamento §8 ✓), `a_pressure` 3.0 estavel (hard pinning OK), `mass_total` 101 → 103 (controlado ✓)

  **Causa raiz da degradacao pos-t=25s — falha de coesao escalada (lição §22 + §23):** com `β=10` + `σ=5.0` + cs saturando em 0.5 no biofilme, o gradiente cs **na fronteira biofilme-agar** e brutalmente afiado (cs cai de 0.5 → 0.075·cs_agar em ~2h kernel). Marangoni resultante na borda e gigantesca (`|∇cs|/h ~ 5`, `a_marangoni_borda ~ β·|∇cs| ~ 50` no kernel). As particulas do rim sao ejetadas radialmente; a coesao EOS atual (`tension_ratio=0.08`, `B_tension≈0.0028`) NAO segura essa tracao amplificada. Cenario direto lição §23: motor amplificado sem reforco proporcional de coesao reproduz brittle neck distribuido.

  **Diagnostico §3.3 (perfil radial cs):** o pico de cs cai **uniformemente em todo o corpo** (`rho_b > 0.5` satura em `cs_max=0.5` por construcao). Nao ha pico LOCALIZADO em mid-arm como em M-B.10 (que tinha `motile_boost` discriminando swarmers ativos). Implicacao: a forca outward atua em **toda a fronteira biomassa-agar simultaneamente** — explica expansao radial uniforme em t=0-20s, sem selecao competitiva. Os 12-15 dendritos perfeitos do frame 005 sao gerados pela **perturbacao inicial cos(8θ)** (modos azimutais) amplificada pelo motor super-criticos antes da degradacao por fragmentacao.

  **Conclusao Pass T1:** mecanismo Trinschek REALMENTE PRODUZ a morfologia certa transitoriamente (revelacao importante — valida teoricamente a abordagem do painel (b)), mas a amplificacao σ=5/β=10 esta acima do que a coesao SPH atual suporta. Tres caminhos para Pass T2 (a decidir):
  - **T2a (mais provavel — caminho seguro):** preservar T1 + reforcar coesao `tension_ratio: 0.08 → 0.20-0.30` para domar tracao na fronteira. Risco: extrair material do agar para dentro da colonia, re-aparecer halo nas baias.
  - **T2b (revisao do motor):** reduzir `β: 10 → 5` ou `σ: 5 → 2.5` para diminuir gradiente cs amplificado. Risco: motor sub-critico, regredir para M-B.10-like.
  - **T2c (hibrido com M-B.10):** restaurar `motile_boost` (localiza cs em swarmers) ATUANDO em conjunto com saturacao `(1-cs/cs_max)`. Combinaria selecao azimutal (M-B.10) com perfil radial saturado (T1). Mais complexo, requer recalibracao de σ/β.

- **Pass T2a — `tension_ratio: 0.08 → 0.20` (2026-05-18, SUCESSO PARCIAL — direcao validada, gap residual):** unica alavanca em [src/scheme.py:134](src/scheme.py#L134). Preserva T1 (σ=5, β=10, cs_max=0.5). Predicao ex-ante (lição #24): razao `|F_mar_borda|/B_tension` cai de ~23000× para ~9000× — ainda acima do limiar (~10×) mas suficiente para reduzir fragmentacao distribuida.

  **Resultado (t=0→48.7s, 13 frames):**

  **Metricas (T2 vs T1, valores em t=49s):**
  - `mass_total` 101.0 → **101.7 (+0.7%)** vs T1 +2.0% — coesao segura material ✅
  - `mean_v` pico 0.0008 vs T1 0.0013 (-38%) — colonia mais "presa" pela coesao maior
  - `n_fast` pico 16 vs T1 50 (-68%) — menos particulas escapando ✅
  - `a_marangoni` pico 10.3 vs T1 13.7 (-25%) — gradiente cs menor (menos drenagem)
  - **`contrast_cs` ESTABILIZOU em 13.0** vs T1 colapso 348→6.4 — **primeira evidencia de plateau** ✅
  - `a_pressure` spikes **6.2 em t=44s** (vs T1 3.0 estavel) — regressao: rogue particles isoladas
  - `max_v` spike **0.52 em t=44s** (vs T1 0.31) — particula isolada com alta v

  **Frames (Protocolo §11):**
  - **T2 frame 005 (t≈20s):** 12-15 dendritos com bordas mais definidas que T1 — nucleo coeso, halo cs nitido.
  - **T2 frame 010 (t≈40s):** dendritos persistem como linhas radiais 1-2 particulas de largura com **nucleo ainda visivel**. T1 nesse momento ja tinha colapso central completo.
  - **T2 frame 012 (t≈48s):** **nucleo AINDA presente** (small green/yellow center), dendritos preservados. T1 frame 012 mostrava centro essencialmente vazio.

  **Diagnostico:** direcao validada (coesao maior preserva nucleo + estabiliza contrast_cs), porem **nao e o ponto de operacao**:
  - ❌ Dendritos ainda finos demais (1-2 particulas vs 3-4 em M-B.10) — coesao reduziu fragmentacao mas nao restaurou bracos solidos.
  - ❌ Rogue particles persistentes (`a_pressure` 6.2, `max_v` 0.52) — escape pontual sem ruptura sistemica.
  - ❌ `mean_v` muito baixo (0.0008) — `tension_ratio=0.20` comeca a se opor ao motor; risco de estagnacao.

  **Caminho de continuidade (uma alavanca por vez):**
  - **T2b:** `tension_ratio: 0.20 → 0.30-0.40` para fechar gap. Risco: motor estagna (mean_v ja cai 38%).
  - **T2c:** se T2b estagnar, reduzir `β=10 → 6-7` + `tension_ratio=0.25` — tradeoff amplitude/coesao.
  - **T2d (audacioso):** restaurar `motile_boost` (lição #3) localizando cs em swarmers, COMBINADO com `(1-cs/cs_max)` — fundir selecao azimutal (M-B.10) com perfil radial Trinschek.

- **Pass T2b — `tension_ratio: 0.20 → 0.30` (2026-05-18, MELHOR PONTO DE OPERACAO DA SEQUENCIA T1):** unica alavanca em [src/scheme.py:134](src/scheme.py#L134). Predicao ex-ante: razao `|F_mar|/B_tension` cai de ~9000× (T2a) → ~6000× — ainda acima do alvo (~10×) mas tendencia consistente.

  **Metricas (T2b vs T2a, valores em t=49.8s):**
  - `mass_total` 101.7 → **101.67** — ≈ identico (+0.7% total, Bloqueio E preservado ✅)
  - `mean_v` pico 0.0008 → 0.00076 (-5%) — coesao maior continua a se opor ao motor, queda marginal
  - `a_marangoni` pico 10.3 → **10.7 (+4%)** — motor preservado ✅
  - **`contrast_cs` (t=49s) 13.0 → 14.3 (+10%) — plateau visivel 14.0-14.7 desde t=37s** ✅ (primeira vez sustentado >15s sob motor T1)
  - **Spikes atenuados drasticamente:**
    - `a_pressure` pico 6.2 → **3.84 (-38%)** ✅
    - `max_v` pico 0.52 → **0.24 (-54%)** ✅
    - `a_drag` pico 40.4 → 29.4 (-27%) ✅

  **Frames (Protocolo §11):**
  - **Frame 005 (t≈20s) — MELHOR FRAME 005 DA SEQUENCIA T1/T2:** 15-18 dendritos com largura **2-3 particulas** (vs 1-2 em T2a), nucleo claramente coeso, halo cs uniformemente distribuido.
  - **Frame 010 (t≈40s):** nucleo CLARAMENTE PRESERVADO (centro green/yellow visivel). Dendritos persistem mas algumas pontas afinam para 1-2 particulas. Sem colapso central.
  - **Frame 013 (t≈50s):** nucleo ainda presente. Algumas pontas mostram **curvaturas terminais** (acumulo de material).

  **Diagnostico:**
  - ✅ Spikes de fragmentacao distribuida atenuados 38-54% — rogue particles em menor quantidade
  - ✅ `contrast_cs` plateau sustentado primeira vez sob motor T1
  - ✅ Massa quasi-constante, Bloqueio E preservado
  - ⚠️ Dendritos voltam a afinar para 1-2 particulas pos-t=30s — coesao insuficiente para manter espessura conforme elongam
  - ⚠️ `mean_v` 0.00076 — abaixo do alvo Pass N (>0.001), queda marginal mas persistente

  **Pass N continua bloqueado** (mean_v < 0.001, dendritos ainda 1-2 particulas em t>30s), porem T2b esta mais perto do baseline que Pass N requer.

  **Caminho de continuidade:**
  - **T2c:** `tension_ratio: 0.30 → 0.40`. Predicao: se mean_v cair <5% adicional, espaco para fechar gap. Se cair >15%, motor comeca estagnar — parar e ir para T2d.
  - **T2d:** reduzir `β: 10 → 7` mantendo `tension_ratio=0.30` — reduz forca amplificada na fronteira proporcionalmente.

- **Pass T2c — `tension_ratio: 0.30 → 0.40` (2026-05-18, PLATEAU DE RETORNOS — REGRESSAO DOS SPIKES):** unica alavanca em [src/scheme.py:134](src/scheme.py#L134). Predicao ex-ante: razao `|F_mar|/B_tension` cai de 6000× (T2b) → 4500× — fechamento adicional do gap; risco de estagnacao do motor (mean_v ja cai).

  **Metricas (T2c vs T2b, valores em t=46s):**
  - `mass_total` 101.66 → **101.63** — ≈ identico (Bloqueio E preservado ✅)
  - **`mean_v` pico 0.00076 → 0.00079 (+4%) — refuta risco de estagnacao** ✅
  - `a_marangoni` pico 10.7 → 10.0 (-7%) — motor preservado mas marginalmente
  - `contrast_cs` (t=46s) 14.0 → 14.0 — **plateau identico**, sem ganho
  - **`a_pressure` spike 3.84 → 5.89 (t=33s) — REGRESSAO +53%** ❌
  - **`max_v` spike 0.24 → 0.42 (t=33s) — REGRESSAO +75%** ❌
  - `a_drag` spike 29.4 → 24.9 (-15%)

  **Frames (Protocolo §11):** dendritos "petala-formato" 2-3 particulas (frame 005), nucleo preservado em frame 012 (t≈48s). **Visualmente quase identicos a T2b** — sem ganho morfologico mensuravel.

  **Diagnostico — PLATEAU DE RETORNOS DA ALAVANCA `tension_ratio`:**
  - ✅ `mean_v` parou de cair (recuperou 4%) — refuta risco estagnacao
  - ❌ Spikes regrediram significativamente (a_pressure +53%, max_v +75%)
  - ⚠️ Frames identicos a T2b — sem ganho morfologico visivel

  **Interpretacao fisica:** coesao extra introduziu **novo modo de instabilidade** — zonas internas com coesao muito alta acumulam energia elastica (analogia mola comprimida) que libera violentamente em eventos esparsos (`a_pressure` 5.89 em t=33s). Fratura discreta ao inves de fragmentacao distribuida. **Função morfologia(tension_ratio) atingiu plateau** — ganho marginal ≈ 0, custo marginal positivo.

  **Decisao: T2b (`tension_ratio=0.30`) e o otimo local desta alavanca isolada.** Reverter para T2b e avancar para T2d (`β: 10 → 7`) para atacar a outra ponta do problema — reduzir tracao amplificada na fronteira em vez de continuar reforcando coesao. Predicao T2d: razao `|F_mar|/B_tension` cai 6000× → ~4200× (proporcional a β); menos energia para gerar rogue particles.

- **Pass T2d — `β: 10 → 5` + `σ: 5 → 20` mantendo `tension_ratio=0.30` (2026-05-18, EQUIVALENTE A T2b — trade-off AR vs mean_v):** mudanca dupla em [main.py:53-54](main.py#L53-L54). β reduzido pela metade (reduz tracao amplificada), σ quadruplicado (acelera saturacao em cs_max). Predicao ex-ante: `|F_mar|/B_tension` cai de ~6000× (T2b) → ~3000× (T2d).

  **Resultado contraintuitivo:** `a_marangoni` caiu apenas 8% (vs predicao -50%). Razao: σ=20 acelera saturacao `cs→cs_max`, aumentando a fracao temporal em que gradiente esta no maximo. Gate `grad_rho_b` ativou em mais particulas (campo cs mais difundido azimutalmente). **Motor operacional preservado apesar de tracao de pico reduzida.**

  **Metricas (T2d vs T2b, valores em t=49.7s):**
  - `mass_total` 101.167 → **101.163** ≈ identico ✅
  - `mean_v` 0.000721 → **0.000638 (-12%)** ⚠️ — colonia mais presa globalmente
  - `a_marangoni` pico 10.75 → 9.93 (-8%) — preservado ✅
  - `contrast_cs` 14.25 → **14.68 (+3%)** ligeiramente melhor ✅
  - `n_fast` 16 → 18 (+13%)
  - `a_pressure` pico 3.84 → **3.56 (-7%)** ≈ similar
  - `max_v` pico 0.24 (t=12s, unico) → **0.31 (t=24s) + 0.21 (t=49s) — dois eventos** ⚠️
  - `a_drag` 29.4 (t=49s) → 21.35 (t=24s) + 19.87 (t=49s) — spike final novo

  **Frames (Protocolo §11):**
  - Frame 005 (t≈20s): 16-18 dendritos similares a T2b mas mais compactos
  - **Frame 010 (t≈40s): dendritos VISIVELMENTE MAIS LONGOS (~r=2.5 vs T2b r=2.0) — extensao radial +25%**
  - Frame 013 (t≈49s): pontas alongadas, nucleo presente, spike final coincide com pontas alcancando ~r=2.8

  **Diagnostico:**
  - ✅ **Primeira melhora de AR pela calibracao T-series** (+25% extensao radial em t=40s)
  - ✅ Spikes intermedios reduzidos
  - ⚠️ mean_v cai 12% — colonia mais "presa" globalmente apesar de pontas alongadas
  - ⚠️ Novo spike terminal (t=49s) — pontas alongadas comecando a fragmentar
  - ⚠️ Largura dendritos nao melhorou (ainda 1-2 particulas)

  **Pass N continua bloqueado:** `mean_v=0.000638 << 0.001` (alvo Pass N). Dendritos ainda 1-2 particulas largura.

  **Conclusao:** T2d e **equivalente a T2b em estabilidade global** com **trade-off bem definido**: dendritos +25% mais longos vs mean_v -12% e novo spike terminal. O σ=20+β=5 mudou regime de "motor forte mas instavel" para "motor difundido com extensao prolongada".

  **Caminhos:**
  - **T2e:** combinar T2d (β=5, σ=20) com `tension_ratio: 0.30 → 0.40` — testa se o plateau de coesao de T2c se devia ao motor T1 amplificado; sob motor T2d reduzido, espaco extra de coesao pode suprimir o spike final.
  - **T2f:** restaurar `motile_boost` mantendo T2d — adiciona localizacao azimutal para tentar engrossar dendritos (1-2 → 2-3 particulas). **Risco documentado** ([feedback_motility_is_tip_discriminator.md]): contradicao com Marangoni physics (inverte gradiente cs) e Xavier/rhlAB CCR (producao deveria ser inversa a c_n).

- **Pass T2d — extensao para validacao em t=100s (2026-05-18, EM ANDAMENTO):** decisao estrategica antes de aplicar Pass N — `total_sim_time: 50 → 100` em [main.py:79](main.py#L79). Justificativa: T2d e o melhor candidato a baseline para Pass N (morfologia extendida +25%, nucleo preservado, massa controlada, contrast_cs plateau), mas o spike terminal em t=49.7s (`max_v=0.21, a_drag=19.87`) sugere que t=100s pode revelar modo de falha tardio. Aplicar Pass N sobre baseline nao-validado em t=100s reproduziria armadilha lição §22 ("Pass N deve coexistir com tension_ratio calibrado, nao substituir") — Pass N estaria mascarando fragmentacao residual em vez de refinar morfologia estavel.

  **Criterios de aceitacao para desbloquear Pass N pos-validacao:**
  - `mean_v` t=50-100s sustentado ≥ 0.0006 (relaxado de 0.001 — gap chicken-and-egg: mean_v baixo PORQUE dendritos finos; Pass N pode resolver)
  - `mass_total` t=100s < 105 (ideal 102-103) — Bloqueio E preservado
  - `contrast_cs` t=50-100s plateau ~14-16 — sem afogamento global
  - `max_v` spikes pos-50s ≤ 0.21 (T2d t=49s spike, nao deve acelerar)
  - `a_pressure` pos-50s < 5 — sem brittle neck novo
  - Frames 14-26: nucleo preservado, dendritos AR ≥ 1:5, sem colisao com fronteira (r < 4 em dominio [-5,5])

  **Decisoes em aberto pos-validacao T2d:**
  - Se PASSAR todos criterios → Pass N desbloqueado; discutir 5 itens pendentes (posicionamento, massa, frequencia, rate limit, interacao com use_splitting)
  - Se FALHAR em mean_v ou massa → T2e (tension_ratio 0.30→0.40) primeiro
  - Se FALHAR em max_v/a_pressure spikes → T2g (reduzir σ de 20 para 10, mantendo β=5) — tradeoff transitorio vs estabilidade
  - Se FALHAR em morfologia (colapso) → reverter para T2b como baseline conservador

- **Pass N v2.4 sobre baseline T2d — FALHA DE PROPOSITO (numericamente OK, morfologicamente no-op) (2026-06-01):** primeiro run com `use_pass_n=True` (v2.4: trigger `σ_a<0.85`, hexagonal n_d==7, `α=ε=0.35`) sobre o motor T2d (`σ=20, β=5, cs_max=0.5, tension_ratio=0.30`), `total_sim_time=50`. **Pre-requisito pulado:** baseline T2d nunca foi validado limpo ate t=100s antes de aplicar Pass N (armadilha lição §22).

  **Sucesso numerico (vs v2.3):** spawns sempre multiplos de 7 (regra n_d==7 funcionando), 1827 spawns total ≈ 261 maes, `mass_total` 101→103 (+2%, de BiomassGrowth nao de Pass N — conservacao OK). `α=ε=0.35` evitou o over-pack sistemico do nucleo da v2.3; run **completou 50s** (vs v2.3 travou em t=13s).

  **Falha de proposito (refinamento e no-op morfologico):** taxa de spawn ~1-2 maes/call e ordens de magnitude lenta demais — dendritos esticam r≈0.5→3 (area ~36×) enquanto Pass N injeta ~5 particulas/~5s. Frame t≈44s (painel rho/rho0): gaps espalhados por **toda** a colonia — refinamento nao sustenta os braços.

  **Degradacao do motor persiste (causa-raiz NAO e resolucao):** narrativa temporal compacto(t4) → **15-18 dendritos belos(t21)** → **fios fragmentados + cs uniformizado(t44)** — identica a lição §24 (T1). `contrast_cs` colapsa 372→9.3, `mean_cs` sobe 0.005→0.0535 (= halo cs enchendo tudo, `contrast_cs ≈ max_cs/mean_cs`). `a_pressure` patologico: pico **123 em t=3s** (over-pack do burst inicial de splits) + spikes 8-33 + subida monotonica terminal 16→21 (instabilidade localizada, NAO motor — `mean_v` parado, `n_fast=0`). `mean_v≈0.001`, regime tip-only/quase-congelado.

  **CUSTO DE dt (descoberta critica):** 580200 iteracoes para 50s (vs 4600 iter para 98s sem Pass N — ver T2g). As filhas com `h=0.35·h_mãe` apertam CFL_force → **dt colapsa ~250×**. Custo de wall-time do Pass N e proibitivo; precisa ser gerenciado (`PASS_N_FREQ` maior, `PASS_N_MAX_PARENTS` menor) antes de qualquer run longo com Pass N ligado.

  **Conclusao:** Pass N v2.4 e sucesso de engenharia numerica (estavel, conservativo, n_d=7 simetrico) mas inutil enquanto o motor T2d/T1 ainda fragmenta sozinho. **Decisao: desligar Pass N (preservar codigo via flag) e validar baseline limpo primeiro → Pass T2g.**

- **Pass T2g — `σ: 20 → 10` + Pass N OFF + `total_sim_time=100` (2026-06-01, SUCESSO — PRIMEIRO BASELINE ESTAVEL ATE t=100s):** alavanca unica no motor ([main.py:55](main.py#L55)), Pass N preservado mas desligado (`use_pass_n=False` [main.py:113](main.py#L113), codigo v2.4 intacto). Objetivo: isolar o motor do baseline (uma alavanca por vez §2.3) e validar o pre-requisito pulado.

  **Descoberta dt:** run chegou a t=98s em **4600 iteracoes** (vs 580200 para 50s com Pass N) — confirma que o colapso de dt ~250× era 100% causado pelas filhas do Pass N (h pequeno → CFL minusculo), nao pelo motor.

  **Metricas (t≈98s):**
  - `a_pressure` = **3.0 plateau** (vs spikes 8-123 com Pass N) — orcamento §8 ✓✓, sem subida terminal (instabilidade T2d eliminada)
  - `a_marangoni` = 8-9 estavel (vs subida terminal a 27) — motor vivo e estavel
  - `contrast_cs` 359 → **12 (quasi-plateau**, t=65→98s: 14→12, deriva ~2 em 33s) — **nao colapsa** (T2d+PassN era 9.3 e caindo)
  - `mean_cs` = 0.040 (subindo lentissimo 0.037→0.040 em 8s) — estavel; `max_cs/mean_cs=11.9 > 4` ✓ (criterio §2.4 #7)
  - `mean_v` = 0.001, `n_fast` 5-11 pulsando — tip-only vivo (lição §21), NAO sub-critico (a_mar=8>>5)
  - `mass_total` 101→106 (+5% em 100s) — controlado, tip pumping leve

  **Frames (§11):** t≈21s ~15-18 dendritos curtos/medios + halo cs; t≈52s dendritos **mais longos e finos** (r≈2.5), baias limpas, **sem fragmentacao**; **t≈98s ~18-20 braços COERENTES persistentes** a r≈2.5-3, ramificacao incipiente, nucleo compacto. **Contraste decisivo vs T2d+PassN:** aquele fragmentava em fios isolados ja em t≈44s; T2g mantem dendritos continuos ate t=98s (>2× o tempo). **Fragmentacao T1/T2d RESOLVIDA.**

  **Validacao vs predicoes:** `a_pressure<6` ✓✓ (3.0); `a_mar` sem spike ✓; fragmentacao reduzida ✓✓ (resolvida); motor nao sub-critico ✓. Parciais: `contrast_cs` plateau ~12 (alvo era ≥16 — abaixo mas **estabilizou** em vez de divergir); `mean_cs` 0.040 (alvo ~0.030 — acima mas estavel). Objetivo central (motor estavel sem fragmentacao ate t=100s) ATINGIDO.

  **Gaps remanescentes vs reference.jpg:** (1) largura dendritos ainda 3-4 particulas (alvo PA14: 1-2) — agora genuinamente problema de **resolucao** = trabalho do Pass N (braços, gate rho_b∈[0.3,0.7]); (2) nucleo esvaziando (gaps no centro, painel 3) — esperado (pin K.17 + Pass N off + gate rho_b<0.7 excluiria nucleo de qualquer forma, decisao v2.3.1).

  **Estado:** T2g e o **baseline validado** para reintroduzir Pass N. Proximo passo: religar Pass N (`use_pass_n=True`) sobre morfologia que NAO fragmenta mais, atacando a largura dos braços — porem com custo de dt (250×) gerenciado via `PASS_N_FREQ`/`PASS_N_MAX_PARENTS` ou runs curtos.

**Invariantes morfologicos descobertos (K.5-K.14):**
1. `smoothstep` em `[a,b]` satura em 1.0 no pico → **max(metric) e cego ao estreitamento do gate**. K.1, K.2 pareceram no-op por isso; diagnostico correto requer `n_active` e `mean_active`, nao `max`.
2. Tip-boost via `|∇rho_b|` e **uniforme no rim** (pontas, baias e trechos retos tem magnitude similar). Nao discrimina pontas sozinho.
3. **Motilidade (|v|) e o discriminador correto**: swarmers ativos estao so nas pontas avancando. Acoplar producao a |v| localiza cs em pontas (K.7 provou — primeira inversao de `contrast_cs`).
4. Feedback `tip_boost × motile_boost` e **multiplicativo** — amplifica exponencialmente. Saturar escala via `sigma` (nao quebrar feedback).
5. `D_ext ∈[0.008, 0.012]` e a janela funcional: menos mata Mullins-Sekerka, mais permite winner-steal.
6. **Multi-seeding sozinho nao seleciona dendritos** (K.11/K.12): com D_ext=0.01 e `lambda_ext_ratio=3.0`, plumes de cs de tips vizinhos se sobrepoe via difusao lateral — 8-10 candidatos nao amadurecem em dendritos distintos. Precisa **confinamento adicional via decaimento**, nao so pre-semeadura.
7. **`lambda_ext_ratio=5.0` confina plumes laterais** (K.13): `L_D_ext 0.149→0.115` impede fusao, produziu dendritos visualmente discerniveis. Porem **nao suficiente**: sem sumidouro fisico, saturacao global do interior destroi gradiente em t>30s.
8. **Ausencia de sumidouro biomassa-dependente e falha terminal** (diagnosticada K.13 → correcao K.14): sem termo `-k·rho_b·cs`, `(1.2-rho_b)` + `motile_boost` geram producao contínua no core saturado. `mean_cs` cresce sem limite, `contrast_cs` colapsa monotonicamente, `mean_v` **cresce aparentemente** por ativacao bulk espuria (n_fast 4x maior com max_v degradando). **Licao critica de diagnostico**: `mean_v` e `n_fast` podem crescer enquanto motor colapsa — sempre validar com tendencia temporal de `contrast_cs` e crescimento de `mean_cs`.
9. **Morfologia `reference.jpg` de P. aeruginosa e FRACTAL, nao radial simples**: ~15 dendritos de primeira ordem + ramificacao secundaria/terciaria por **tip-splitting**. Modelo atual (Marangoni + flagelar + EOS + difusao bi-escala + sumidouro em K.14) gera apenas instabilidade primaria de Mullins-Sekerka → **lobulos bulbosos sem bifurcacao**. Tip-splitting requer mecanismo adicional: substrato consumivel, ruido estocastico forte, ou razao capilar β/f0 otimizada.
10. **Core pinning failure e pre-requisito invisivel para tip-splitting** (diagnosticado K.15, 2026-04-22): ativacao em massa (`n_fast` 900-1300 em ~10⁴ particulas) e `a_pressure` 15-30 (vs orcamento <3) impedem bifurcacao porque pontas nao tem **trajetoria propria** — sao arrastadas pela expansao radial do blob. Na morfologia real (literatura P. aeruginosa PA14), nucleo EPS e imovel e pontas sao a unica frente de avanco. **Pinning tem que vir antes de tip-splitting** — sem nucleo ancorado, nenhuma competicao capilar (β/f0) consegue isolar uma ponta que se bifurque: a proxima iteracao do blob arrasta as duas metades juntas. **Tres alavancas hidrodinamicas**: (a) `gamma_mature` drasticamente maior (5-10× gamma vs 1.5× atual), (b) `c0` reduzido (0.8→0.5 para B/a_pressure drop 2.56×), (c) revisao de `edge_fade` para zerar pressao EOS justamente no nucleo (`rho_b>0.8`) ao inves de maximiza-la.

11. **Bloqueios mecanico e quimico sao FISICAMENTE INDEPENDENTES** (diagnosticado pos-K.19, 2026-04-29): K.16-K.19 demonstrou empiricamente que resolver core pinning mecanico (K.17 hard pinning) **nao e suficiente** para morfologia dendritica. Mesmo com `u=v=0` no nucleo, o motor de producao cs (`σ · qs · (1.2-rho_b) · tip_boost · motile_boost`) na **borda ativa** (rho_b ∈ [0.1, 0.8], que continua se movendo) gera saturacao global. Em K.19 `mean_cs` cresceu 0.13→2.20 com nucleo congelado. **Os dois bloqueios precisam de alavancas diferentes**: (a) mecanico via integrador/drag/EOS, (b) quimico via `k_consume`/coeficiente do `motile_boost`/redistribuicao da producao. Erro original §2.4 (K.15): unir os dois sintomas em "core pinning failure" e bloquear surfactante. Erro corrigido em §2.4 atual.

12. **Reduzir `λ_ext` abaixo de `λ_int` destroi Pass J** (lição K.19): a drenagem que mantem `mean_cs` interior bounded depende do **gradiente de cs na fronteira biofilme-agar**. Quando `λ_ext > λ_int`, cs decai rapido no agar → `cs_ext ≈ 0` → gradiente alto → fluxo difusivo drena cs do interior. Quando `λ_ext = λ_int`, cs acumula no agar → gradiente cai → drenagem cessa → `mean_cs` interno explode. **A drenagem nao e via decaimento direto; e via gradiente que move massa para fora**. Lição metodologica: qualquer mudanca em `λ_ext` deve ser acompanhada de validacao do `cs_ext / cs_int` ratio em log; se subir > 0.3, drenagem comprometida.

13. **Predicoes geometricas (Mullins-Sekerka) requerem bulk com cs bounded primeiro** (lição K.19): a tese "pontas roubam gradiente das baias via L_D_ext lateral" pressupoe que o cs interior e finito e localizado, nao saturado uniformemente. Com `mean_cs >> mean_cs(agar)`, todo o gradiente na borda e dominado pelo salto vertical biofilme-agar, **nao pela curvatura da interface**. Implicacao: bloqueio quimico (B) **deve** ser resolvido antes de qualquer ajuste em `D_ext`/`λ_ext` para focalizacao geometrica (C). Tentar resolver C antes de B (K.18, K.19) so gera predicoes invertidas.

16. **Influxo isotropico no rim NAO seleciona pontas** (lição Pass M-A.1, 2026-05-05): qualquer mecanismo da forma `dm ∝ gate(rho_b) · |∇c_o|` aplicado por particula do rim falha em criar selecao competitiva, porque `|∇c_o|` na fronteira biomassa-agar e dominado pelo **salto vertical** entre interior saturado e exterior virgem — magnitude aproximadamente uniforme ao longo do perimetro, independentemente de a localizacao ser ponta ou baia. O resultado e expansao tipo balao com bumps todos crescendo igualmente. Para criar assimetria tip-baia o campo precisa de **estrutura azimutal no agar** (acumulo de `c_o` nas baias confinadas vs agar virgem nas pontas), o que exige `L_D_o ~ d_tip-tip` e mecanismo de confinamento (decaimento bi-escala ou reflexao em fronteiras). Sem isso, qualquer mecanismo de osmolito producao-acumulo gera borbulhamento ao inves de selecao. **Implicacao para Pass M-A.2:** `D_o` precisa ser comparable a `D_ext` de cs (~0.04-0.08), nao 1e-3.

17. **Integracao de massa em SPH exige dt-scaling explicito** (lição Pass M-A.1, 2026-05-05): `d_m[d_idx] += termo` em `post_loop` SEM multiplicar por dt acumula a `1/dt` por segundo (com dt=5e-5, isso e 20 000×/s). Sintoma: `mass_total` runaway (em M-A.1, 11× em 30s). Padrao correto e (a) `d_m[d_idx] += dt * termo` ou (b) integrar via accumulator `d_am[d_idx]` analogo a `BiomassGrowth.am` que e somado ao `d_m` pelo `CustomEulerStep` com dt-scaling. Generalizacao: **toda equacao SPH que modifica diretamente `d_m`, `d_x`, `d_v` em `post_loop` precisa multiplicar por dt; equacoes que populam acumuladores (`d_au`, `d_am`, `d_a_c_s`) NAO multiplicam por dt** porque o integrador o faz no stage1.

18. **Acoplar c_n ao `BiomassGrowth` e CRITICO para suprimir baias — nao opcional** (lição M-B.3, 2026-05-07): o campo `c_n` acoplado apenas a `SurfactantEquation` e insuficiente para controlar a expansao do nucleo. O mecanismo de falha e o "nucleus maturation": swarmers das baias com c_n esgotado param de produzir surfactante, mas continuam crescendo biomassa (rho_b) a taxa completa. Eventualmente atingem rho_b=0.8 e sao congelados pelo hard pinning — fisicamente dentro da baia. O nucleo avanca para dentro das baias. Simultaneamente, a massa dessas particulas continua crescendo, gerando pressao EOS uniforme (`mass_total` cresce 15× em 100s). **Regra geral:** qualquer campo de recurso limitante (c_n, nutriente, osmolito) cujo esgotamento deva PARAR o crescimento celular PRECISA ser acoplado a `BiomassGrowth`, nao apenas a equacoes de sinalizacao (cs). O sinal (cs) pode ser localizado nas pontas enquanto o substrato (biomassa) continua crescendo uniformemente — os dois mecanismos sao ortogonais.

19. **Gate Michaelis-Menten falha se o campo nao estiver realmente esgotado no local gateado** (lição M-B.3→M-B.4, 2026-05-11): `c_n_factor = c_n/(c_n+K_n)` so suprime crescimento quando `c_n << K_n` (ou seja, `c_n` genuinamente perto de zero). Se a difusao reabastece o campo mais rapido que o consumo, `c_n` nunca cai o suficiente. Com `D_n` uniforme em 0.02, `L_D_n = sqrt(0.02/0.3) = 0.258` — comparavel ao raio da colonia — o agar exterior reabastece nucleus boundary com `c_n ≈ 0.57` e baias com `c_n ≈ 0.3-0.5`. O gate produz apenas 1.12× seletividade tip/baia — insuficiente. **Regra:** antes de implementar qualquer gate por campo escalar, verificar se `L_D = sqrt(D/(consumo*rho_b))` e pequeno comparado ao comprimento caracteristico que se quer criar. Se `L_D` e comparavel ao raio, difusao bi-escala (D_int << D_ext) e necessaria para criar a barreira real. Aplicavel a qualquer campo (c_n, osmolito, nutriente) onde EPS deve funcionar como barreira fisica.

20. **Smoothstep com faixa curta e MAIS permissivo que Michaelis-Menten no regime alto** (lição M-B.5→M-B.6, 2026-05-11): substituir MM `c_n/(c_n+K)` por smoothstep `[a, b]` parece "mais estrito" porque tem cutoff hard em c_n<a, mas **satura em 1.0 para c_n > b** — enquanto MM nunca chega a 1.0 (assintotico). Se a maioria das particulas vive na regiao c_n > b, o smoothstep deixa todas crescerem a 100% vs MM que limitaria a ~b/(b+K). Em M-B.5 (smoothstep [0.2, 0.5] vs MM K=0.1), c_n=0.5 cresceu a 100% vs 83% em MM; com `mean_c_n ≈ 0.89`, mass cresceu 14× em vez de 4×. **Regra:** ao substituir MM por smoothstep com cutoff estrito, verificar o valor de smoothstep no `mean(c_n)` operacional — se for ≈ 1.0, a faixa esta curta demais. Faixa de transicao deve cobrir a maior parte do regime onde as particulas vivem, nao so o regime de cutoff baixo. Em M-B.6 a faixa foi alargada para [0.2, 0.8] cobrindo o regime de operacao.

22. **Pinning rigido pode mascarar coesao SPH insuficiente** (lição M-B.9a, 2026-05-13): substituir `u=v=0` (pin cinematico, nao-fisico mas estruturalmente protetor) por drag forte com γ_eff~5000-15000 (fisicamente correto, v_term ~ 0.001) NAO preserva integridade estrutural se a coesao SPH (`B_tension`) e baixa. Pin rigido transforma um problema mecanico (separacao de vizinhos) em um problema sem fisica (vizinhos nao mexem); drag deixa o problema mecanico ativo, mas com forcas pequenas. Quando vizinhos tem `c_n` diferentes, seus γ_eff sao diferentes, suas v_term sao diferentes, abre gap, e com `B_tension=0.00069` (M-B.8 baseline) nada segura. **Sintoma**: fragmentacao distribuida (particulas individuais soltas em halo) ao inves de fragmentacao localizada (fraturas internas dos bracos). **Regra geral**: ao remover pin rigido, verificar se `B_tension` consegue resistir ao diferencial de v_term tipico no rim. Se `(γ_max - γ_min)·v_typical / B_tension > 1`, e necessario reforcar coesao simultaneamente. Implicacao: **drag suave so e viavel apos coesao SPH ser proporcional aos diferenciais de mobilidade que ele cria**.

23. **Cohesao SPH e tradeoff direto contra fragmentacao interna** (lição M-B.9b, 2026-05-13): quadruplicar `tension_ratio` (0.02→0.08, `B_tension` 0.00069→0.00275) em `BiomassEOS` reduziu fraturas internas dos bracos E melhorou dramaticamente a localizacao do cs (`contrast_cs` 96→283, **3×**). Mecanismo plausivel: bracos mais coesos mantem swarmers do rim em movimento agrupado por mais tempo, sustentando producao localizada de cs nas pontas (motile_boost via |v| coerente) em vez de difundir lateralmente. Confirma que **lição #15 (brittle neck) era o mecanismo dominante de fragmentacao**, e que a coesao 0.02 (legado de I.5 quando o objetivo era "permitir estiramento") era subóptima por margem ampla pos-introducao de hard pinning. Predicao de risco — "atracao extra puxa material das baias e re-introduz halo" — **invalidada** porque pinning quimico hard ainda controla baias. **Regra:** sempre que motor flagelar/Marangoni for amplificado (lição §15) OU pin quimico for usado em conjunto com motor forte, re-avaliar `tension_ratio` no mesmo passo. Caminho seguro de teste pos-M-B.9b: `tension_ratio = 0.12` em M-B.11+ para extrair ganho marginal.

24. **Mecanismo Trinschek-like (saturacao local em cs_max) PRODUZ a morfologia (b) transitoriamente — mas amplificacao σ/β escala alem da coesao SPH** (lição Pass T1, 2026-05-18): substituir o motor multiplicativo (motile_boost · tip_boost · c_n_factor) por `production ∝ (1 - cs/cs_max)` a la Trinschek 2018 produziu em **frame 005 (t≈21s) a melhor morfologia transitoria do projeto** — 12-15 dendritos finos coerentes matching reference_result.png painel (b). Porem amplitude `σ=5, β=10` necessaria para escalar com `cs_max=0.5` gera Marangoni na fronteira (|F_mar_borda| ≈ β·cs_max/h ≈ 70) muito acima do que `tension_ratio=0.08` (`B_tension≈0.003`) suporta — particulas do rim sao ejetadas radialmente em t>25s, dendritos se fragmentam em linhas de particulas isoladas, nucleo colapsa (frames 010-012). **Tres consequencias para futuro:** (a) o mecanismo de saturacao Trinschek **e validado empiricamente** (nao apenas teoricamente) — primeira vez que a morfologia (b) emerge no SPH, mesmo que transitoriamente; (b) sem `motile_boost` ou equivalente, a perda de assimetria azimutal de cs torna a morfologia DEPENDENTE da perturbacao inicial cos(Nθ) — fingerings sao artefato de seeding amplificado, nao de selecao competitiva sustentavel; (c) regra de escala: para qualquer T2 que mantenha `cs_max ~ 0.5`, `tension_ratio` deve escalar como ~`(β/β_ref)·tension_ratio_ref` — para `β=10`, `tension_ratio ≈ 0.20-0.30` (multiplicar tension_ratio_M-B.10 por β/β_M-B.10 = 10x). Lição metodologica: **calcular `|F_mar_borda| / B_tension` antes de mexer em σ ou β** — se razao > 10, brittle neck garantido em t < 30s.

35. **O vacuo tem DUAS regioes fisicamente distintas — estrutural (frozen, preenchivel) e frontier-ativa (motil, NAO-preenchivel por frozen) — separadas em `rho_b≈0.5`** (lição Rota C3.3 + diagnostico do usuario, 2026-06-02): apos C3.3 preencher o vacuo CENTRAL, o usuario notou vacuo remanescente num anel intermediario + nos dendritos. Esses estao na **frontier ativa (`rho_b<0.5`)** — a zona motil (gate flagelar [0.1,0.6], swarmers, tips). C3.3 os EXCLUI (`INSERT_RHO_B_MIN=0.5`) de proposito. **Por que nao da p/ preencher a frontier com a mesma ferramenta:** filler frozen ali CONGELA o motor (#31, shifting v1 provou); biomassa ativa ali realimenta crescimento (#34). O vacuo estrutural (`rho_b>0.5`, junto ao pin congelado) e preenchivel por filler inerte frozen porque NAO e o motor; o vacuo da frontier e o PROPRIO motor em expansao — preenche-lo exige (a) KGC (ja on — corrige o operador SEM preencher, ∇cs valido apesar do vacuo, Liu §3.3); (b) refinamento ATIVO que se move com o arm (Pass N velocidade-herdada + Rota B timesteps individuais [T9], OU inserção ativa C4 a h=h0 com u,v herdados); (c) A3 finer-res global (orcamento de particulas). **Regra:** classifique o vacuo por `rho_b` antes de escolher a ferramenta — frozen-filler so no estrutural (`>0.5`); na frontier ativa, ou se aceita KGC (operador valido) ou se usa refinamento que herda velocidade (nunca frozen). A mesma co-localizacao vacuo↔frontier da lição #32 (shifting) reaparece: na frontier, "mover/congelar particula" sempre conflita com a morfologia; so "corrigir operador" (KGC) ou "refinar herdando velocidade" escapam.

34. **Particula inserida para preencher vacuo TEM que ser filler inerte — se for biomassa ativa, realimenta o vacuo (runaway)** (lição Rota C3.2 t=100s, 2026-06-02): a inserção C3 (lição #33) preencheu o vacuo, mas com as inseridas como BIOMASSA ATIVA (crescem + produzem cs) o run de t=100s teve RUNAWAY: inserções cravaram no cap (102/call de t~67), `mass` +9%, `contrast_cs` colapsou 360→8, estrutura assimetrica. Causa: inseridas na junção (rho_b 0.5-0.6) MATURAM (BiomassGrowth) → expandem o núcleo pinado → **nucleus maturation (#18)** re-disparada → junção maior → mais inserção (feedback +); e produzem cs → inflam o interior → contrast colapsa. **A validacao de t=50s ENGANOU** (auto-limitante ate t~57, runaway depois) — mais um caso §11 de validacao curta insuficiente. **Fix (C3.3):** flag `is_filler`; inseridas NAO crescem (gate em BiomassGrowth), NAO produzem cs (qs=0 em SurfactantEquation), e sao PINADAS (CustomEulerStep) — puro suporte de densidade/kernel, congelado. **Regra geral:** particula adicionada para corrigir um defeito NUMERICO (suporte de kernel) deve ser inerte para a FISICA (nao crescer, nao secretar, nao mover) — senao o remendo numerico vira fonte fisica e realimenta o defeito. Refinamento Vacondio (Pass N) preserva isso dividindo a massa da mae (conserva); inserção C3 ADICIONA massa, entao a inercia fisica tem que ser explicitamente desligada.

33. **Inserir particulas na rede dx (NAO dividir/mover) preenche o vacuo central escapando de TODAS as armadilhas** (lição Rota C3, 2026-06-02): apos Pass N (dt #29), shifting (#31/#32) e KGC (nao cria particula) falharem no vacuo central, a **inserção de particulas frescas na rede dx** nos buracos (rho_b>0.5 & rho/rho0<trig) PREENCHEU o vacuo central pela 1a vez no projeto. Por que escapa de tudo: insere em espaçamento **dx** (densidade local ≈ rho0 → sem over-pack #27); **h=h0** (nao encolhe → dt intacto, escapa #29 que matou o Pass N); **nao move** particula existente (escapa #31/#32 do shifting); inseridas herdam rho_b alto → pinam → congelam estavel. **Distincao-chave vs Pass N:** Pass N DIVIDE (1 mãe→7 filhas coladas a ε·h, over-pack + h pequeno); C3 INSERE na rede (espaçamento dx, h0). O vacuo era deficit de particulas numa colonia em expansao com núcleo pinado — so CRIAR particulas preenche, e criar SEM encolher h nem colar e o que torna barato. **Custo:** `mass_total` cresce (~2× o BiomassGrowth) — inserir = adicionar celulas; controlar via trigger (so vacuo profundo <0.6) e cap. **Regra:** para vacuo que e deficit de contagem (nao de distribuicao nem de operador), inserção-na-rede e a ferramenta; reservar splitting (Pass N) so quando se precisa de RESOLUCAO sub-dx real (feature menor que dx), pareado com timesteps individuais [T9].

32. **O vacuo co-localiza com o frontier motil em `rho_b` — gating por `rho_b` nao separa "preencher" de "proteger"; shifting (Rota A) esgotado** (lição Rota A.3, 2026-06-02, ancorada em Liu §6.4/§6.5): apos a v1 congelar a colonia (#31), a A.3 subiu o gate para `rho_b∈[0.6,0.8)` para deslocar so o interior. Resultado: dendritos voltaram (✅ morfologia, `mean_v`=0.0008≈T2g, dt custo ZERO) MAS o **vacuo da junção núcleo↔rim PERSISTIU** (painel 3: gaps `<0.7` ainda no centro). Causa: as particulas esparsas do anel oco tem `rho_b<0.6` (mesma faixa do frontier que A.3 protege) → excluidas do shift. **`rho_b` nao discrimina vacuo↔frontier — coabitam <0.6.** v1 preenche+congela, A.3 preserva+deixa vacuo: nao ha janela de `rho_b` que faca as duas coisas. Camada profunda: o vacuo e CINEMATICO (hard pin congela o núcleo enquanto o rim expande → junção regenera o gap continuamente); shift de banda fina nao acompanha um vacuo recriado a cada passo. **Regra:** quando o defeito numerico (vacuo/baixo σ_a) co-localiza com a estrutura fisica sensivel (frontier de Mullins-Sekerka), tecnicas que MOVEM particulas (shifting) ou ADICIONAM particulas no frontier (Pass N gate [0.3,0.7]) ficam presas no mesmo dilema. A saida e corrigir o OPERADOR sem mexer na distribuicao — KGC/CSPM [T10] (Rota C): torna ∇cs consistente APESAR do vacuo, sem mover/criar particula, sem depender de separar vacuo↔frontier. Para vacuos CINEMATICOS de pin, "fazer o operador valido apesar do vacuo" supera "tentar preencher o vacuo".

31. **Particle Shifting isotropico no rim MATA a morfologia dendritica — a homogeneizacao que preenche o vacuo e a mesma que aniquila o fingering** (lição Rota A v1, 2026-06-02, ancorada em Lind et al. 2012 [T8] free-surface treatment, Liu §6.5): PST aplicado a TODA a colonia (gate rho_b∈[0.05,0.8)) entregou tudo numericamente — dt custo ZERO (0.0256 ≈ T2g, vs 85× do Pass N), `a_pressure`=3, massa conservada, painel 3 limpo (σ_a→1) — MAS congelou a colonia num disco lumpy r≈1.0 sem dendritos (`mean_v`≈0, vs T2g r≈2.5). Causa: PST e um HOMOGENEIZADOR (criado para bulk uniforme); a morfologia de Mullins-Sekerka e ANTI-homogenea. O shift no frontier (a) amortece as perturbacoes de curto λ que SEMEIAM os dedos, (b) redistribui de volta o push da Marangoni. `contrast_cs`=37 foi falso-positivo (disco compacto de borda nitida, nao pontas — lição §11). **Regra:** qualquer regularizacao geometrica (shifting, XSPH, density diffusion δ-SPH) que toque a INTERFACE EM MOVIMENTO suprime o fingering. Mitigacao canonica (Lind 2012): suprimir a componente do shift NORMAL a superficie livre, OU (simplificado, A.3) gatear o shift SO ao interior estrutural (rho_b>0.6), nunca ao frontier motil (rho_b<0.6). Corolario: o vacuo do INTERIOR (junção núcleo↔rim) pode ser preenchido por shift; o vacuo do RIM ATIVO (braços esticando) NAO — ali so refinamento real (Pass N + timesteps individuais [T9]) ou KGC [T10] (corrige operador sem mover particula).

30. **O vacuo (perda de suporte de kernel) e primariamente problema de DISTRIBUICAO e de OPERADOR, nao de CONTAGEM de particulas — refinamento e so 1 de ≥4 mecanismos** (analise comparativa 2026-06-02, ancorada em Liu §3.3/§6.4/§6.5, Violeau §3.6, [T8]/[T9]/[T10]): o objetivo cientifico e `σ_a≈1` / operadores SPH consistentes na zona ativa (obrigatorio — sem isso o `∇cs` da Marangoni nos braços e espurio e a comparacao com Trinschek invalida). Mas isso NAO implica splitting. A literatura oferece: **(A) Particle Shifting** (Lind/Xu/Adami — redistribui particulas, `δr=-D∇C`, drena σ_a→1 por geometria, **custo dt ZERO**, criada para free-surface/kernel truncado = Liu §6.5); **(B) timesteps individuais** (Springel GADGET-2 + limiter Saitoh-Makino — dissolve a lição #29 deixando h pequeno integrar em passo proprio); **(C) Kernel Gradient Correction** (Bonet-Lok/CSPM — corrige o operador a custo ZERO sem preencher o vacuo). **Regra:** antes de refinar (custo de dt alto), pergunte se o objetivo nao e atingivel por shifting (distribuicao) ou KGC (operador) a custo zero. Splitting so quando a RESOLUCAO FISICA (feature menor que dx) e genuinamente necessaria — e nesse caso pareie com timesteps individuais (B), nao com over-pack. Conflundir objetivo↔mecanismo custou os 85× de dt da v2.4 (lição #29).

29. **O dt adaptativo e fixado pela MENOR partícula, INDEPENDENTE da contagem — refinamento com h pequeno tem custo fixo mesmo com poucos splits** (lição Pass N v2.4 sobre T2g, 2026-06-01, ancorada em Liu §4.2 CFL): o run v2.4 (α=0.35, gen≤1) splitou só **98 mães (686 filhas = +1.7% de partículas)** em 50s, cobertura negligivel — e mesmo assim pagou **~85× de penalidade de dt** (203400 iter vs ~2400 do baseline). Causa: o `dt = CFL·min_a(h_a/...)` e governado pela partícula de menor `h` (filha gen-1, h=0.35·h_mãe=0.034), e 98 filhas pequenas cravam esse mínimo tanto quanto 10000 cravariam. **Consequencia de design:** com `α` pequeno o custo de dt e quase fixo e o benefício escala com a contagem → custo/benefício péssimo para spawn esparso. Duas saídas: (a) `α` grande (h_filha ~ h_mãe → não crava o dt → a contagem deixa de importar, viabilizando enchimento agressivo — Pass N v2.5), aceitando over-pack (lição #27/§10 amendment, Liu §6.5); (b) abandonar refinamento. **Corolário validado:** Pass N entrega RESOLUÇÃO (suporte de kernel, σ_a), NÃO largura física dos dendritos — v2.4 deu `contrast_cs` e morfologia IGUAIS ao baseline T2g (sem Pass N) a 85× o custo. Estreitamento físico dos dedos e problema do MOTOR (comprimento de onda Mullins-Sekerka), ortogonal ao refinamento.

28. **Partição da unidade discreta (σ_a) é o trigger correto para refinamento adaptativo — `ρ_rel < 0.7` é enviesado por massa** (lição Pass N v2.4, 2026-05-28, ancorada em Violeau §3.6 e Liu §3.3.3): o trigger `ρ_rel = ρ_a/ρ_0 < 0.7` (v2.2-v2.3.1) parece medir gap mas mistura duas coisas: a partição da unidade `σ_a = Σ_j V_j W_aj` (consistencia de ordem zero do SPH) E o fator de escala de massa das filhas (`m_d = m_m / n_d`). Como `ρ_a = Σ_j m_j W_aj`, em regime de uniformidade de massa (gen 0) `ρ_a / ρ_0 ≈ σ_a` e os dois criterios coincidem. **Pos-refinamento**, filhas em região esticada têm `σ_a < 1` (gap real, kernel mal suportado) MAS `ρ_a` artificialmente alto se cercadas de irmãs da mesma mãe a curta distância (`ε·h`) — o trigger por densidade silencia o sinal real do erro do operador SPH. **Regra geral**: para refinamento adaptativo, sempre disparar com base em medida invariante sob refinamento — Liu §3.3.3 e Violeau §3.6 identificam `σ_a < 0.85` como criterio de ~15% erro nos operadores SPH, válido independente de geração. Implementacao em [equations.py](src/equations.py) `KernelSum` em Group separado apos `SummationDensity` (necessario para `s_rho` finalizado).

27. **Pinning cinemático + refinamento simétrico são incompatíveis — refinar no núcleo é academic** (lição Pass N v2.3 → v2.3.1, 2026-05-26, ancorada em Liu §4.5 e §6.4): o hard pin K.17 (`u=v=0` em `ρ_b ≥ 0.8`) zera velocidade APÓS o cálculo da força no integrador, quebrando reciprocidade Newton-3: a partícula vizinha não-pinada recebe push integral de pressão mas a reação é zerada artificialmente. Quando v2.3 desbloqueou refinamento no núcleo (gen-1 via tree_mask), os splits se concentraram lá porque `ρ_b ≥ 0.8` + relaxacao Gaussiana inicial cria gaps no centro. Resultado: `a_pressure = 5.015` sustentado (vs orcamento §8 < 3), `max_v` cravado, **dt collapse 12×** (5e-3 → 4.3e-4), simulação travou em t=13s/iter 7200. Causa raiz dupla: (a) com `ε=0.35` e `α=0.6`, o offset 0.378·dx das filhas-vértice cai no pico do kernel `W(r/h ≈ 0.35)` — densidade pós-split ≈ 1.7×densidade pré-split (Violeau §7.4 prova que `ε/α=1` minimiza esse erro; manter razao ≠ 1 é over-pack garantido); (b) Liu §4.5 alerta que a equação de momentum não é integrada em pin cinemático — refinar lá adiciona vizinhos cujas forças vão ao lixo, mas cuja contribuição ao kernel inflaciona a densidade dos vizinhos NÃO-pinados, criando instabilidade de pressão por toda a borda. **Regra geral**: refinamento adaptativo deve ser EXCLUÍDO de regiões com pin cinemático ativo (Liu §6.4 — tensile instability se agrava). Em Pass N v2.3.1 isso virou o gate `ρ_b < PASS_N_RHO_B_MAX = 0.7`. Trade-off explicito: nucleo continua perdendo vizinhos com o tempo sem reposicao — aceito por design porque partição da unidade no núcleo congelado é teoricamente irrelevante.

26. **A mãe sob teste NAO pode estar na KDTree do proximity guard — gen-1 viram "intocaveis"** (lição Pass N v2.2, 2026-05-21): em refinamento hexagonal Vacondio/Feldman, a verificacao de proximidade `d(vk, vizinho) >= prox_min` precisa excluir a propria mae da arvore espacial, porque por construcao geometrica `d(vk, mãe) = r_off = ε·h_mãe`. Se `r_off < prox_min`, a mae aparece como o vizinho mais proximo e **todos os 6 vertices sao rejeitados silenciosamente** (`n_d < 4` → split abortado sem print de falha). Com `ε=0.35`, `α=0.6` e `prox_min=0.4·dx_0`, isso acontece para `h_mãe < 1.143·dx_0` — ou seja, **todas as filhas gen-1 (h=1.08·dx_0)** ficam permanentemente bloqueadas para re-splittar. Sintoma: log mostra splits apenas em iter 200 (gen-0 da relaxacao inicial) e pausa de ~32s antes de retomada esporadica; nucleo esvazia sem ser repovoado; mass_total quasi-invariante (correto por design). **Fix v2.3:** `tree_mask[split_idx] = False` antes de `cKDTree(positions[tree_mask])` ([main.py:451-465](main.py#L451-L465)). **Regra geral:** em qualquer guard de proximidade que verifica vertices geometricamente derivados de uma particula, essa particula precisa ser excluida da arvore de busca, porque a distancia vertice-particula original e DETERMINISTICAMENTE menor que a distancia vertice-qualquer-outra-particula (na regiao de refinamento). Equivale a verificar colisao com um vizinho que ja sabemos que vai ser deletado.

25. **Refinamento "1 filha por gap angular" FALHA estruturalmente; padrão Vacondio/Feldman (7 filhas hexagonais) é obrigatório** (lição Pass N v1, 2026-05-20): Pass N implementado como "detectar gap angular > π/2 e adicionar 1 filha em direção do gap" falhou em 3 modos simultaneamente: (a) **núcleo permanece oco** porque gaps no núcleo são isotrópicos (perda em múltiplas direções por hard pin K.17 + drift do rim), e filtro `max_gap > π/2` os rejeita; (b) **braços permanecem espaçados** porque ganho de resolução por split é 2× (mãe + 1 filha) ao invés de 7× (substituição hexagonal Vacondio) — não acompanha taxa de alongamento dos braços; (c) **velocidade NÃO herdada** (filhas nascem com u=v=0) viola conservação de momentum prescrita por Soleimani 2017 §3.2.6 — gera gradiente artificial e concentra propriedades dinâmicas nas partículas originais. **Solução obrigatória**: implementar Vacondio 2013 / Feldman 2006 (cit. Soleimani §3.2.6): substituir 1 mãe por 7 filhas em padrão hexagonal 2D (1 centro + 6 vértices a 60°), `ε = α = 0.6` (offset e smoothing), massa dividida igualmente (`m_filha = m_mãe/7`), velocidade IDÊNTICA herdada (`u_filha = u_mãe`), todas propriedades (rho_b, cs, c_o, c_n) copiadas para cada filha. Mãe é DELETADA. Erro de densidade < 5% comprovado (Feldman 2006). **Regra geral**: qualquer mecanismo de refinamento SPH neste projeto deve seguir o padrão Vacondio/Feldman — substituição N×, não adição 1×. Antes de propor variação, computar: ganho de resolução = N (não fração) e verificar conservação simultânea de mass + linear momentum + angular momentum. Se algum for violado, voltar ao paper. Diagnóstico em §12 Pass N v1.

36. **Inserção de filler ATIVO na frontier (C4 hibrido) REIGNITE o runaway #34 — inserção de partícula NAO consegue preencher o vácuo do motor sem afogá-lo** (lição C4 híbrido + KGC off, 2026-06-15): testou-se a Opção 2 do roadmap — preencher o vácuo da frontier por inserção (não KGC): KGC OFF, gate estendido a `rho_b>0.1`, filler congelado no anel interno (`rho_b≥0.5` ou `|v|<0.02`) + filler ATIVO (is_filler=0, velocidade herdada) nos dendritos móveis. **Falha em TODOS os eixos vs baseline C3.4 (frozen-só + KGC):** `mass_total` +28.5% (vs +7.5%) — runaway; `contrast_cs` colapso 361→4-6 (vs platô ~14) — motor afogado; `mean_cs` 0.135 (vs 0.035) — inflação química §2.4-B; dt colapsou ~4× (19400 iter p/ t=97 vs 5000) — custo-zero perdido; `frac_lowsig_arms` 0.17-0.30 (NAO melhorou — vácuo não foi preenchido). Frames: blob inchado, cs saturado em todo o domínio, zero dendritos. **Causa raiz (falsifica a hipótese otimista da Opção 2):** o filler ativo no braço (a) PRODUZ cs → infla mean_cs → afoga (#34, exatamente como o Bloqueio Químico); (b) CRESCE biomassa → incha a colônia → cria MAIS vácuo ao expandir → inserções cravam no cap → feedback +; (c) MOVE-se (herda velocidade) → não densifica, alimenta a expansão (por isso frac_lowsig não caiu). **Regra geral:** inserção de partícula só serve para vácuo ESTRUTURAL estático (frozen filler inerte, C3.3). Para o vácuo da FRONTIER ATIVA (= o próprio motor em expansão), NAO existe "filler bom": frozen mata o fingering (#31), ativo reignite o runaway (#34). O vácuo da frontier ou se ACEITA + corrige o operador com KGC (Opção 1 — ∇cs válido apesar do vácuo, e a medição #2 mostrou pontas σ_a=0.98 = sem problema científico), ou exige refinamento Vacondio com conservação de massa + timesteps individuais [T9] (Opção 3, caro). Inserção ad hoc (C4) está DESCARTADA para a frontier. Confirma definitivamente lição #35.

37. **Limite de geração de refinamento por MASSA é derrotado pelo crescimento — usar contador `gen` explícito** (lição Pass N v2.5 Fase 1, 2026-06-15): o Pass N limitava a geração via `m > m_floor=m₀/7` (assumindo que dividir só reduz massa: m_d=m_m/7 ≤ floor ⇒ gen-1 não re-splita). Mas o `BiomassGrowth` CRESCE a massa (`d_m += dt·d_am`) → uma filha gen-1 cresce de volta acima do floor → fica "splittable" de novo → **gen-2,3… espúrios** (confirmado no HDF5: `min m`=m₀/49, `h_min`=0.75²·h0=0.56·h0). Combinado com o over-pack (ε/α=0.47 → `rho/rho0` até 2.94), o h menor colapsou o **dt 97×** e o run parou em t=22.5s (23%) — diagnosticado erroneamente como "afogamento" pelo usuario (mas `mean_cs`=0.023, `contrast_cs`=21, massa conservada → quimica saudavel; era puramente dt). **Regra geral:** qualquer limite de geração/profundidade de refinamento adaptativo DEVE usar um contador explícito que só incrementa (propriedade `gen`, filha=mãe+1, gate `gen<MAX_GEN`), NUNCA um proxy (massa, h, volume) que outra física pode reverter. Massa/h não são monotônicos sob crescimento/EOS. Fix: propriedade `gen` (init 0, herdada +1 nas filhas), removido `m_floor`. Distinção diagnóstica crítica (§11): dt-collapse (sim vira lesma, mas motor vivo) ≠ afogamento químico (mean_cs/contrast colapsam) — sempre cruzar antes de nomear a falha.

38. **Refinar 7× sob dt GLOBAL sempre colapsa o dt — over-pack↔h-pequeno é um bind estrutural, não calibração** (lição Pass N v2.5 Fase 1, 2 runs 2026-06-15): substituir 1 mãe por 7 filhas coloca 7 partículas na área de ~1 → espaçamento ~0.38·dx. A razão de Feldman ε/α (offset/smoothing) governa o trade-off: para NÃO over-packar (ε/α≈1, densidade pós-split correta) seria preciso `h_filha ≈ 0.21·h0` (α≈0.21) → h minúsculo → `dt` colapsa por CFL/viscoso (#29). Para preservar h (α=0.75) → ε/α<1 forçado → filhas dentro do kernel umas das outras → `rho/rho0` até 3.31 → over-pack → `dt` colapsa por pressão. Comprovado: ε=0.35 colapsou em t=19.5s, ε=0.5 só adiou para t≈34s (`rho` ainda 3.31). **Não existe ε/α que escape num integrador de dt único** — é o teorema implícito por trás de [T9] (Springel/Saitoh-Makino): refinamento espacial real EXIGE timesteps locais, senão a partícula mais fina dita o dt de todas (#29). Regra: APR (Vacondio/Feldman) só é viável pareado com timesteps individuais; em PySPH (dt global, `compute_accelerations` monolítico) isso exige forkar o `acceleration_eval`. Alternativa de custo-zero que NÃO refina mas valida o operador: KGC (#2 mostrou pontas σ_a=0.98 → suficiente). Corolário viz: se o vácuo é aceito (Opção 1), o problema VISUAL é de RENDERIZAÇÃO (campo SPH reconstruído/Shepard — Price 2007 splash; ou acumulação temporal "sombra"/footprint), NUNCA de adicionar partículas ao solver (= C4 #34 / freeze #31).

**Proximos diagnosticos obrigatorios (atualizado 2026-04-30 pos-K.21):**

**Estado atual (K.22 aplicado):**
- K.17 hard pinning ✅ (`u=v=0` em `rho_b≥0.8`)
- K.18 `D_ext=0.04` ✅
- K.19 REVERTIDO ✅ (`lambda_ext_ratio=5.0` restaurado)
- K.20 `k_consume=2.0` ✅ (Bloqueio B resolvido)
- K.21 `r_growth=0.15` ✅ (K.17 trap resolvido — swarmer ring sustentado ate t>140s)
- K.22 `f0: 0.5→3.0` ✅ RECEM APLICADO (restaura ignição flagellar)

**Diagnostico K.21 (2026-04-30) — "Motor Starvation / Regime Subcritico":**
Blockers A (pinning) e B (saturacao quimica) confirmados resolvidos pelo log:
- `a_pressure` plateau 0.7-8.6 ✅ (orcamento <10)
- `mean_cs` 0.08-0.18 ✅ (nao explodindo)
- `contrast_cs` 80-240 oscilando ✅ (sinal existe nas pontas)

**Causa raiz K.21:** f0=0.5 produz v_term=0.008 << limiar motile_boost (v_sat=0.1). Sistema em bistabilidade subcritica:
```
v_term(flag) = f0/γ = 0.5/60 = 0.008
motile_boost @ v=0.008: 1 + 50·0.008 = 1.4×  (vs 6× quando v=0.1)
cs_∞(ponta, ρ_b=0.35) ≈ σ·qs·0.85·2.2·1.4 / 1.04 ≈ 2.9  (vs 12 necessario)
```
Resultado: pulsacoes episodicas (n_fast pico=114 em t=1.5s via perturbacao inicial, depois collapse para 0-20). mean_v=0.001-0.002 — insuficiente para elongacao dendritica.

**Pass K.22 — `f0: 0.5 → 3.0` (2026-04-30):** Restaura ignição flagellar ao valor teorico (§8: f0 ~ γ·v_term/2 = 60·0.1/2 = 3). Alavanca unica, parametro flagellar — nao toca cs chemistry nem pinning.

**Predicoes K.22:**
- v_term(flag) = 3/60 = 0.05 → motile_boost = 1+50·0.05 = 3.5×
- cs_∞(ponta, ρ_b=0.35) ≈ 7.3 (vs 2.9 K.21, 12 em K.15)
- cs_∞(interior, ρ_b=1) ≈ 0.10 (k_consume=2.0 inalterado — Blocker B preservado)
- a_flag = 3·gate ≈ 3 (alvo §8 ✓)
- max_cs/mean_cs ≈ 73 >> 4 (criterio #7 ✓)
- mean_v previsto: 0.010-0.030
- n_fast previsto: 50-300 sustentado (vs 0-46 K.21)

1. **Validar K.22:** criterios §2.4-B #5-8 + morfologia dendritica (AR ≥ 1:3 pelo menos). Se mean_cs subir para > 0.45, k_consume ja nao e suficiente — alavanca alternativa: reduzir `motile_boost` factor 50→20.

2. **Pass K.23 (se K.22 estabilizar) — `D_ext: 0.04 → 0.08`:** Mullins-Sekerka geometrico (Blocker C). `L_D_ext = √(0.08/0.75) = 0.327` (+41%). Mantem `λ_ext_ratio=5.0`, drenagem preservada. Criterio: contrast_cs > 150, baias visiveis entre dendritos.

3. **Se K.22 + K.23 produzirem baias mas pontas curtas:** Etapa 2 tip-splitting: razao β/f0 para regime Mullins-Sekerka capilar.

4. **Pass L (rugosidade) permanece BLOQUEADO** ate frames mostrarem morfologia fractal com tip-splitting matching `reference.jpg` de P. aeruginosa.

**Licao #15 — "Brittle Neck" — instabilidade de tração SPH com hard pinning + motor forte (2026-05-01):** Hard pinning (K.17, `u=v=0` em `rho_b≥0.8`) + força flagelar alta (f0=3.0) cria tração de magnitude 3.0 no pescoço entre nucleo congelado e ponta avancando. A coesão EOS (B_tension = rho0·c0²·tension_ratio/7) deve resistir a essa tração; com c0=0.35 e tension_ratio=0.02, B_tension ≈ 0.00035 — irrisório. A viscosidade (mu=0.012, alpha_mon=0.06) também não amortece a separação rápida de partículas. **Resultado:** fluido "quebradiço" — pescoço fragmenta antes de elongar, dendrito não se forma. **Fix:** aumentar mu e alpha_mon (revert parcial I.2) para dar ao fluido propriedades viscoelásticas que aguenham tração. Regra: sempre que f0 for aumentado, verificar se a coesão SPH (mu + alpha_mon + B_tension) escala proporcionalmente.

**Licao #14 — "Motor Starvation" por subcriticidade flagellar (2026-04-30):** f0 foi reduzido de 2.0 (Pass K) para 0.5 (K.5) para "domar" o motor, mas isso empurrou o sistema abaixo do limiar do motile_boost. O motile_boost=1+50·|v| so amplifica producao quando v≥v_sat=0.1; com v_term(flag)=f0/γ << v_sat, o ciclo de feedback nunca arranca. **Diagnostico correto:** comparar v_term(flag) = f0/γ com v_sat. Se v_term < v_sat/3, o sistema esta em regime subcritico e f0 deve ser aumentado antes de qualquer ajuste de cs chemistry.

---

## 10. Protocolo de Trabalho

### Ao receber uma tarefa de calibracao:
1. Ler `log.csv` e identificar a fase atual do motor (bootstrap, ativo, declinante, morto).
2. Consultar o orcamento de aceleracoes (Secao 8) e verificar quais termos estao fora do alvo.
3. Formular predicao quantitativa do efeito da mudanca proposta.
4. Implementar a mudanca (preferencialmente uma variavel por vez).
5. Apos execucao, comparar predicao vs. resultado e atualizar este documento.

### Ao receber uma tarefa de implementacao:
1. Verificar alinhamento com os 3 objetivos de longo prazo (Secao 1).
2. Checar se a mudanca afeta invariantes criticos (ex: ordem das equacoes no `scheme.py`).
3. Implementar com testes minimos de sanidade (ex: verificar que forcas somam zero em equilibrio).
4. Documentar a mudanca neste arquivo na secao apropriada.

### Politica de Comentarios em Codigo

**Regra geral:** Codigo limpo, sem poluicao de comentarios explicativos. Comentarios sao permitidos APENAS em dois casos:

1. **Comentarios de seguranca** (obrigatorios): validacoes de limite ou guarda contra NaN / overflow que nao sao obvias da logica do codigo. Exemplo: `rho_safe = s_rho[s_idx]  # densidade da particula vizinha` + `if rho_safe < 1e-6:  # minimo para evitar NaN`.

2. **Comentarios de contexto historico** (apenas se nao-obvio): Solver customizado, integradores especiais, ou fixes para bugs SPH nao-triviais. **Devem ser movidos para CLAUDE.md**, nao ficar no codigo. O codigo propriamente nao explica o "por que", apenas o "o que". Ver exemplos em §3.0 (referencias da literatura), §2.4 (bloqueios mecanico/quimico), §9 (historico de Passes).

**Proibido:** Comentarios explicativos de algoritmo, detalhe de equacoes, logica de gates, ou historico de decisoes. Tudo isso vai para CLAUDE.md alem do numero de linha (ex: "MarangoniForce.loop" + "gate de interface"). O codigo valida-se por si — bem nomes, estrutura clara, parametros significativos.

**Violacoes encontradas (2026-06-17):** 40+ comentarios explicativos removidos de [src/equations.py](src/equations.py) (docstrings, algoritmo, historico de Passes, estagio de desenvolvimento). Arquivo mantido limpo com APENAS comentarios de seguranca em KernelSum.

### Proibicoes explicitas:
- **Nao** alterar parametros fisicos sem antes ler `log.csv`.
- **Nao** declarar sucesso de um Pass sem evidencia visual da morfologia comparada a `reference.jpg`.
- **Nao** introduzir dependencias externas sem justificativa.
- **Nao** refatorar a ordem das equacoes em `scheme.py` sem verificar invariantes.
- **Nao** remover equacoes desabilitadas (como `OsmoticForce`) que sao parte do roadmap.
- **Nao sugerir implementar rugosidade (Pass L) enquanto a morfologia nao reproduzir `reference.jpg`.** Rugosidade e extensao fisica, nao remedio para motor insuficiente ou selecao competitiva ausente. Se Marangoni + Flagelar + EOS nao geram dendritos finos separados (`AR >= 1:5`), a causa-raiz esta em um desses mecanismos — investigar e refinar antes de adicionar nova fisica.
- **Nao mudar parametros de surfactante (`σ`, `k_consume`, `D_ext`, `λ_ext`) sem verificar criterios duplos §2.4.** Bloqueio mecanico (A) resolvido por K.17. Bloqueio quimico (B) ATIVO — alavanca permitida atual: `k_consume` (lever direta sobre sumidouro). Nao mexer em `D_ext` e `λ_ext` simultaneamente (lição K.18-K.19). Nao reduzir `λ_ext` abaixo de `λ_int` enquanto motile_boost ativo (destroi Pass J).
- **Nao sugerir Pass L (rugosidade)** enquanto morfologia nao reproduzir `reference.jpg` (dendritos AR ≥ 1:5, baias estacionarias, tip-splitting visivel). Bloqueio C (Mullins-Sekerka geometrico) provavelmente requer abordagem apos B resolvido.
- **Nao remover hard pinning mecanico `rho_b>=0.8`** (K.17). Confirmado nao negociavel por K.25a (falha catastrofica) e M-B.9a (fragmentacao distribuida quando substituido por drag forte, lição §22). Pin cinematico e estruturalmente protetor — coesao SPH atual `B_tension≈0.0028` nao resiste ao diferencial de v_term ~ 0.1 entre core e tip sem o pin.
- **Ao introduzir mecanismo que aumenta tracao no rim** (f0 maior, motile_boost maior, etc), **re-avaliar `tension_ratio` no mesmo passo** (lição §23). Cohesao SPH deve escalar com motor — ignorar isso reproduz brittle neck (lição §15).
- **OBRIGATORIO — executar o protocolo §3.3.6 antes de qualquer mudanca em fator de producao ou sumidouro de cs.** Forcas `MarangoniForce` e `FlagellarForce` apontam na direcao `−∇cs` (de alto cs para baixo cs). Para push outward, cs deve **decrescer monotonicamente para fora** ao longo da biomassa — o pico cs deve cair dentro do corpo da colonia, com biomassa contigua entre o pico e o agar. Mudancas em `c_n_factor`, `growth_headroom`, `sigma`, `tip_boost`, `motile_boost`, `qs`, `k_consume`, `lambda` que movam o pico para o agar ou para o rim externo isolado sao PATOLOGICAS — geram push inward na maior parte da colonia. **A falacia "cs concentrado nas pontas puxa para fora" e proibida** — calcular cs_∞ em 4 zonas antes de propor.
- **OBRIGATORIO — qualquer refinamento adaptativo (Pass N e variantes) deve seguir Vacondio 2013 / Feldman 2006 (Soleimani 2017 §3.2.6).** Substituir 1 mãe por N filhas (hexagonal 2D, N=7) com massa dividida IGUALMENTE (`m_filha = m_mãe/N`), velocidade IDÊNTICA herdada (`u_filha = u_mãe, v_filha = v_mãe` — conservação de momentum). **Sobre `α`/`ε` (AMENDADO 2026-06-02, Pass N v2.5):** Feldman 2006 deriva `ε/α=1` como otimo de densidade para particula **isolada**; em dominio **empacotado** isso e inviável (ε grande → filhas colidem com vizinhos a ~dx → proximity guard aborta). Politica do projeto: `α` pode EXCEDER `ε` (over-pack transitorio) quando o objetivo for preservar dt (h_filha grande) E houver coesao+viscosidade suficientes para relaxar (compromisso Liu §6.5). Trade-off explicito — v2.4 usou `α=ε=0.35` (sub-amostragem, dt 85× pior); v2.5 usa `α=0.75, ε=0.35` (over-pack aceito, dt preservado). **Toda escolha `α≠ε` DEVE citar Liu §6.5 e validar o critério "dt avg t>30s ≥ 0.5× inicial"** (teste de relaxacao do over-pack; lição #27 mostra que over-pack não-relaxado trava o run). **Proibido "adicionar 1 filha"** (lição §25 — falha estrutural em 3 modos comprovada por Pass N v1). **Proibido n_d < 7** (licao §27/§28 — splits parciais quebram momento angular e amplificam tensile instability Liu §6.4). Antes de propor variação, computar: (a) ganho de resolução por split = N (não fração), (b) conservação simultânea de mass + linear momentum + angular momentum. Se algum for violado, voltar ao paper. Filhas com `rho_b ≥ 0.8` herdam pin automaticamente do scheme.py.
- **OBRIGATORIO — citar Liu [T6] ou Violeau [T7] em qualquer mudanca de fisica SPH (2026-05-28).** Toda alteracao em (a) refinamento adaptativo, (b) integrador customizado (CustomEulerStep, pinning, drag implicito), (c) equacoes SPH (forcas, difusao, EOS), (d) trigger ou criterio de qualidade do kernel, (e) tratamento de fronteira, **deve referenciar o capitulo/secao especifico de Liu 2003 ou Violeau 2012 que fundamenta a mudanca**. Exemplos validos: "Violeau §3.6 — particao da unidade discreta", "Liu §6.4 — tensile instability sob p<0", "Liu §3.4 — conservacao de momento linear com formulacao antisimetrica". Solucoes "tentativa e erro" sem ancoragem teorica nestas referencias sao **rejeitadas em revisao**. Ver §3.0 [T6]/[T7] para mapa de capitulos relevantes. **A IA deve consultar Liu/Violeau ANTES de propor qualquer modificacao em fisica SPH**, exatamente como ja deve consultar Trinschek/Srinivasan/Bru antes de modificar a fisica biologica.

---

## 11. Protocolo Mandatorio para Analise Morfologica Comparativa

Antes de afirmar que uma morfologia de simulacao esta "de acordo" ou "consistente" com uma imagem de referencia, e imperativo e nao opcional executar a seguinte analise profunda em tres etapas. Uma semelhanca visual superficial e considerada uma **falha analitica**.

### Etapa 1: Analise Morfologica Quantitativa e Qualitativa

Nao se limite a "olhar". Descreva e compare as caracteristicas geometricas fundamentais do padrao:

- **Proporcoes dos bracos:** Analise a razao `comprimento/largura` dos dendritos. Sao finos e alongados (alta razao) ou curtos e grossos (baixa razao)?
- **Morfologia da ponta:** As pontas de crescimento sao afiadas, arredondadas (bulbosas), ou estao se dividindo (bifurcando)?
- **Estrutura de ramificacao:** O padrao exibe ramificacoes secundarias ou terciarias? Ha ausencia total de bifurcacoes (`tip-splitting`)?
- **Densidade e espacamento:** Os bracos sao densamente compactados ou bem espacados?

### Etapa 2: Validacao Fisica Cruzada com Dados (`log.csv`)

Esta e a etapa mais critica. **Toda afirmacao sobre a morfologia deve ser justificada com dados da simulacao.**

- **Valide o crescimento:** Se a morfologia parece estavel ou em crescimento, isso deve ser confirmado por uma `mean_v` (velocidade media) estavel ou crescente no log. Se os bracos parecem estagnados, a `mean_v` deve estar em declinio ou proxima de zero. **Cuidado:** `mean_v` pode crescer espuriamente por ativacao bulk — sempre cruzar com `contrast_cs` tendencia temporal.
- **Valide o motor principal:** A forca motriz deve ser consistente com a forma. Para um crescimento dendritico saudavel via Marangoni, a `contrast_cs` (contraste de surfactante) deve ser **alta e estavel**. Se a morfologia e um "blob" ou os bracos sao bulbosos e estagnados, voce **deve** verificar se a `contrast_cs` esta em colapso, indicando a falha do motor.
- **Valide o equilibrio de forcas:** A forma dos bracos e explicada pelo balanco de forcas? Pontas arredondadas e grossas sugerem que a `a_pressure` (aceleracao de pressao) esta dominando uma `a_marangoni` enfraquecida. Bracos finos requerem um dominio claro da `a_marangoni` na ponta.

### Etapa 3: Sintese e Diagnostico Final

Somente apos cruzar a analise visual (Etapa 1) com a validacao dos dados (Etapa 2), emita um diagnostico.

**Exemplo de diagnostico INCORRETO (superficial):**
> *"Sim, a simulacao parece consistente com a referencia, pois ambos mostram um padrao de crescimento com multiplos bracos."*

**Exemplo de diagnostico CORRETO (rebuscado):**
> *"Nao, a simulacao **falha** em replicar a morfologia de referencia. Enquanto a referencia exibe dendritos finos e ramificados, a simulacao produz bracos curtos e bulbosos que estagnam. A analise do `log.csv` confirma essa falha: a `contrast_cs` entra em colapso ao longo do tempo, enfraquecendo a `a_marangoni`. Isso causa a queda da `mean_v` e permite que a `a_pressure` domine, resultando na morfologia arredondada e estagnada, que e fisicamente e visualmente inconsistente com o alvo."*

### Regra Geral

**Nunca declare conformidade morfologica sem apresentar dados quantitativos do log que justifiquem fisicamente a estabilidade, a dinamica e a forma da estrutura observada.**

Uma violacao deste protocolo ja ocorreu (Pass K.13, 2026-04-20 → 2026-04-22): frame 110 foi declarado como "morfologia indistinguivel de `reference.jpg`" baseado em inspecao visual isolada. Analise temporal subsequente revelou `contrast_cs` em colapso monotonico (131→32) e `mean_cs` em crescimento sem limite (0.26→1.17) — afogamento global do motor. **O Pass foi revertido.** Incidentes deste tipo devem ser impossibilitados pelo cumprimento rigoroso deste protocolo.

---

## 12. Roadmap: Passes I-L (Correcao Morfologica e Extensoes Fisicas)

Diagnostico realizado em 2026-04-16 comparando frames da simulacao (Pass G/H, branch `ram-8827-v1`) com referencia experimental (Michiels et al., reference.jpg). A simulacao produz um **blob ameboide com ~30 bumps curtos**, enquanto a referencia mostra **10-15 dendritos longos e finos** com separacao clara. Cinco causas-raiz identificadas (P1-P5).

### Pass I — Correcao morfologica (dedos longos e separados)

**Objetivo:** Transformar o blob ameboide em uma colonia com 8-12 dendritos finos (razao aspecto >= 1:5), separados por espaco vazio, visualmente similar a referencia.

**I.1 — Unificar e reduzir modos azimutais na condicao inicial**
- **Problema (P1):** N=16 em `rho_b` ([particles.py:32](src/particles.py#L32)) e N=12 em `noise` ([main.py:85](main.py#L85)) interferem, criando ~30 bumps em vez de 10-15 ramos dominantes.
- **Mudanca:** Unificar ambos em N=8. Aumentar amplitude coerente de 0.25→0.4 em `noise`, reduzir ruido de 0.1→0.03. Em `particles.py`: N=16→8, amplitude 0.3→0.5.
- **Predicao:** 8 protuberancias dominantes claras no frame 000 em vez de ~30 ruido-dominadas. Competicao mais limpa → dedos selecionados mais cedo.

**I.2 — Reduzir amortecimento viscoso para permitir filamentos finos**
- **Problema (P2):** `mu=0.025` + `alpha_mon=0.15` suavizam gradientes de velocidade, impedindo formacao de estruturas finas.
- **Mudanca:** `mu: 0.025 → 0.012`, `alpha_mon: 0.15 → 0.06`.
- **Predicao:** Escala de dissipacao viscosa reduzida ~2x. Perturbacoes de comprimento de onda curto crescem em vez de serem amortecidas. Risco: instabilidade numerica se `dt` nao acompanhar (CFL adaptivo deve compensar).

**I.3 — Reduzir difusao para concentrar ∇cs nas pontas**
- **Problema (P3):** `D=4e-3` → `L_D=0.163≈1.5h`. Gradiente de `cs` e suave e largo.
- **Mudanca:** `D: 4e-3 → 1.5e-3`. Novo `L_D = sqrt(1.5e-3/0.15) = 0.1 ≈ 0.93h`.
- **Predicao:** Gradiente ~2.7x mais afiado. `a_mar` pico deve subir de ~50 para ~135 (proporcional a `1/L_D`). Risco: se `L_D < h`, o gradiente nao e resolvido pelo kernel SPH. Com 0.93h estamos no limite — monitorar artefatos de resolucao.

**I.4 — Aumentar contraste de mobilidade core/borda**
- **Problema (P4):** `gamma_mature=0.3*gamma` → razao core/edge = 1.3x (quase nenhuma diferenciacao).
- **Mudanca:** `gamma_mature: 0.3*gamma → 1.5*gamma`. Novo `gamma_eff(rho_b=1) = 60+90 = 150`. Razao core/edge = 150/60 = **2.5x**.
- **Predicao:** Nucleo efetivamente imobilizado. Expansao forcada para as pontas apenas. Combinado com I.2, cria frente de mobilidade afiada.

**I.5 — Reduzir coesao EOS para permitir estiramento**
- **Problema (P5):** `tension_ratio=0.05` puxa material de volta ao corpo, impedindo elongacao.
- **Mudanca:** `tension_ratio: 0.05 → 0.02` em [equations.py:283](src/equations.py#L283).
- **Predicao:** Ramo atrativo 2.5x mais fraco. Dedos podem se esticar sem serem "puxados" de volta. Risco: se muito fraco, particulas na ponta podem se desconectar (fragmentacao).

**Resultado I.1-I.5 (2026-04-16):** Motor +28% mais forte (`a_mar` pico 103 vs 81, `mean_v` 0.027 vs 0.021, `n_fast` 965 vs 769). Morfologia AINDA blob-like — expansao uniforme sem selecao de dedos. Causa-raiz reclassificada: nao e parametrica, e **mecanistica**. O campo de `cs` penetra apenas `L_D=0.1` (~1.7dx) no exterior — insuficiente para focalizacao geometrica Laplaciana nas pontas. Sem campo de longo alcance exterior, todas as secoes do perimetro recebem gradiente identico.

**I.6 — Difusao bi-escala (mecanismo de Mullins-Sekerka)**
- **Problema:** `L_D_int = 0.1` (apenas ~2 particulas de penetracao no exterior). Sem campo de longo alcance, nao ha amplificacao geometrica nas pontas.
- **Mudanca:** `SurfactantEquation.loop` agora usa `D_eff` variavel: `D_int = 1.5e-3` dentro da colonia (`rho_b > 0.5`), `D_ext` no agar (`rho_b < 0.1`), com smoothstep na transicao. Justificativa biologica: ramnolipideo difunde rapido no agar livre e lento no biofilme/EPS.
- **Resultado I.6 (`D_ext=0.03`):** Motor ENFRAQUECIDO 50% (`a_mar` pico 45 vs 103). `D_ext` excessivo (20x `D_int`) drenava `cs` para o exterior mais rapido que a producao. Porem, uma protuberancia elongada formou-se no quadrante inferior — sinal de que a instabilidade geometrica funciona, mas o motor e fraco demais para alimenta-la.
- **Correcao I.7:** `D_ext: 0.03 → 0.01` (7x `D_int` em vez de 20x, `L_D_ext = 0.26 ≈ 2.4h`) + `sigma: 1.2 → 2.0` (+67% producao para compensar drenagem).
- **Resultado I.7 (2026-04-16):** Motor restaurado e mais forte que nunca: `a_mar` pico **148** (vs 103 em I, 45 em I.6), `mean_v` pico 0.025, `n_fast` pico 888. `max_cs = 7.9` (quase 2x Pass I). **MARCO MORFOLOGICO:** frame final mostra **transicao de blob para lobular/dendritico** — tendril dominante no quadrante inferior (selecao competitiva confirmada), perimetro assimetrico com protuberancias de razao aspecto ~1:3-1:4. Primeira evidencia clara de focalizacao geometrica tipo Mullins-Sekerka. Motor ainda morre em t≈36s — piso de producao `(1.2 - rho_b)` insuficiente para sustentar indefinidamente. Gaps restantes: dedos ainda curtos/largos (resolucao), motor nao sustentado (Pass J), ausencia de propulsao ativa (Pass K).

**I.8 — Aumentar resolucao da grade**
- **Mudanca:** `x_dim, y_dim: 100 → 150`. Novo `dx=0.04`, `h=0.072`. Dedo de largura 0.2 tera ~5 particulas (era ~3).
- **Predicao:** Melhor resolucao de estruturas finas. Custo: tempo ~3.4x.

### Pass J — Motor sustentado indefinido

**Objetivo:** Manter motor ativo alem de t=50s. Meta: `mean_v > 0.005` e `n_fast > 100` sustentados ate t=100s.

**Mudanca implementada:** Decaimento espacialmente dependente em `SurfactantEquation.post_loop`:
- Interior (`rho_b > 0.5`): `lambda_eff = lambda = 0.15` — preserva reservatorio de `cs`.
- Exterior (`rho_b < 0.1`): `lambda_eff = lambda_ext = 3 * lambda = 0.45` — remove `cs` rapido no agar.
- Transicao: smoothstep em `rho_b` in `[0.1, 0.5]`.

**Justificativa biologica:** Ramnolipideo no agar livre esta exposto a degradacao ambiental (UV, oxidacao, diluicao por difusao radial). Dentro do biofilme, a matriz EPS protege o surfactante.

**Predicao:** No exterior, `L_D_ext` efetivo cai de `sqrt(D_ext/lambda)=0.26` para `sqrt(D_ext/lambda_ext)=0.15 ≈ 2h`. O gradiente na borda fica mais afiado (`cs` cai mais rapido para zero fora da colonia). O `cs` interior persiste com `tau_int = 1/lambda = 6.7s`, alimentando continuamente a difusao para fora → motor sustentado.

### Pass K — Motilidade flagelar (Frente 5)

**Objetivo:** Adicionar propulsao ativa que sobreviva a saturacao do reservatorio de `cs` e adicione seleção direcional nas pontas.

**Implementacao:** `FlagellarForce(Equation)` em `equations.py:408-483`:
- Padrão `loop`/`post_loop`: `loop` acumula `grad_cs_x/y` pelo gradiente SPH simetrico, `post_loop` normaliza e aplica `-f0 * gate * n̂(∇cs)`.
- Magnitude constante (`f0 * gate`, independente de `|∇cs|`) — escolha-chave para sobreviver a saturacao.
- Gate smoothstep em `rho_b ∈ [0.1, 0.6]` com pico em 0.35.
- `f0 = 2.0`, `r_growth: 0.8 → 0.4` (estende janela de swarmers de borda).
- `grad_cs_x`, `grad_cs_y` registrados como propriedades em `main.py`.
- Posicionada apos `SurfactantEquation` no `scheme` para usar `cs` atualizado.

**Resultado parcial (`log.csv` ate t=24.5s, 2026-04-17):** Motor aparentemente sustentado — `mean_v` estavel 0.0075, `n_fast` 150-300. Numeros de aceleracao reportados como `a_mar` pico 182 sao **invalidos** (bug do `LOG_HEADER`, ver §9). Morfologia: colonia lobular com ~15-20 bumps curtos (AR ~1:2).

**Resultado completo (`log.csv` ate t=99.6s, 2026-04-20):** Motor **NAO** sustentado. Ciclo ativo 0-67s, colapso abrupto 70-72s por saturacao global de `rho_b` zerando o gate flagelar. Frame 390 (t=99s): colonia preenche dominio com ~40 bumps curtos, morfologia ainda nao dendritica. Ver §9 "Em andamento" para analise das 3 fases e causa-raiz.

**Instrumentacao aplicada (2026-04-20):**
- `plot.py` expandido para 1×3 paineis (`rho_b`, `cs`, `|a_flag|`) com `vmax=2.0=f0` no painel de motilidade — permite ver diretamente o gate saturando.
- `au_flag`, `au_mar` adicionados a `add_output_arrays` em [main.py:92](main.py#L92).

**Proximos diagnosticos pendentes (hidrodinamicos — Pass L bloqueado, ver §10):**
1. **Corrigir `LOG_HEADER`** (prerequisito — sem isso, todas as metricas de aceleracao mentem).
2. **Widen flag gate** `rho_b ∈ [0.1, 0.6] → [0.05, 0.8]` para sobreviver a saturacao.
3. Se (2) nao restaurar dendritos: refinar mecanismos hidrodinamicos (estreitar gate Marangoni, aumentar `f0` ou `beta`, revisar ruido inicial). Ver §9 Em andamento para lista detalhada.

### Pass K.1–K.12 — Busca por selecao competitiva de dendritos

Sequencia de intervencoes pos-K iniciada em 2026-04-20 com objetivo de produzir a morfologia de `reference.jpg` (10-15 dendritos AR ≥ 1:5). `LOG_HEADER` corrigido no inicio desta sequencia.

**K.1 (refutado) — estreitar gate Marangoni:** `grad_rho_b_mag ∈ [0.05, 0.6] → [0.15, 0.5]` em `MarangoniForce`. Sem efeito mensuravel. Lecao: `max(a_mar)` e **cego ao estreitamento** porque `smoothstep` satura em 1.0 no pico.

**K.2 (refutado) — estreitar gate flagelar:** `rho_b ∈ [0.1, 0.6] → [0.2, 0.6]` em `FlagellarForce`. `a_flag` permaneceu cravado em 2.0. Mesma razao que K.1.

**K.4 (head-start azimutal) — sem efeito morfologico:** Perturbacao inicial `cos(5θ)` com amplitude 1.0 em rho_b + `sin(5θ)` com amplitude 0.7 em noise, ruido reduzido a 0.01. Motor +65% (`a_mar` pico 306), mas morfologia permaneceu blob porque motor forte **afoga** assimetria inicial.

**K.5 (domar motor) — domou mas nao selecionou:** `beta: 4→1`, `f0: 2→0.5`. `a_mar` pico 63, `mean_v=0.002`, `n_fast=1-23` — motor apagou. Morfologia: blob com ~15 bumps AR ~1:2. Confirmou que "motor forte afoga selecao" era verdade parcial; motor fraco tambem nao seleciona.

**K.6 (tip-boost via `|∇rho_b|`) — amplifica mas nao discrimina:** `production *= (1 + 3*grad_rho_b_mag)`. `max_cs` subiu 9.6 → 27.9 (+190%), `a_mar` 63 → 127. Mas cs forma **anel uniforme**, nao picos nas pontas — `|∇rho_b|` e aproximadamente uniforme em todo o rim (pontas, baias e trechos retos tem magnitude similar).

**K.7 (motility-coupled production) — BREAKTHROUGH + explosao numerica:** `production *= motile_boost` onde `motile_boost = 1 + 50*min(|v|, 0.1)`. **Primeira inversao de `contrast_cs`** (estabilizou em 120-160 apos sempre cair). Frame 030 mostrou **8-10 dendritos radiais** com AR ~1:6. Porem `a_mar` pico **745** (124x orcamento), `max_v=6.1`, dt adaptivo colapsou em t=5.6s. Mecanismo correto, amplitude destrutiva.

**K.8 (dominar via sigma) — one-tip tyranny:** `sigma: 2.0 → 0.6` (3.3x reducao para normalizar producao absoluta). Todas as 7 predicoes quantitativas bateram (`a_mar` 225, `max_v` 1.7, `mean_cs` 0.28, `contrast_cs` 130). **Morfologia: UM dendrito** longo e fino (AR ~1:10) emergindo para cima em frame 034 (t~28s) — qualidade individual identica a `reference.jpg`, mas so um. `n_fast` caiu de 670 → 130 em 10s (winner-takes-all — winner consome gradiente de cs via difusao lateral, extinguindo competidores).

**K.9 (few-tip) — meio caminho:** `sigma: 0.6 → 0.8`. `a_mar` pico 300, `contrast_cs` 95-155, `n_fast` 300-550 sustentado. Frame 039 (t~30s) mostra **3-4 dendritos** com AR 1:3-1:4 + bumps menores. Melhor que K.8 em multiplicidade mas pior em AR individual.

**K.10 (reduzir D_ext) — refutado:** `D_ext: 0.01 → 0.005` tentou localizar cs para prevenir winner-steal. Resultado: **5-6 bumps curtos** AR ~1:2 — **pior** que K.9. Causa: D_ext reduzido matou o campo long-range necessario para focalizacao geometrica Mullins-Sekerka. **Invariante confirmado: `D_ext=0.01` e otimo, nao reduzir.**

**K.11 (multi-seeding N=10 rho_b / N=8 noise) — regressao:** `D_ext` revertido para 0.01, `azimuthal_perturb = 0.8*cos(10θ)` em particles.py. Motor reforcado: `a_mar` pico 368 (+23% vs K.9). Mas `contrast_cs` caiu monotonicamente **188→83** em 16s, `max_cs` declinou apos pico (51.5→33.4). Winner-steal via sobreposicao de plumes: 10 tips semeados interferiram entre si via difusao lateral. Mismatch N=10/N=8 violou licao de coerencia I.1.

**K.12 (unificar N=8 em ambos) — neutra:** `cos(8θ)` tambem em particles.py. `a_mar` pico identico (368), `contrast_cs` caiu igual (204→81). Unificacao de modos **nao resolveu homogeneizacao** — a sobreposicao de plumes independe de modos serem coerentes. **Dado novo**: `n_fast` pulsou 627→165→567 em 12s (oscilacao tipica de ondas sucessivas de ativacao), possivelmente indicando side-branching em t>20s (nao confirmado sem frames).

**Lecoes consolidadas K.1-K.12 (ver §9 Em andamento para lista completa):**
- Metricas `max(·)` sao cegas a gates smoothstep — usar `n_active`/`mean_active`.
- `|∇rho_b|` nao discrimina pontas de resto-do-rim; `|v|` discrimina (motilidade localiza automaticamente).
- Feedback multiplicativo (tip×motile) divergence exponencial — saturar escala via `sigma`.
- `D_ext ∈ [0.008, 0.012]` e janela funcional; fora dela o mecanismo morre.
- **Multi-seeding sozinho nao e suficiente** (K.11/K.12): plumes de cs de tips vizinhos se fundem via difusao lateral no ágar. Precisa **confinamento de plume via `lambda_ext`** (K.13) complementar a pre-semeadura.

### Pass M — Substrato consumível `c_n` (Frente 6) — **PROXIMO PASSO OBRIGATORIO**

> 🎯 **Diagnosticado em K.27 (2026-05-01):** calibracao parametrica esgotada. Os 8 dendritos formam corretamente (K.26-K.27), mas as baias entre dendritos NAO sao suprimidas — particulas no rim recebem push radial uniforme, gerando halo isotropico que sobrepoe ao padrao dendritico. **Mass flow das baias para os tips (esperado biologicamente) nao acontece sem mecanismo ativo de supressao das baias.**
>
> ⚠️ **Atualizacao 2026-05-05 (pos-Pass M-A.1):** primeira tentativa de Pass M-A (osmolito producao + influxo `dm`) FALHOU produzindo borbulhamento isotropico de ~30 bumps uniformes (frames 000-062). Tres falhas independentes diagnosticadas (ver §9 Pass M-A.1): (1) bug de integracao de massa sem dt-scaling causou `mass_total` runaway 11×; (2) `c_o` saturou globalmente em 2.86s; (3) **falha conceitual**: mecanismo `dm ∝ |∇c_o|` por particula do rim nao discrimina pontas de baias (lição #16). Implicacao: a proxima iteracao (Pass M-A.2 ou M-B) deve atacar primeiro a **assimetria geometrica do campo escalar entre baia e ponta**, nao a magnitude do influxo.

**Justificativa biologica:** P. aeruginosa em ágar consome glicose/amino acidos do substrato. Bactérias na baia esgotaram o nutriente local; bactérias nos tips alcançaram ágar virgem com nutriente fresco. Resultado natural:
- Baias param (sem combustivel para metabolismo / producao de surfactante)
- Tips continuam (encontram nutriente em cada avanço)
- Mass flow funnels do bulk para os tips via gradiente azimutal de cs

**Implementacao:**

1. **Novo campo `c_n`** em [particles.py](src/particles.py): `pa.add_property("c_n")` inicializado em 1.0 (agar virgem) para todas as particulas.

2. **Nova equacao `NutrientConsumption(Equation)`** em [equations.py](src/equations.py):
   ```python
   # loop: difusao SPH (similar a SurfactantEquation, sem producao)
   d_a_c_n[d_idx] += 2.0 * D_n * Brookshaw_term
   # post_loop: consumo biomassa-dependente
   d_a_c_n[d_idx] -= k_n * d_rho_b_grown[d_idx] * d_c_n[d_idx]
   ```

3. **Acoplamento em `SurfactantEquation.post_loop`** — modificar termo de producao:
   ```python
   c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + K_n)  # Michaelis-Menten
   production = sigma * qs * (1.2 - rho_b) * noise * tip_boost * motile_boost * c_n_factor
   ```

4. **Acoplamento em `BiomassGrowth.loop`** (✅ implementado em M-B.3 — NAO opcional, CRITICO):
   ```python
   c_n_factor = d_c_n[d_idx] / (d_c_n[d_idx] + 0.1)
   rate = self.r_growth * (1.0 - rho_b / rho_max) * c_n_factor
   ```
   Sem este acoplamento, swarmers das baias maturam para biofilme (rho_b → 0.8) independentemente de c_n, expandindo o nucleo e preenchendo as baias. Ver diagnostico M-B.3 em §9 e licao #18.

5. **Integrar `c_n` em `CustomEulerStep.stage1`**:
   ```python
   d_c_n[d_idx] += dt * d_a_c_n[d_idx]
   d_c_n[d_idx] = max(0.0, d_c_n[d_idx])  # nao pode ser negativo
   ```

**Parametros iniciais propostos:**
- `D_n = 5e-4` (nutriente difunde lento, ~1/3 de D_int para cs)
- `k_n = 1.0` (taxa de consumo: bacteria com rho_b=1 esgota c_n local em ~1s)
- `K_n = 0.1` (Michaelis-Menten: producao cs cai a 50% quando c_n=0.1)

**Predicoes Pass M:**
- Frame 030 (t≈30s): c_n no interior da colonia ja esgotado (~0.05); c_n no rim em transicao; c_n no agar fresco virgem (~1.0)
- Particulas das baias: cercadas por c_n esgotado (atras) e c_n moderado (frente — mas tambem sendo consumido pela frente da baia) → producao cs cai
- Particulas dos tips: avancam em c_n virgem → producao cs sustentada
- **Resultado morfologico:** halo das baias DESAPARECE; mass flow funnels para tips via gradiente azimutal de cs; dendritos finos com agar limpo entre eles.

**Criterio de sucesso Pass M:**
- Frames mostram dendritos finos (largura 1-2 particulas), AR ≥ 1:5
- Agar limpo (sem halo radial) entre dendritos
- mean_v concentrado nos tips (n_fast localizado, nao distribuido em todo rim)
- max_cs/mean_cs > 10 (cs altamente localizado em tips)
- Possivelmente tip-splitting emergente (Mullins-Sekerka secundario)

**Apos Pass M validado:** desbloqueado o Pass L (rugosidade).

---

### Pass N — Refinamento adaptativo de particulas (preenchimento de gaps em tips) — **DIAGNOSTICO DE FALHA + PROPOSTA HEXAGONAL VACONDIO/FELDMAN**

> **Status (2026-05-20):** Pass N v1 (1 filha por gap) IMPLEMENTADO mas FALHOU em t=0→50s. Pass N v2 (hexagonal 7-filhas Vacondio/Feldman) recomendado como substituicao baseado em Soleimani 2017 §3.2.6.

**Objetivo:** preencher espacos vazios entre particulas (gaps estruturais) gerando novas particulas com propriedades herdadas. Em colônia inteira (núcleo + braços + transição), não apenas tips.

**Motivacao biologica:** P. aeruginosa em swarming se divide ativamente para acompanhar a expansão radial. No modelo SPH M-B.10 os dendritos alongam-se de comprimento característico ≈ 0.5 até ≈ 3.0 em ~50s (6× área) sem que partículas novas sejam criadas. Resultado: espaçamento entre partículas cresce de dx≈0.054 inicial para 3-4·dx no rim → kernel sem suporte → forças SPH degradadas. Analogia direta: gap espacial = espaço para divisão celular.

**Motivacao numerica:** `rho/rho0 < 0.7` localmente quebra a hipótese fundamental da EOS quase-incompressível. Preencher gaps restaura o suporte do kernel SPH (alvo: ~35 vizinhos com `h_factor=1.8`). Lições §22-§23 demonstraram que coesão por `tension_ratio` sozinha não resolve gaps emergentes em dendritos esticados.

#### Pass N v1 — IMPLEMENTACAO ATUAL (1 filha por gap angular) — FALHA (2026-05-20)

**Implementacao em [main.py:380-473](main.py#L380-L473):**
- Detecta partículas com `rho_b > 0.05 AND rho_sph < 0.82·rho0`.
- Calcula maior gap angular entre vizinhos no raio 1.5·dx.
- Se `max_gap > π/2`: cria 1 filha em `pos_parent + dx·(cos(gap_mid), sin(gap_mid))`.
- Massa: `m_parent *= 0.5; m_daughter = m_parent` (Option B — conserva massa).
- Propriedades herdadas: `rho_b, cs, c_o, c_n, noise` da mãe; `u=v=0` (NÃO herda velocidade).
- Frequência: 200 iter (~5s); limite 80 por call.

**Resultado (frame 013, t≈48s — `frame_i7_013.png`):**
- ❌ **Nucleo permanece oco** (gap central visível em rho_b e cs panels).
- ❌ **Braços continuam espaçados** (panel rho/rho0 mostra blue circles em gap crítico < 0.7).
- ❌ **`pass_n_spawned` registrou apenas 1-11 partículas por call** (esperado ≫50 para preencher gaps reais).
- ❌ `mass_total` cresce de 101 → 103 em 50s — quasi-conservada (Option B funcional), mas estrutura não recupera.

**Diagnostico das 3 falhas (cruzando com Soleimani 2017 §3.2.6):**

1. **Falha A — "Núcleo oco" (gaps isotrópicos rejeitados):**
   - Hard pin K.17 (`u=v=0` em `rho_b ≥ 0.8`) impede núcleo de fluir para preencher gaps.
   - Marangoni puxa rim para fora → partículas adjacentes ao núcleo migram → perdem vizinhos do núcleo.
   - Núcleo desenvolve gap **isotrópico** (perdeu 3-4 vizinhos em direções variadas), não angular.
   - Filtro `max_gap > π/2` REJEITA esses candidatos porque nenhum gap individual atinge 90°.
   - Mesmo se aceito, 1 filha em 1 direção não compensa déficit de 3-4 vizinhos.

2. **Falha B — "Braços espaçados" (ganho de resolução insuficiente):**
   - Algoritmo Vacondio 2013 (cit. Soleimani §3.2.6): substitui 1 mãe por **N filhas** em padrão hexagonal → ganho 7× de resolução por split (2D).
   - Implementação atual: mãe + 1 filha = **2 partículas onde antes era 1** → ganho 2× apenas.
   - Para um braço que cresceu 6× de área, precisa de 6× mais partículas, não 2×.
   - Próximas calls (a cada 5s) tentam compensar incrementalmente — não acompanham a taxa de alongamento dos braços.

3. **Falha C — "Propriedades concentradas em partículas em movimento" (velocidade NÃO herdada):**
   - Filha inicializada com `u=0, v=0` independente da velocidade da mãe.
   - Soleimani §3.2.6: "daughter particles have the same velocity as their mother particle" — conservação de momentum.
   - Implementação atual VIOLA isso → gradiente de velocidade artificial no nascimento da filha → forças viscosas distorcem estrutura local.
   - Para o usuário: "propriedades concentradas em partículas em movimento" — porque novas partículas nascem PARADAS, e só as ORIGINAIS continuam carregando velocidade + propriedades dinâmicas.

4. **Falha estrutural metodológica (Soleimani §3.2.6 + Feldman 2006):**
   - α e ε não calibrados. Feldman 2006: α=ε=0.6 minimiza erro de densidade a <5% com hexagonal lattice (7 filhas em 2D).
   - Atual: ε efetivo = dx/h ≈ 0.56 (próximo por coincidência), mas `h_daughter = h_mother` (não escala com α=0.6).
   - Consistency degradada porque kernel das filhas tem h grande demais para densidade local pós-split.

#### Pass N v2 — REFINAMENTO HEXAGONAL VACONDIO/FELDMAN (IMPLEMENTADO 2026-05-20)

**Mecanismo correto (Soleimani 2017 §3.2.6, citando Vacondio 2013 + Feldman 2006):**

1. **Substituir 1 mãe por 7 filhas** em padrão hexagonal 2D:
   - 1 filha no centro (posição da mãe original).
   - 6 filhas em vértices a 60°, 120°, 180°, 240°, 300° e 360°.
   - Distância centro-vértice: `ε · h_mae` com **ε = 0.6**.

2. **Massa dividida igualmente:** `m_filha = m_mae / 7` (todas as 7 filhas).

3. **Velocidade herdada IDÊNTICA:** `u_filha = u_mae, v_filha = v_mae` (conservação de momentum linear E angular).

4. **Smoothing length reduzido:** `h_filha = α · h_mae` com **α = 0.6**.

5. **Propriedades herdadas por TODAS as 7 filhas:** `rho_b_grown, cs, c_o, c_n, noise` cópia direta da mãe.

6. **Mãe DELETADA** após gerar as 7 filhas (não permanece — substituição, não adição).

7. **Critério de detecção mais agressivo:** `rho/rho0 < 0.7` (gap real estrutural, não borderline 0.82).

**Riscos remanescentes e mitigações:**

- **Pinning do núcleo:** filhas herdam `rho_b ≥ 0.8` → pin automaticamente ativo via scheme.py. **Sem mudança no scheme**.
- **Multi-geração:** cada split reduz `m` por 7 e `h` por 0.6. Permitir **até 2 gerações** (m_min = m_inicial/49, h_min = 0.36·h_inicial). Marcar como "não-refinável" além disso.
- **Discontinuidade de h no kernel:** kernel das filhas tem h menor → vizinhos antigos (h maior) ainda fazem soma SPH com elas. SPH é robusto a isso (h_ij = 0.5·(h_i + h_j)), mas pode haver transiente. **Mitigação**: chamar `solver.nnps.update()` após split (já feito).
- **Custo computacional:** 7 filhas/split × ~50 splits/call = 350 partículas/call. Para PASS_N_FREQ=200 iter em t=50s (10 calls): ~3500 partículas adicionadas. `n_total` cresce 35000 → 38500 (+10%). Wall time +10-15%.
- **Tip pumping (Bloqueio E):** v2 conserva massa rigorosamente (m_total não muda, apenas redistribui). **mass_total invariante** por construção — risco zerado.

**Parâmetros propostos:**
- `PASS_N_RHO_THRESH = 0.7` (mais estrito que v1=0.82).
- `PASS_N_FREQ = 100 iter` (mais frequente — gaps progridem rapido em braços).
- `PASS_N_MAX_PARENTS = 50` (× 7 filhas = 350 novas por call).
- `PASS_N_ALPHA = 0.6` (smoothing length escala).
- `PASS_N_EPSILON = 0.6` (offset hexagonal).
- `PASS_N_MAX_GENERATIONS = 2` (m_min = m_inicial/49).

**Predicoes Pass N v2:**

| Metrica | M-B.10 (atual sem Pass N) | Pass N v1 (atual) | Predicao Pass N v2 |
|---------|:---:|:---:|:---:|
| `mass_total` em t=50s | 101.7 | 103 (mantém) | **101.7 (invariante)** |
| `n_total` (fluid) | ~35 000 | ~35 050 | **~38 500** |
| Núcleo oco | sim | sim | **não** |
| Braços espaçados (rho/rho0 < 0.7) | sim | sim | **não** |
| Velocidade herdada | n/a | NÃO (viola momentum) | **SIM (idêntica)** |
| `contrast_cs` em t=40-50s | 14-15 | similar | 14-20 (sustentado) |
| Wall time t=0→50s | 387s | ~400s | **~440s (+13%)** |

**Criterio de sucesso Pass N v2:**
- Frames mostram núcleo SEM gap central (rho_b ≈ 1.0 contíguo no centro).
- Frames mostram braços com partículas com espaçamento ~dx (≈ inicial).
- `mass_total` em t=50s invariante (101.7 ± 0.5).
- `contrast_cs` em plateau 14-20 sustentado em t>30s.
- Sem regressão do Bloqueio E (n_fast > 0, mean_v > 0.0005 em t > 30s).
- Sem oscilação numérica (max_v < 1.0 em todo t, sem partículas voando isoladas).

**Implementacao (2026-05-20)** em [main.py:97-104](main.py#L97-L104) (constantes) e [main.py:384-475](main.py#L384-L475) (logica em `SwarmApp.post_step`):
- Substituiu a logica v1 (cKDTree + gap angular + 1 filha) por hexagonal Vacondio/Feldman puro.
- `cKDTree` removido (trigger por volume, sem busca espacial — O(N) ao inves de O(N log N)).
- `x_spawn_ref`/`y_spawn_ref` mantidos no array de partículas mas nao usados (leftover v1).
- Conservação verificada matematicamente:
  - **Massa**: `(N_total - n_mae)·m_normal + n_mae·7·(m_m/7) = (N_total - n_mae)·m_normal + n_mae·m_m` = total invariante.
  - **Momentum linear**: 7 filhas com velocidade idêntica `u_m` → `Σ m_d·u_d = 7·(m_m/7)·u_m = m_m·u_m` = mãe.
  - **Momentum angular**: 6 vértices simétricos em torno do centro da mãe + 1 filha no centro → contribuição rotacional nula = mãe (translação pura).
- `solver.nnps.update()` chamado após append+remove para revalidar lista de vizinhos.

**Nao foi necessario modificar `scheme.py`:**
- `CustomEulerStep.stage1` é genérico — opera por particula sem assumir N constante.
- Hard pin K.17 (`rho_b ≥ 0.8` ou `c_n < 0.6`) propaga automaticamente — filhas herdam rho_b e c_n da mãe, entao herdam o status de pin.
- `SummationDensity` no grupo pre-step recomputa `rho` das filhas a cada step.
- Equações usam `h_ij = 0.5·(h_i + h_j)` automaticamente — variação de h por gerações tratada pelo kernel.

**Validacao pendente — proximos passos:**
1. Rodar `make run` com `total_sim_time = 50` (validacao curta).
2. Verificar nos frames: nucleo SEM gap central, braços com espaçamento ~dx, sem partículas isoladas voando.
3. Verificar no log: `mass_total` invariante (101.7 ± 0.5), `n_total` crescente, `pass_n_spawned` >> v1 (centenas/call).
4. Se validar em t=50s, estender para `total_sim_time = 100s` (corresponde à condição que motivou o trabalho).
5. Se passar em t=100s: documentar como Pass N v2 validado e desbloquear Pass L (rugosidade).

#### Pass N v2.0 — VALIDACAO PARCIAL t=0-13s + 3 FALHAS DIAGNOSTICADAS (2026-05-20)

**Resultado experimental** (50s solicitado, ~25% executado antes de diagnostico):
- `pass_n_spawned`: 287 → 77 → 21 → 14 → 7 → 7 → 14 → **0 ... 0** (zera apos t=7.6s e permanece zero).
- **Anomalia numerica em t=4.73s (call 3):** `max_v = 8.25` (vs ~0.1 normal), `a_pressure = 332` (vs orcamento <10). Pressao explode.
- `mass_total`: 101.08 → 101.43 em ~13s (cresceu ~0.35% — devagar e nao por splits que conservam massa).
- **Frame 016 (t≈12s, painel 3 rho/rho0):** quasi toda colonia em "gap critico" (rho < 0.7·rho0), pior que sem Pass N.

**Diagnostico das 3 falhas:**

1. **Limite gen ≤ 1 muito restritivo** (causa do `pass_n_spawned = 0` perpetuo):
   - `m_floor = m_0 / 7 ≈ 0.000409`. Apos primeiro split, filhas tem `m = m_0/7 = m_floor` exatamente.
   - Condicao `fluid.m > m_floor` e **FALSE** (igualdade nao satisfaz `>`). Filhas gen 1 ficam BLOQUEADAS.
   - Resultado: braços alongando-se no rim pos-t=7.6s nao tem como manter densidade. Splits cessam.

2. **Overlap geometrico ε=0.6 com vizinhos a ~dx** (causa da explosao numerica):
   - `r_offset = ε·h = 0.6 · 1.8·dx = 1.08·dx`. Vertices a 1.08·dx do centro da mae.
   - Particulas vizinhas (gen 0) ja a ~dx do centro da mae.
   - **Distancia filha-vizinho = 0.08·dx — MUITO PERTO**. SPH pressure forces explodem.
   - Feldman 2006 calibrou ε=0.6 para particulas **isoladas**; em dominio empacotado isso colide.

3. **Trigger inicial falso-positivo na borda da Gaussiana** (causa do pico 287 splits no t=2.77s):
   - Gate `rho_b > 0.05` inclui borda dilute da Gaussiana inicial (`R = 0.3 + 0.06·cos(8θ)`).
   - Particulas com `rho_b ∈ [0.05, 0.3]` tem rho SPH naturalmente baixo por kernel truncado (sem vizinhos do lado do agar).
   - Trigger `V > V_limit` dispara nelas mesmo sem necessidade real de refinamento.

#### Pass N v2.1 — 4 CORRECOES APLICADAS (2026-05-20)

**Mudanca de parametros** em [main.py:97-114](main.py#L97-L114):

| Parametro | v2.0 | v2.1 | Razao |
|---|---|---|---|
| `PASS_N_M_FLOOR_RATIO` | `1/7` | `1/49` | Permite gen ≤ 2 (m floor = m₀/49). Particulas no rim alongando podem splitar 2x. |
| `PASS_N_EPSILON` | `0.6` | `0.35` | Vertices a 0.63·dx (vs 1.08·dx). Dentro da "exclusion zone" da mae. Sem overlap com vizinhos. |
| `PASS_N_RHO_B_MIN` | (não havia, gate=0.05) | `0.3` | Gate `rho_b > 0.3` exclui borda dilute da Gaussiana. Cobre swarmers ativos + nucleo. |
| `PASS_N_PROXIMITY_MIN` | (não havia) | `0.4·dx` | Filhas a < 0.4·dx de vizinho existente sao descartadas. Massa redistribuida em N_actual filhas. |

**Mudanca de logica** em [main.py:384-499](main.py#L384-L499):
- Reintroduzido `cKDTree` para busca espacial O(log N) de vizinhos.
- **Proximity guard**: cada um dos 6 vertices testa distancia a (a) particulas existentes via tree, (b) filhas ja adicionadas nesta call via lista `new_positions`. Vertices que falham sao descartados.
- **Mass conservation per mother**: `m_d = m_m / n_d` onde `n_d = 1 + n_vertices_validos`. Se n_d < 4, mae nao splita (ganho de resolucao baixo demais).
- **Quebra de simetria angular aceita**: se vertices descartados nao sao simetricos (e.g., 0°, 60°, 120° mas nao 180°, 240°, 300°), centro de massa das filhas se desloca da posicao da mae. Drift residual aceito como tradeoff.

**Conservacao validada matematicamente:**
- Massa: n_d × (m_m/n_d) = m_m ✓
- Momentum linear: Σ m_d · u_d = n_d · (m_m/n_d) · u_m = m_m · u_m ✓
- Momentum angular: parcial (deslocamento do CM se vertices descartados nao sao simetricos).

**Predicoes v2.1 vs v2.0:**

| Metrica | v2.0 observado (t=12s) | v2.1 predicao |
|---|:---:|:---:|
| `pass_n_spawned` em t=0-13s | 287, 77, 21, 14, 7, 7, 14, 0, 0, ... | 50-100 sustentado/call |
| `max_v` pico | 8.25 (t=4.73s) | ≤ 0.5 |
| `a_pressure` pico | 332 (t=4.73s) | ≤ 10 |
| `mass_total` em t=13s | 101.43 (+0.35%) | 101.4 ± 0.5 (invariante) |
| Nucleo oco | sim | nao |
| "gap critico" no painel 3 | quase tudo | minoritario |

#### Pass N v2.1 — VALIDACAO PARCIAL t=0-40s (~77% de 50s) + 1 FALHA DIAGNOSTICADA (2026-05-20)

**Resultado experimental** (40s de 50s solicitados):
- ✅ **Sem explosao numerica**: `max_v` pico 0.43 (vs v2.0=8.25), `a_pressure` 3-4 (vs v2.0=332). Correcao ε=0.35 funcionou.
- ✅ **Massa quasi-invariante**: 101.08 → 102.49 em 40s (+1.4%, lento e por BiomassGrowth nao por Pass N).
- ✅ **Morfologia starfish excelente em t=38s** (frame 019): ~20 dendritos longos e finos, nucleo compacto pequeno bem visivel.
- ❌ **`pass_n_spawned` baixo apos t=6s**: 151 (call 1) → 0 por 14s → 5, 38, 50, 14, 13, 5, 6, 5, 6, 5 (intermitente).
- ❌ **Anel azul (gap critico) ao redor do nucleo + braços majoritariamente em "gap critico" no painel 3** (rho/rho0 < 0.7). Vazios persistem visualmente.
- ❌ `contrast_cs` colapso monotonico: 375 → 11 em 40s (sintoma classico §2.4).

**Diagnostico da falha (1, fundamental):**

**Trigger absoluto V_a > 1.5·dx² fica restritivo por geracao:**
- Trigger: `V_a = m_a / rho_a > 1.5·dx²` ↔ `rho_a < m_a / (1.5·dx²)`.
- Para gen 0 (m = m_0 = dx²): `rho_a < 0.667·rho_0` dispara → casa com painel visual (0.7).
- Para gen 1 (m = m_0/7): `rho_a < (m_0/7)/(1.5·dx²) = 0.095·rho_0` dispara → particula tem que estar quase totalmente isolada!
- **Apos primeiro split o trigger fica 7× mais estrito.** Particulas gen 1 com rho/rho_0 ∈ [0.1, 0.7] aparecem como "gap critico" no painel 3 mas o trigger NAO dispara.

**Causa raiz**: trigger por volume absoluto era incompativel com refinamento adaptativo multi-geracao. Toda particula refinada (gen ≥ 1) viraria invisivel ao detector exceto em casos extremos.

#### Pass N v2.2 — TRIGGER RELATIVO RHO/RHO_0 (2026-05-20)

**Mudanca em [main.py:101](main.py#L101) (constantes) e [main.py:399-422](main.py#L399-L422) (logica):**

| Parametro | v2.1 | v2.2 | Razao |
|---|---|---|---|
| `PASS_N_V_RATIO = 1.5` | absoluto (V_a > 1.5·dx²) | REMOVIDO | Trigger absoluto ficava 7× restritivo apos gen 1 |
| `PASS_N_RHO_TRIG = 0.7` | (não havia) | rho_a/rho_0 < 0.7 | Trigger RELATIVO — casa com painel visual "gap critico" |

Trigger novo: `colony_mask AND (rho_a/rho_0 < 0.7) AND (m > m_floor)`.

**Propriedades do novo trigger:**
- Independente da geracao — gen 0, 1, 2 disparam igualmente quando rho cai abaixo de 70% do nominal.
- Casa exatamente com a definicao visual de gap (rho < 0.7·rho_0 do painel 3).
- Prioriza mais esvaziadas (menor rho_rel) quando excede limite de mães por call.

**Predicoes v2.2 vs v2.1:**

| Metrica | v2.1 observado (t=40s) | v2.2 predicao |
|---|:---:|:---:|
| `pass_n_spawned` em t=20-40s | 5-50 intermitente | 50-300 sustentado/call |
| Anel azul ao redor do nucleo | sim | nao (anel preenchido por splits gen 2) |
| "gap critico" no painel 3 (t=40s) | maioria dos braços | minoria, isolado |
| `mass_total` em t=50s | 102.5 | 102.5 ± 0.5 (invariante) |
| Morfologia | ✅ ~20 dendritos starfish | ✅ preservada com partículas mais densas |
| `contrast_cs` em t=40s | 11 (colapso) | depende — talvez melhore com kernel mais denso |

#### Pass N v2.2 — VALIDACAO t=0-50s + FALHA DIAGNOSTICADA (2026-05-21)

**Resultado experimental** (50s solicitados, executou ate t=49.95s, 9200 iter):
- ❌ **`pass_n_spawned` cataclismicamente baixo**: 165 em iter 200 (t=3.1s) → **0 sustentado durante iter 400-3000 (~32s)** → retomada esporadica iter 3200+ com picos 7-21/call.
- ❌ **Nucleo permanece oco** em frame 23 (t≈50s): painel rho/rho_0 mostra circulos azuis (rho < 0.7·rho_0) persistentes no centro.
- ❌ **Braços como fios isolados**: rho/rho_0 < 0.7 espalha pelos dendritos sem reposicao.
- ✅ `mass_total` 101.08 → 103.01 (+1.9%) — consistente com conservacao por construcao (m_d = m_m/n_d); crescimento residual vem de BiomassGrowth.am.

**Diagnostico (2026-05-21) — falha geometrica: filhas gen-1 sao ESTRUTURALMENTE IMPOSSIVEIS de re-splittar.**

A cKDTree em [main.py:452](main.py#L452) (v2.2) era construida **sobre todas as particulas existentes, inclusive a propria mae sob teste**. Para cada vertice `vk` no padrao hexagonal:

```
d_existing = distância(vk, vizinho mais próximo na tree)
           ≤ distância(vk, mãe)
           = r_off = ε · h_mãe = 0.35 · h_mãe (por construção)
```

Filtro de proximidade `d_existing >= prox_min = 0.4·dx_0` exige:
```
0.35 · h_mãe > 0.4 · dx_0  ↔  h_mãe > 1.143 · dx_0
```

Comportamento por geracao (`h_0 = 1.8·dx_0`, `α = 0.6`):

| Geracao | `h_mãe / dx_0` | Passa `> 1.143·dx_0`? | Resultado |
|---|:---:|:---:|---|
| 0 | 1.80 | ✅ | Splita normalmente |
| 1 | 1.08 | ❌ | **Todos 6 vertices abortados → split silenciosamente recusado** |
| 2 | 0.65 | ❌ | (já bloqueada por m_floor = m_0/49) |

**Implicacao:** as ~700-1000 filhas gen-1 acumuladas no iter 200 ficam permanentemente "invisiveis" ao refinamento. Quando o nucleo continua perdendo vizinhos por migracao do rim (Marangoni + hard pin K.17), nada repoe. Os splits gen-0 esporadicos pos-t=39s sao apenas particulas do rim em estiramento — one-shot, nao progressivo.

**Vetor 2 (secundario) — Trigger `rho_rel < 0.7` é late-acting:** so dispara apos a particula perder ~3-4 vizinhos. Combinado com `r_growth=0.02` (motor lento), durante t=0-37s a colonia esta quasi-estatica e o trigger nao dispara — comportamento fisicamente correto, mas combinado com Vetor 1 significa que gen-1 do iter 200 permanecem nao-refinaveis enquanto o nucleo se esvazia silenciosamente.

**Vetor 3 (não-causal) — Hard pinning interage mas nao e o conflito direto:** pin em `rho_b ≥ 0.8` congela nucleo enquanto rim migra → nucleo perde vizinhos → `rho/rho_0 < 0.7` dispara em principio. Mas se as candidatas viraram gen-1 antes (iter 200), elas ficam barradas por Vetor 1. Hard pin gera o gap; Vetor 1 impede a solucao.

#### Pass N v2.3 — EXCLUIR MAES DA KDTREE (2026-05-21, IMPLEMENTADO)

**Fix unico** em [main.py:451-465](main.py#L451-L465): construir cKDTree apenas com particulas que NAO sao mães do split atual. As mães serão removidas em `remove_particles` ao final da call — mantê-las na árvore era inconsistente fisicamente (elas são "transparentes" para o filtro de colisão).

```python
tree_mask = np.ones(len(fluid.x), dtype=bool)
tree_mask[split_idx] = False
positions = np.column_stack([fluid.x[tree_mask], fluid.y[tree_mask]])
tree = cKDTree(positions)
```

**Predicoes v2.3 vs v2.2:**

| Metrica | v2.2 observado (t=50s) | v2.3 predicao |
|---|:---:|:---:|
| `pass_n_spawned` em t=10-37s | 0 sustentado | 50-200 por call (gen-1 do iter 200 desbloqueadas) |
| Nucleo oco em t=50s | sim | nao |
| "gap critico" rho_rel < 0.7 (t=50s) | núcleo + braços | minoria, isolado |
| `mass_total` em t=50s | 103.01 (+1.9%) | ≈ 103 ± 0.5 (invariante por construcao) |
| Morfologia ~20 dendritos | preservada | preservada com partículas mais densas |
| Risco: gen-2 com h=0.65·dx_0 | n/a | sub-amostragem do kernel — **monitorar max_v < 1.0** |

#### Pass N v2.3 — VALIDACAO PARCIAL t=0-13s + EFEITO COLATERAL: OVER-PACK DO CENTRO (2026-05-26)

**Resultado experimental** (50s solicitados, executou apenas t=12.98s em 7200 iter antes do diagnostico — dt colapsou):

**Lado positivo — fix v2.3 funcionou empiricamente:**
- iter 200 (t=3.19s): 139 splits ✅
- iter 400 (t=4.64s): **374 splits** ✅ (gen-1 do iter 100/200 entraram como mães — desbloqueio confirmado)
- iter 600 (t=5.29s): 39 splits
- iter 800 (t=6.41s): 4 splits
- Total 556 splits em t=0-6s (vs 165 em todo v2.2 — 3.4× mais).
- Tamanho HDF5: 5078→5141 KB (+63 KB ≈ ~500 particulas novas).
- Conservacao de massa: 101.08 → 101.43 em 13s (consistente com BiomassGrowth.am, Pass N nao adiciona massa).

**Lado negativo — splits concentrados no centro causaram travamento:**
- iter 1000+: 0 splits sustentado.
- Causa raiz: trigger `rho/rho_0 < 0.7` em t=0-6s só dispara onde ha gap real, que e SO no centro relaxando a Gaussiana inicial (rim ainda coeso, motor fraco com r_growth=0.02). 556 mães identificadas todas perto do centro → ~2200 daughters reocupam o centro → rho local sobe → trigger se autoamortece.
- Efeito colateral devastador no dt: filhas com h=1.08·dx_0 empacotadas a 0.378·dx_0 → kernel vê excesso de vizinhos próximos → pressão EOS sobe → CFL_force aperta dt.
- `a_pressure` cravado em **5.015** desde iter 800 (vs baseline 3.0 em v2.2 — alavanca §8 violada).
- Pace temporal: dt avg de 5e-3 (iter 0-1000) → 4.3e-4 (iter 2000-7200) — **12× menor**. Simulação rastejando.
- Extrapolacao: chegar a t=50s exigiria ~75 000 iter adicionais (~2-3h wall time).

**Bug secundário diagnosticado:** `mean_c_n > max_c_n` (1.0177 vs 1.0000 cravado) em iter 6000+. Matematicamente impossivel — implica particulas com c_n > 1.0 não-clampadas. Provavel: `add_particles` em PySPH alocando buffer com lixo > 1 para propriedade `c_n` não explicitamente no `data` dict. Baixa prioridade, mas vale verificar inicialização de daughters.

**Frame i7_032 (iter 6400, t≈12s) — confirmacao visual:** painel 1 (rho_b) mostra colonia starfish com ~12-15 protrusoes radiais (motor de Marangoni inicial). Painel 3 (rho/rho_0) mostra **mancha azul intensa no centro** correspondendo as ~556 filhas empacotadas la — refinamento visualmente localizado mas NAO morfologicamente util (centro nao precisava de mais resolucao, dendritos sim).

#### Pass N v2.3.1 — GATE SUPERIOR rho_b < 0.7 (Opção 2, 2026-05-26, IMPLEMENTADO)

**Fix em [main.py:117-130](main.py#L117-L130) + [main.py:427-434](main.py#L427-L434):**

Adicionado `PASS_N_RHO_B_MAX = 0.7`. `colony_mask` agora exige `rho_b ∈ [0.3, 0.7]` (zona de transicao/rim ativo) ao inves de `rho_b > 0.3` apenas. Particulas em `rho_b ≥ 0.7` (nucleo + shell adjacente ao hard pin K.17) **nao sao candidatas a split**.

**Trade-off explicito:** o gap inicial do nucleo da relaxacao da Gaussiana NAO sera repovoado por Pass N. Aceita-se este custo porque:
1. O nucleo nao precisa de resolucao adicional — esta hard-pinado (K.17), nao se move, nao gera morfologia.
2. Over-pack do centro em v2.3 inviabilizou o run por colapso de dt (12× menor).
3. Refinamento morfologicamente util e nos braços que esticam (rho_b ∈ [0.3, 0.7]).

**Predicoes v2.3.1 vs v2.3:**

| Metrica | v2.3 observado (t=13s) | v2.3.1 predicao |
|---|:---:|:---:|
| `pass_n_spawned` em iter 200-1000 | 556 (no centro) | 0-30 (rim ainda coeso, sem gaps) |
| `a_pressure` em t=10s | 5.015 cravado | ~3.0 (orcamento §8 ✓) |
| dt avg em iter 2000-7200 | 4.3e-4 | ~5e-3 (recuperado) |
| `pass_n_spawned` em t > 30s | n/a (run nao chegou) | 50-200 sustentado (braços esticando) |
| t=50s atingido em | ~75 000 iter (extrapolado) | ~9000 iter (baseline v2.2) |
| Nucleo oco no painel 3 | preenchido (over-pack) | sim — trade-off aceito |
| Braços com rho_rel < 0.7 em t=50s | n/a | minoritarios (preenchidos por splits) |

**Validacao pendente v2.3.1:**
1. Rodar `make run` com `total_sim_time = 50`.
2. Criterios de aceitacao:
   - `pass_n_spawned` em t=0-10s próximo de 0 (não há gaps no rim ainda).
   - `pass_n_spawned` em t=30-50s sustentado em 20-200/call (braços esticando).
   - `a_pressure` plateau ~3.0 (sem over-pack).
   - dt avg ~5e-3 ao longo do run (pace temporal restaurado).
   - Frame 23 (t≈50s) painel 3: BRAÇOS sem azul (preenchidos por splits); núcleo PODE ter azul (aceito por design).
3. Bug secundario `mean_c_n > max_c_n` — checar se persiste; se sim, adicionar inicializacao explicita de `c_n` (e outras propriedades flutuantes) no `data` dict de `daughters.add_particles`.

#### Pass N v2.4 — RECONSTRUCAO TEORICA Liu/Violeau (2026-05-28, IMPLEMENTADO — aguardando validacao)

**Motivacao:** sob analise teorica rigorosa baseada em Liu 2003 (*SPH: A Meshfree Particle Method*) e Violeau 2012 (*Fluid Mechanics and the SPH Method*) — agora referencias obrigatorias do projeto [T6] e [T7] — diagnosticou-se que v2.3.1 ainda violava tres principios de consistencia SPH (licoes #27 e #28). v2.4 corrige todos simultaneamente.

**Tres mudancas estruturais:**

**(a) Trigger por particao da unidade σ_a (Violeau §3.6, Liu §3.3.3) — substitui ρ_rel < 0.7:**

Adicionada equacao `KernelSum` em [src/equations.py](src/equations.py) e novo Group `equations_kernel_sum` em [src/scheme.py](src/scheme.py) entre `equations_pre` e `equations_main`:
```
σ_a = Σ_j (m_j/ρ_j) · W(r_aj, h_aj)
```
Group separado e obrigatorio: se `KernelSum` estivesse no mesmo Group de `SummationDensity`, `s_rho[s_idx]` seria parcial durante o loop (PySPH acumula durante o loop e finaliza apos post_loop do Group). Trigger `σ_a < 0.85` corresponde a ~15% erro nos operadores SPH (Violeau §3.6) e e INVARIANTE sob refinamento — gen 0/1/2 disparam pelo mesmo limiar. Diferente de `ρ_rel = ρ_a/ρ_0`, que mistura particao da unidade com fator de massa das filhas (licao #28).

**(b) α = ε = 0.35 (Feldman 2006 — razao otima):**

`PASS_N_ALPHA: 0.6 → 0.35`, mantido `PASS_N_EPSILON = 0.35`. Feldman 2006 deriva que `ε/α = 1` minimiza erro de densidade pos-split (< 5% em regime ideal). v2.3 usou `α=0.6, ε=0.35` (razao 0.58), causando over-pack do nucleo: filhas-vertice a `0.378·dx` caem no pico do kernel `W(r/h≈0.35)` → densidade pos-split ≈ 1.7×densidade pre-split → `a_pressure=5.015` sustentado → dt collapse 12× (licao #27).

Custo aceito: kernel das filhas agora com `h_filha = 0.35·h_mae = 0.63·dx` → suporte ~1.26·dx → ~12 vizinhos efetivos (vs ~35 em gen-0). Liu §3.3 reconhece que < 20 vizinhos degrada consistencia de ordem 1, mas Liu §6.5 argumenta que sub-amostragem e menos catastrofica que over-pack porque **nao viola conservacao** — apenas reduz acuracia local nas filhas.

**(c) n_d = 7 ESTRITO (Violeau §7.4.3 — simetria hexagonal):**

Substituida exigencia `n_d ≥ 4` (v2.3.1) por `n_d == 7` (v2.4) em [main.py:529](main.py#L529). Se qualquer vertice falhar no proximity guard, o split inteiro e adiado para proxima call. Violeau §7.4.3 prova que apenas distribuicao hexagonal **simetrica** preserva simultaneamente:
1. Erro de densidade < 5% (Feldman 2006)
2. Centro de massa na posicao da mae (conservacao espacial)
3. Tensor de inercia local (conservacao de spin)
4. Momento angular exato (Liu §3.4)

Splits parciais (n_d ∈ [4, 6]) com vertices descartados assimetricamente deslocam o CM e introduzem torque espurio que amplifica instabilidade de tracao (Liu §4.2.4).

**Implementacao:**
- [src/equations.py](src/equations.py) `KernelSum`: novo loop SPH para σ_a.
- [src/scheme.py](src/scheme.py): import `KernelSum`, novo `equations_kernel_sum` Group, retorno em 3 grupos.
- [main.py](main.py): propriedade `sigma_a` no array, inicializada em 1.0; constantes `PASS_N_SIGMA_TRIG=0.85`, `PASS_N_ALPHA=0.35`; trigger via `sigma_a < 0.85`; priorizacao por menor σ_a; `n_d == 7` estrito; `sigma_a` copiada para filhas.

**Conservacao validada (Liu §3.4, Violeau §5.3):**
- **Massa:** `7 · (m_m/7) = m_m` por mae ✓ (exato)
- **Momento linear:** `Σ_d m_d · v_d = 7 · (m_m/7) · v_m = m_m · v_m` ✓ (exato porque v_d = v_m)
- **Momento angular:** 6 vertices simetricos em torno do centro + 1 filha no centro → contribuicao rotacional nula em torno do CM da mae ✓ (exato pela exigencia n_d=7)

**Predicoes v2.4 vs v2.3.1 (a validar):**

| Metrica | v2.3.1 observado (t=13s) | v2.4 predicao |
|---|:---:|:---:|
| `pass_n_spawned` em t=0-10s | 556 (over-pack centro) | 0-10 (nucleo excluido + σ_a controlado) |
| `pass_n_spawned` em t=30-50s | n/a (run nao chegou) | 50-200 sustentado nos braços (σ_a cai com alongamento) |
| `a_pressure` plateau | 5.015 (vs orcamento <3) | ≤ 3.5 (densidade pos-split com erro ~5-8%) |
| dt avg | 4.3e-4 (12× degradado) | ~5e-3 (recuperado) |
| t=50s atingido em | ~75 000 iter (extrapolado) | ~9 000 iter (baseline v2.2) |
| Vizinhos por filha | ~35 (α=0.6) | ~12 (α=0.35) |
| Frame 23 (t≈50s) painel σ_a | n/a | braços σ_a > 0.85 (refinados); nucleo σ_a < 0.85 (excluido por gate) |

**Criterios de aceitacao v2.4:**
1. Sem regressao numerica: `max_v < 1.0` em todo run, `a_pressure < 4.0` plateau.
2. `mass_total` invariante (101.7 ± 0.5) — conservacao exata por construcao.
3. Spawn sustentado em t > 30s nos braços (σ_a cai quando dendritos esticam).
4. Frames mostram braços com densidade restaurada (σ_a > 0.85 visivel no plot, se instrumentado).
5. Motor preservado: `mean_v >= 0.0005`, `n_fast > 0` em t > 30s.

**Riscos identificados:**
- Sub-amostragem das filhas (h pequeno) pode degradar acuracia do `MarangoniForce` e `FlagellarForce` localmente — monitorar `a_marangoni` nas filhas via inspecao de outliers no log.
- Exigencia n_d=7 estrita pode reduzir spawn rate inicial (rim muito empacotado adia split ate gap > 2·prox_min). Aceito por design — splits prematuros causam over-pack.
- Bug secundario `mean_c_n > max_c_n` de v2.3 nao endereçado em v2.4 — checar se persiste no log.

**Validacao pendente:**
1. Rodar `make run` com `total_sim_time = 50`.
2. Verificar criterios de aceitacao acima.
3. Se passar t=50s, estender para `total_sim_time = 100`.

#### Pass N v2.4 — VALIDADO em t=50s sobre baseline T2g: ESTAVEL mas NET-NEGATIVE (2026-06-01)

Run completo t=50s (`α=ε=0.35, σ_trig=0.85, gate rho_b∈[0.3,0.7], gen≤1` via `m_floor=1/7`, `MAX_PARENTS=25`). **Numericamente saudavel, morfologicamente inutil, computacionalmente proibitivo.**

- **dt:** 203400 iter para 50s (vs ~2400 do T2g sem Pass N) — **~85× overall, ~130× no regime permanente** (dt ~1.6e-4 vs 0.021). Confirmado proibitivo (3.5h→5.4h wall).
- **Fisica saudavel:** `a_pressure` 3.0 plateau (1 spike de 40 em t=17s num split), `a_marangoni` 8-9 estavel, `mass_total` 103 conservada. Morfologia ~18-22 dendritos coerentes, SEM fragmentacao (diferente do T2d+PassN).
- **MAS sem ganho vs T2g:** `contrast_cs` termina ~15 = IGUAL ao T2g; `mean_v` 0.0005 < T2g 0.001 (mais congelada). Morfologia visualmente equivalente.
- **Raiz do zero-benefício:** só **98 mães splitaram (686 filhas) = +1.7% de partículas** — cobertura negligivel para densificar braços. Gate `σ_a<0.85` + `n_d=7 estrito` + `rho_b∈[0.3,0.7]` dispara pouquissimo.
- **BIND FUNDAMENTAL (lição #29):** o dt e fixado pela MENOR partícula (gen-1 h=0.034), **INDEPENDENTE da contagem**. 98 filhas cravam o dt tanto quanto 10000. Paga-se a penalidade de 85× inteira por ~0 refinamento útil — custo/benefício péssimo por construção. Confirma: **Pass N entrega resolução, não largura física** (estreitamento dos dedos e problema do MOTOR, não de refinamento).

#### Pass N v2.5 — REFINAMENTO DE PREENCHIMENTO REAL ("enchimento agressivo") (2026-06-02, IMPLEMENTADO — aguardando validacao)

**Decisao do usuario (NAO abandonar Pass N):** braços sem suporte de kernel violam a Particao da Unidade (`σ_a≈1`, Violeau §3.4) e a Consistencia de Interpolacao (Liu §3.3) — gradientes de Marangoni nos dendritos ficam numericamente espurios. A base cientifica da tese exige resolver isto. A v2.4 falhou por ser **excessivamente conservadora**, nao por estar errada.

**Quatro mudancas (vs v2.4):**
1. **`α: 0.35 → 0.75`** (`PASS_N_ALPHA`) — h_filha proximo da mae, **preserva dt-por-h** (penalidade ~1.5× em vez de ~85×). **ε mantido em 0.35** (decoupling deliberado de Feldman, ver §10 amendment). **Aceita-se o over-pack inicial**, confiando na EOS coesiva (`tension_ratio=0.30`) + Monaghan (`alpha_mon=0.12`, K.23) para relaxar (compromisso Liu §6.5).
2. **`σ_trig: 0.85 → 0.95`** (`PASS_N_SIGMA_TRIG`) — gatilho preemptivo: agir ao primeiro sinal de estiramento (5% de perda de consistencia, nao 15%).
3. **gate `rho_b: [0.3,0.7] → [0.15,0.95]`** (`PASS_N_RHO_B_MIN/MAX`) — Pass N livre para fechar buracos em quase toda a colonia, incluindo a junção núcleo-dendrito.
4. **`MAX_PARENTS: 25 → 100`** — enchimento agressivo, agora VIÁVEL porque α=0.75 mantem h grande (lição #29: dt e min-h, não contagem; com h grande a contagem deixa de cravar o dt).
- `n_d=7 estrito` MANTIDO (Violeau §7.4.3 — CM + momento angular). `gen≤1` mantido (`m_floor=1/7`).

**RISCO PRINCIPAL (lição #27 re-apostada):** α=0.75 com ε=0.35 dá `ε/α=0.47` — MAIS longe do otimo Feldman do que a v2.3 (0.58) que travou em t=13s. Efeito contraintuitivo: α maior melhora dt-por-h mas PIORA dt-por-força via over-pack. Efeito líquido no dt e incerto — pode re-travar. **Critério #2 (dt avg t>30s ≥ 0.5× inicial) e o teste direto de se o over-pack relaxou.** Fallback se falhar: ε→0.5 (reduz over-pack mantendo parte do ganho de h) ou α→0.6.

**Criterios de Sucesso v2.5 (obrigatorios para desbloqueio):**
1. **Consistencia geometrica:** `σ_a ≥ 0.90` sustentado na zona de transição e braços.
2. **Continuidade de dt:** dt avg em t>30s ≥ 0.5× o dt avg inicial (t<5s) — guarda contra over-pack não-relaxado.
3. **Morfologia:** braços densos e contínuos (fluido), NÃO correntes de partículas isoladas.
4. **Conservacao:** `mass_total` varia < 0.1% por splitting (crescimento via `r_growth` permitido).

**Validacao pendente:** `make run` com `total_sim_time=50`. Se relaxar (criterio #2 ✓) e densificar (criterio #3 ✓), estender t=100s. Se re-travar (lição #27), aplicar fallback ε=0.5.

#### Pass N v2.5 — AUDITORIA contra referencias + 3 correcoes ANTES de religar (2026-06-15)

A pedido do usuario, auditou-se o codigo do split v2.5 ([main.py:547-714](main.py#L547-L714)) contra Vacondio 2013 / Feldman 2006 / Soleimani §3.2.6 / Liu / Violeau ANTES de `use_pass_n=True`. **Nucleo do algoritmo FIEL** (1 mãe→7 filhas hexagonal, n_d=7 estrito Violeau §7.4.3, massa/momento linear+angular conservados Liu §3.4, velocidade herdada, escalares copiados, mãe deletada, mãe excluida da KDTree lição #26). **Tres problemas corrigidos:**

1. **BUG `is_filler` (corrigido):** o dict das filhas NAO setava `is_filler` → herdava lixo do realloc → filha podia nascer `is_filler>0.5` e ser PINADA como filler congelado (scheme.py), matando o refinamento. Mesmo risco do bug `mean_c_n>max_c_n` da v2.3. Fix: `"is_filler": [0.0]*7` (filhas = biomassa real refinada). Outros campos novos (`shift_*`, `L*`, `M*`) sao recomputados a cada passo antes de usar → lixo inofensivo; **so `is_filler` e estado persistente** lido por BiomassGrowth/SurfactantEquation/CustomEulerStep.
2. **Comentarios/print defasados (corrigido):** diziam `gen≤2 / m₀/49` mas a constante real e `1/7` = `gen≤1`. Print dizia `ε=α={eps}` (mentira — α=0.75≠ε=0.35) → corrigido p/ `ε={eps}, α={alpha}`.
3. **Comentario do void_mask C3 (corrigido):** ainda descrevia o C4 revertido.

**Dois DESVIOS de referencia sinalizados com ⚠️ no codigo (em alerta p/ seguir a referencia se nao atender):**
- `PASS_N_SIGMA_TRIG = 0.95` (engenharia preemptiva, ~5% erro) vs **canonico 0.85** (Violeau §3.6, ~15% erro). SE over-refino/dt travado → mudar p/ 0.85.
- `ε/α = 0.47` (ε=0.35, α=0.75) vs **otimo Feldman ε/α=1**. Autorizado pelo §10 amendment (over-pack p/ preservar dt, Liu §6.5) MAS condicionado ao criterio runtime "dt avg t>30s ≥ 0.5× inicial". RISCO #27: 0.47 e MAIS distante do otimo que a v2.3 (0.58) que travou. SE falhar → ε→0.5 ou α→0.6 (rumo à referencia); alvo Feldman ε=α custa dt (#29), so com timesteps individuais (Fase 2).

**Estado:** codigo auditado e corrigido, `use_pass_n` AINDA False. Plano de 2 fases em `/Users/isa/.claude/plans/linear-dazzling-comet.md` — Fase 1 = validar v2.5 (com `use_insert=False` p/ isolar, §2.3); Fase 2 = fork block-timestep [T9] so se v2.5 over-packar. **Achado de viabilidade:** timesteps individuais fieis sao incompativeis com o PySPH (dt global = min(h); `compute_accelerations` monolitico; `one_timestep` transpilado) → exigiriam forkar o `acceleration_eval`. v2.5 (α grande, h grande) e a saida nativa.

#### Pass N v2.5 — Fase 1 (use_pass_n=True, use_insert/kgc=False) FALHOU Critério #2: dt colapsou 97× (2026-06-15)

Run da Fase 1 parou em **t=22.5s (23%)** — o usuario achou que "afogou", mas **NAO e afogamento**. Quimica saudavel: `mean_cs`=0.023 (baixo; afogamento seria ~0.13), `contrast_cs`=21 (sem colapso), `mass_total`=101.7 (+0.6%, conservada). Frame: halo cs limpo e radial, ~15-20 dendritos incipientes (estagio inicial t=22).

**É COLAPSO DE dt (~97×):** dt 1.8e-2 (t<15s) → 1.8e-4 (pós-t≈19.5s). A sim virou lesma; t=100 exigiria ~500k iter. **Critério #2 do plano (dt avg ≥ 0.5× inicial) FALHOU** (dt ficou 0.01× inicial). `a_pressure`=3.0 e `max_v`=0.07 normais (NAO foi a explosao de pressao da v2.3).

**Duas causas raiz (assinatura no HDF5 final):**
1. **Over-pack (#27, exatamente o que o alerta ε/α=0.47 previu):** `rho/rho0` ate **2.94** (p99=1.72), 3.8% das particulas >1.3. As 7 filhas a `ε·h_m=0.35·h_m` com `h=0.75·h_m` → distancia `ε/α=0.47·h_filha` → dentro do kernel umas das outras → densidade inflada.
2. **BUG NOVO — limite gen≤1 por massa derrotado pelo BiomassGrowth (lição #37):** `h_min`=0.0544=`0.75²·h0` e `min m`=5.9e-5=`m₀/49` → **gen-2 espúrio** apesar do "gen≤1". O gate `m > m₀/7` assumia que a massa so diminui ao dividir, mas `BiomassGrowth` cresce `m` de volta acima do floor → filha re-splita → gen-2,3… → h encolhe → dt cai (#29). **Limite de geração por massa NÃO funciona quando há crescimento.**

**Fix aplicado (item 1, NAO religado ainda):**
- **Contador `gen` explícito** (propriedade nova, init 0; filha = mãe+1; gate `gen < PASS_N_MAX_GEN=1`). Imune ao crescimento — substitui o `m > m_floor`. `PASS_N_M_FLOOR_RATIO` e `V_0` removidos.
- **`ε`: 0.35 → 0.5** (ε/α 0.47→0.67, rumo a Feldman) p/ reduzir over-pack. RISCO #26: r_off=0.9·dx pode abortar splits no domínio empacotado → monitorar `pass_n_spawned`; se cair a ~0, reduzir `α` (0.75→0.6) em vez de subir `ε`.
- Alertas ⚠️ e comentarios-resumo atualizados; v2.4→v2.5 no print.

**Bind estrutural reconhecido:** over-pack (ε/α baixo) vs abort (#26, ε alto) vs dt-cost (#29, α baixo) é um trilema. Se o re-teste do item 1 ainda colapsar dt, é a evidência definitiva de que a rota nativa (α grande) não escapa → Fase 2 (fork) OU Opção 1 (KGC, ja suficiente cientificamente pela medição #2).

#### Pass N v2.5 — Fase 1 RE-TESTE (gen-counter + ε=0.5): gen-fix OK mas over-pack↔dt é BIND FUNDAMENTAL (2026-06-15)

Re-rodada com o fix do item 1 (contador `gen` + ε 0.35→0.5). **gen-fix FUNCIONOU:** `h_min`=0.0726=`0.75·h0` EXATO → só gen-1, **ZERO gen-2** (antes 0.56·h0); massa 101.1→102.3 (+1.2%, conservada). Colapso de dt ADIADO t≈19.5→**t≈34** → run foi de 23%→**37%**. Quimica saudavel de novo (`mean_cs`=0.027, `contrast_cs`=18) — NAO afogamento. Frame t=36.7: ~18-20 dendritos limpos, nucleo compacto.

**MAS o over-pack ainda colapsou o dt ~35×:** `rho/rho0` ate **3.31** (p99=1.95), 8.2% >1.3. Com `h_min` só 1.33× menor (gen-1), os 35× sao DOMINADOS pelo over-pack (pressao/CFL nos clusters), nao por h. ε=0.5 (ε/α=0.67) reduziu mas NAO eliminou.

**BIND FUNDAMENTAL comprovado (lição #38):** refinar 7× = 7 filhas na area de 1 mãe → espacamento ~0.38·dx. Para NAO over-packar (Feldman ε/α≈1) precisaria `h_filha≈0.21·h0` → `α≈0.21` → h minusculo → dt colapsa por CFL (#29). Manter h grande (α=0.75) força ε/α<1 → over-pack → dt colapsa por pressao. **Em framework de dt GLOBAL, refinar 7× SEMPRE colapsa o dt — por h pequeno OU por over-pack. Nao ha calibracao de ε/α que escape** (ε/α=1 exigiria α≈0.2). Por isso os timesteps individuais [T9] (Fase 2) sao o fix "fisicamente correto" — deixam filhas de h pequeno integrarem no rung curto delas sem arrastar o resto.

**Veredito:** Fase 1 (rota nativa α grande) FALHOU o Critério #2 de forma IRREMEDIAVEL (2 runs). gen-fix removeu uma causa (gen-2) mas o over-pack↔dt e estrutural. Decisao: **Opção 1 (KGC, recomendado — frontier nao e problema cientifico pela medição #2)** OU **Fase 2 (fork block-timestep, semanas, só se resolucao real exigida)**. Item 3 (ε=α=0.6 Feldman ótimo) provavelmente futil (#26 abort + dt-cost). **Discussao pendente (2026-06-15): literatura sobre o problema + como tratar o VÁCUO VISUAL se optar pela Opção 1 (renderizacao de campo reconstruido / "sombras" temporais — viz, NAO solver).**

#### DECISÃO: Opção 1 ADOTADA (KGC + inserção frozen) — config de baseline (2026-06-15)

Flags: `use_kgc=True`, `use_insert=True` (C3.4 frozen-só), `use_pass_n=False`, `use_shift=False`. Pass N e Shifting preservados via flag (Fase 2 futura). main.py limpo (histórico migrado p/ cá). Divisão do vácuo: **KGC corrige o operador ∇cs na frontier** (tipo-B, pontas σ_a=0.98) + **inserção frozen preenche o vácuo estrutural** (tipo-C, núcleo/junção).

**DEFESA DE LITERATURA — KGC resolve o vácuo e NÃO é erro de fundamentação:**

O vácuo degrada o **operador SPH**, não a física. A estimativa de gradiente SPH `∇A_i ≈ Σ_j V_j (A_j−A_i) ∇W_ij` só reproduz `∇A` exatamente se a **condição de consistência de 1ª ordem** `Σ_j V_j (x_j−x_i)⊗∇W_ij = I` valer. Com suporte completo/regular ela vale; com suporte INCOMPLETO (superfície livre, braços esparsos) ela falha → o gradiente puro tem erro O(1) (inconsistência de ordem zero). KGC/CSPM multiplica por `L = M⁻¹` (M = a matriz acima) e **RESTAURA a condição exata**: o gradiente corrigido reproduz campos lineares EXATAMENTE, independente da completude do suporte. É **restauração de consistência matematicamente rigorosa**, não um remendo.

Por que NÃO "esconde" o vácuo: KGC não finge que há partículas onde não há — calcula o MELHOR gradiente 1ª-ordem-consistente a partir das partículas que existem. A força de Marangoni precisa de `∇cs`; KGC entrega o `∇cs` correto para a distribuição real. A física (`F=−β∇cs`) é inalterada; só a *estimativa numérica* de `∇cs` é corrigida. O auto-gating por `det(M)`: bulk (det≈1) → L≈I (não toca); frontier (det<1) → corrige exatamente onde o suporte é incompleto. Empiricamente (medição #2): pontas σ_a=0.98 → L≈I (KGC quase não altera onde já está bom); a amplificação `a_marangoni` 8→14 é a correção agindo onde o suporte é deficiente.

**Referências bibliográficas:**
- **Bonet & Lok (1999)**, *CMAME* 180:97-115 — formulação KGC; gradiente corrigido `L=(Σ V_j(x_j−x_i)⊗∇W)⁻¹` restaura consistência linear com distribuição irregular.
- **Liu, M.B. & Liu, G.R. (2006)**, *Applied Numerical Mathematics* 56:19-36, "Restoring particle consistency in SPH" — referência canônica: trata DIRETAMENTE a inconsistência por suporte incompleto (fronteira/superfície livre) e prova que formulações corretivas restauram consistência C0/C1. **É exatamente o nosso caso (braços = suporte truncado tipo superfície-livre).**
- **Chen, Beraun & Carney (1999)**, *IJNME* 46:231-252 — CSPM, consistência perto de fronteiras/suporte incompleto.
- **Randles & Libersky (1996)**, *CMAME* 139:375-408 — SPH normalizado/corrigido para consistência.
- **Oger et al. (2007)**, *JCP* 225:1472-1492 — gradiente renormalizado (Bonet-Lok) é exato para campos lineares, melhora convergência.
- **Liu & Liu (2003)** [T6] §3.3 (inconsistência por deficiência de fronteira; CSPM); **Violeau (2012)** [T7] §3.4-3.6 (consistência discreta; σ_a = medida do erro).

**Caveat honesto (documentado):** a forma gather corrigida `L_i·∇W` NÃO é pairwise-antisimétrica → conservação de momento é aproximada. Origem do trade-off: **Bonet & Lok (1999)** — o título é "Variational and **momentum preservation** aspects of SPH" — mostram que a correção ingênua sacrifica conservação e derivam a forma variacional simétrica que a restaura; princípio geral (antissimetria⟺momento) em Liu §3.4 / Violeau §5.3 / Monaghan 1992. Aceitável aqui porque a colônia é hard-pinada + drag-dominada → momento já não é conservado (pin e drag são forças externas). Sem erro escondido.

**VALIDAÇÃO OBJETIVA (teste-ouro, 2026-06-15):** [test_kgc_consistency.py](test_kgc_consistency.py) replica a matemática do solver (KernelGradientCorrection + L·DWIJ em MarangoniForce) com o mesmo CubicSpline e mede ∇ de um campo LINEAR (∇ exato conhecido) em suporte completo vs. truncado. Resultado: **KGC reproduz o gradiente a ~1e-16 (precisão de máquina) em TODO ponto** — interior E bordas sub-resolvidas (det(M) de 0.99 até 0.34) — enquanto o SPH puro erra 1% (interior, ordem-zero M_xx=0.988≠1) → 62% (borda). Isto é a definição matemática de consistência de 1ª ordem restaurada (Liu §3.3); prova que a implementação é FIEL a Bonet-Lok 1999, não remendo post-hoc (fórmula errada não reproduziria campo linear a 1e-16 sob suporte truncado). As 3 escolhas de engenharia (det_min=0.25, escopo no ∇cs, forma gather) não afetam a correção — teste roda com det_min=0.25 e KGC é exata até det=0.34.

**Vácuo VISUAL (separado da física):** mesmo com KGC, o painel `rho/rho0` (scatter) mostra buracos. Solução = RENDERIZAÇÃO (pós-processamento, não solver): (a) campo SPH reconstruído + interpolante de Shepard/σ-normalizado (**Price 2007, splash, PASA 24:159**) — buracos somem por interpolação de kernel; (b) ocupação/"sombra" temporal `∫ρ_b dt` — casa com a pegada acumulada da `reference.jpg` (que é ela própria um footprint integrado da frente de swarming). AMBOS puro pós-processamento — pôr "sombras" no solver = C4 #34/#31.

---

### Rotas para o problema do vacuo — analise comparativa (2026-06-02)

Apos v2.4 (net-negative, lição #29) e antes de apostar tudo na v2.5 (over-pack, lição #27), fez-se um levantamento da literatura SPH ([T8]/[T9]/[T10], §3.0.1). **Decisao do usuario: NAO abandonar o Pass N** — braços sem suporte de kernel violam `σ_a≈1` (Violeau §3.4) e a consistencia de interpolacao (Liu §3.3), tornando o `∇cs` da Marangoni nos dendritos numericamente espurio. Isso invalidaria a comparacao com Trinschek (T1). **Mas o objetivo (`σ_a→1`) ≠ mecanismo (splitting):**

| Rota | Mecanismo | Custo dt | Ancoragem | Status |
|------|-----------|:--------:|-----------|--------|
| **A** | Particle Shifting (redistribui, `δr=-D∇C`) | **ZERO** | [T8] Lind/Xu/Adami; Liu §6.4/§6.5 | **ESGOTADA** — v1 congela (#31), A.3 deixa vacuo (#32) |
| **B** | Timesteps individuais/locais | controlado | [T9] Springel; Saitoh-Makino | futura (reescreve integrador) |
| **C** | Kernel Gradient Correction (KGC/CSPM) | **ZERO** | [T10] Bonet-Lok; Liu §3.3 (CSPM ∈ T6) | **EM IMPLEMENTACAO 2026-06-02** |
| **D** | Pass N "do paper" (mae passiva + grad-h) | alto | Vacondio/Barcarolo/Chiron; Springel-Hernquist | auditoria pendente |

Sucessos preservados: a maquinaria Vacondio/Feldman do Pass N (n_d=7, conservacao exata, trigger σ_a) esta correta; T2g e baseline estavel ate t=100s; lição #29 alinha com a literatura de APR (resolucao fina = timestep curto = caro).

#### Rota A — Fickian Particle Shifting (IMPLEMENTADO 2026-06-02, aguardando validacao)

Preenche o vacuo REDISTRIBUINDO particulas existentes — custo de dt ZERO (nao encolhe h). Ataca o vacuo como problema de DISTRIBUICAO, nao de contagem.

**Implementacao:**
- **[src/equations.py](src/equations.py) `ParticleShift`:** `loop` acumula `∇C = Σ_j (m_j/ρ_j)∇W_ij` em `shift_dC_x/y`; `post_loop` calcula `δr = -shift_coeff·h²·∇C`, com cap `|δr| ≤ shift_cap·h` (Lind 2012 — shifts grandes desestabilizam). **Gate = inverso do hard pin K.17:** shift=0 se `rho_b≥0.8` (nucleo) OU `c_n<0.6` (starvation) OU `rho_b<0.05` (agar). Auto-limitante: uniforme → ∇C=0 → δr=0.
- **[src/scheme.py](src/scheme.py):** Group `equations_shift` separado apos `equations_main` (condicional a `use_shift`); `CustomEulerStep.stage1` aplica `d_x += shift_x; d_y += shift_y` apos o passo de velocidade (shift=0 para pinados → soma incondicional segura).
- **[main.py](main.py):** flags `use_shift=True`, `SHIFT_COEFF=0.5`, `SHIFT_CAP=0.05`; propriedades `shift_dC_x/y`, `shift_x/y`. **`use_pass_n=False` (Pass N preservado, desligado para isolar a Rota A).**

**SIMPLIFICACAO v1 (documentada):** shift aplicado so a posicao; campos (`rho_b, cs, c_n`) acompanham a particula SEM a correcao de Taylor `δr·∇φ` de Lind 2012. Justificado por shift pequeno (cap 0.05h), `rho` recomputado a cada passo (SummationDensity), cs/c_n difusivos. TODO: adicionar `δr·∇φ` se houver drift.

**RISCO PRINCIPAL — erosao morfologica:** PST e um HOMOGENEIZADOR (criado para escoamentos de bulk uniformes); a morfologia dendritica e ANTI-homogenea (braços finos + baias vazias). Mitigantes: (a) baias sao AGAR particle-filled (rho_b<0.05, gated OUT) → ∇C pequeno na fronteira braço-baia → contorno preservado; (b) ∇C grande so onde particulas estao genuinamente esparsas (braços esticados, junção). **Monitorar nos frames se os braços alargam/borram** — se sim, reduzir `SHIFT_COEFF`/`SHIFT_CAP` ou estreitar o gate.

**Quando religar Pass N junto com shift:** `add_particles` das filhas precisara incluir `shift_x/y/dC` no dict `data` (senao herdam lixo, cf. bug `mean_c_n>max_c_n` v2.3).

**Criterios de sucesso Rota A:**
1. **Vacuo preenchido:** painel 3 (rho/rho0) sem gaps `<0.7` sustentados nos braços e na junção núcleo-rim.
2. **dt preservado:** dt avg ~ T2g (~0.02; **sem** os 85× do Pass N) — confirma custo-ZERO.
3. **Conservacao:** `mass_total` identica ao T2g (shift nao cria/destroi massa).
4. **Morfologia NAO erodida:** braços continuam finos e separados (NAO blob homogeneo); `contrast_cs` nao colapsa abaixo do T2g.
5. **Estabilidade:** `a_pressure` ≤ 4, sem particulas ejetadas.

#### Rota A v1 (gate [0.05,0.8)) — SUCESSO NUMERICO / FALHA MORFOLOGICA (2026-06-02, t=50s)

`SHIFT_COEFF=0.5, SHIFT_CAP=0.05, gate rho_b∈[0.05,0.8)` (isotropico, incluindo o frontier motil).

- ✅✅ **dt = 0.0256 ≈ T2g — custo ZERO confirmado** (1932 iter para t=50, vs 203400 do Pass N). Promessa central cumprida.
- ✅ `a_pressure`=3.0, `mass`=102 conservada, **painel 3 LIMPO** (vacuo preenchido, σ_a→1 atingido). Objetivo numerico alcancado.
- ❌❌ **MORFOLOGIA SUPRIMIDA (critério #4):** colonia **congelou em disco lumpy r≈1.0** (vs T2g r≈2.5-3), ~12-15 bumps curtos, **ZERO dendritos**. `mean_v≈0.00008` (10× < T2g) — colonia parada.
- ⚠️ `contrast_cs`=37 = **falso-positivo** (borda nitida de disco compacto, NAO pontas dendriticas; `mean_cs`=0.013 baixo só porque a colonia e pequena). Lição §11 confirmada: nunca ler contrast_cs isolado.

**Causa-raiz (lição #31):** PST e um HOMOGENEIZADOR; fingering e ANTI-homogeneo. O shift isotropico no rim (a) **amortece as perturbacoes curtas** que SEMEIAM os dedos (Mullins-Sekerka nunca nucleia), (b) redistribui de volta o push da Marangoni → trava a expansao. **A homogeneizacao que preenche o vacuo e a mesma que mata o fingering.** Confirmado o risco sinalizado ex-ante. Lind 2012 [T8] previne isto suprimindo a componente do shift NORMAL a superficie livre.

#### Rota A.3 (gate interior [0.6,0.8)) — IMPLEMENTADO 2026-06-02, aguardando validacao

Alavanca unica: `SHIFT_RHO_B_MIN: 0.05 → 0.6` ([main.py](main.py), threaded → scheme → `ParticleShift.rho_b_min`). Mantem `coeff=0.5, cap=0.05` (numericamente OK na v1). Desloca SO a banda interior `[0.6, 0.8)` — a junção núcleo↔rim onde vive o vacuo — **EXCLUINDO todo o frontier motil** (rho_b<0.6, que cobre o gate flagelar [0.1,0.6] e a zona de swarmers). As pontas crescem livres; só o interior estrutural e homogeneizado. Versao simplificada do tratamento de superficie livre de Lind 2012 (em vez de projetar a componente normal, exclui a faixa de rho_b da interface ativa).

**Predicao:** dt continua ~T2g (zero custo); junção núcleo↔rim preenchida (painel 3 limpo SO no interior); **dendritos voltam a formar** (frontier intocado) → `mean_v` recupera a ~0.001 (T2g), morfologia ~T2g ou melhor. Se FALHAR (colonia ainda congela OU vacuo da junção persiste), a homogeneizacao por shift e incompativel com esta morfologia → **partir para Rota C (KGC)** ou A.2 (supressao da componente normal via ∇rho_b).

**RESULTADO A.3 (2026-06-02, t=50s) — FALHOU NO OBJETIVO (vacuo persiste):**
- ✅ dt=0.0218 (custo ZERO), ✅ morfologia restaurada (~18-20 dendritos a r≈2.5, `mean_v`=0.0008 ≈ T2g, 10× > v1), ✅ `a_pressure`=3.0, ✅ massa conservada (101→102.9), ✅ `contrast_cs`=15.6 (saudavel, NAO o falso-positivo 37 da v1).
- ❌ **CRITÉRIO #1 (vacuo) FALHOU:** painel 3 mostra aglomerado de gaps (`rho/rho0<0.7`) **ainda concentrado na junção núcleo↔rim**; painel 1 mostra anel oco entre o núcleo central e as bases dos dendritos. **O vacuo — alvo da Rota A — persistiu.**
- ⚠️ Violacao §11 cometida na 1a leitura: declarei "sucesso" pela morfologia sem checar a metrica-alvo (vacuo). Corrigido apos feedback do usuario.

**Causa-raiz (lição #32):** o vacuo da junção esta na MESMA faixa de `rho_b` (<0.6) que o frontier motil que A.3 protege. As particulas esparsas do anel oco tem `rho_b<0.6` → **excluidas pelo gate [0.6,0.8)**. v1 [0.05,0.8) preenche o vacuo MAS congela; A.3 [0.6,0.8) preserva morfologia MAS deixa o vacuo. **`rho_b` nao separa "vacuo a preencher" de "frontier a proteger" — coabitam <0.6.** Camada mais profunda: o vacuo e artefato CINEMATICO do hard pin K.17 (núcleo congelado + rim expandindo → junção esticada regenera o gap a cada passo); shift de banda fina nao acompanha.

**Shifting (Rota A) ESGOTADO** — duas falhas independentes: v1 (homogeneiza→mata fingering, #31) e A.3 (vacuo co-localizado com frontier, #32). **Pivot autorizado para Rota C (KGC).**

#### Rota C — Kernel Gradient Correction (KGC/CSPM) — IMPLEMENTADO 2026-06-02, aguardando validacao

Bonet & Lok 1999 / CSPM ([T10], Liu §3.3 — CSPM e de Chen-Beraun/Liu-Liu, dentro do T6). **Corrige o operador `∇cs` da Marangoni SEM mover/criar particula** → restaura consistencia de 1a ordem nos braços sub-resolvidos (∇ de campo linear exato mesmo com vizinhanca incompleta). Ataca a justificativa cientifica real (∇cs espurio invalidaria Trinschek) em vez de tentar preencher o vacuo cinematico.

**Por que KGC escapa dos becos da Rota A:**
- **Nao move particula** → impossivel congelar a morfologia (lição #31 nao se aplica).
- **Auto-gateia pelo determinante** → nao depende de separar vacuo↔frontier por `rho_b` (lição #32 nao se aplica). Bulk: `det(M)≈1` → L≈I (sem correcao). Rim/vacuo: `det<1` → corrige.
- **Custo de dt ZERO** (so algebra local).

**Mecanismo:** `M_i = Σ_j V_j (x_j-x_i)⊗∇_iW_ij` (→ I p/ vizinhanca completa); `L_i = M_i^{-1}`; gradiente corrigido `∇cs_i = L_i·Σ_j V_j(cs_j-cs_i)∇W_ij`. Em `MarangoniForce.loop`, `DWIJ` e substituido por `L_i·DWIJ`.

**Implementacao:**
- **[src/equations.py](src/equations.py) `KernelGradientCorrection`:** `loop` acumula M (2×2); `post_loop` inverte → L (Lxx,Lxy,Lyx,Lyy). **Fallback CSPM:** `det(M)<det_min` (0.25) → L=identidade (reverte a SPH padrao; evita amplificar forca em particula quasi-isolada; limita `|L|≲4×`).
- **[src/scheme.py](src/scheme.py):** Group `equations_kgc` SEPARADO entre `kernel_sum` e `main` (L precisa estar invertida antes de MarangoniForce usar). Condicional a `use_kgc`.
- **[src/equations.py](src/equations.py) `MarangoniForce.loop`:** `cdwij = L_i·DWIJ`, usado no lugar de DWIJ. L inicia identidade → KGC off ⇒ SPH padrao (toggle limpo).
- **[main.py](main.py):** `use_kgc=True`, `KGC_DET_MIN=0.25`; props M/L (L init identidade). **`use_shift=False`** (Rota A esgotada, isola KGC).

**Escopo v1:** KGC aplicado SO ao `∇cs` da `MarangoniForce` (motor central, a preocupacao cientifica). `FlagellarForce` (tambem usa ∇cs) e `BiomassGradient` ficam para follow-up se a v1 validar.

**RISCO:** em regioes muito degradadas, `L~1/det` amplifica o gradiente → pode spike de `a_marangoni`/instabilidade. Mitigado pelo fallback `det_min=0.25` (|L|≲4×). Monitorar `a_marangoni` e `a_pressure` — se spikarem, subir `det_min`.

**Criterios de sucesso Rota C:**
1. **Operador valido:** `a_marangoni` consistente (sem o ruido espurio dos braços sub-resolvidos); idealmente instrumentar `σ_a`/erro do gradiente.
2. **dt ZERO custo:** dt avg ~ T2g (~0.02).
3. **Estabilidade:** `a_pressure` ≤ 4, `a_marangoni` sem spikes (fallback funcionando), sem particulas ejetadas.
4. **Morfologia:** dendritos preservados ou MELHORADOS (gradiente correto → Marangoni mais fisica nas pontas); `mean_v`, `contrast_cs` ≥ T2g.
5. **Conservacao:** `mass_total` ~ T2g.

**Nota de conservacao:** a forma gather corrigida (`L_i·DWIJ`) NAO e pairwise-antisimetrica → perde conservacao exata de momento (Bonet-Lok teriam a forma variacional simetrica). Aceitavel aqui: a colonia ja e hard-pinada + drag-dominada (momento nao e conservado de qualquer forma). Documentar se houver drift de CM.

**RESULTADO Rota C / KGC (2026-06-02, t=50s) — SUCESSO no proposito, vacuo central persiste (esperado):**
- ✅ dt=0.0208 (custo ZERO), ✅ `a_pressure`=2.99 plateau (fallback `det_min` segurou — sem spikes), ✅ massa conservada, ✅ morfologia (~18-20 dendritos, braços OK).
- ✅✅ **`a_marangoni` subiu de ~8 (T2g/A.3) → 13-14** — o gradiente corrigido **amplificou ~75% a forca nas pontas sub-resolvidas**: o `∇cs` que o SPH padrao SUBESTIMAVA no rim agora e consistente de 1a ordem (Liu §3.3). **KGC cumpriu seu objetivo** (operador valido nas pontas) de graca e estavel. **KGC e KEEP.**
- ❌ vacuo **CENTRAL** persiste — mas KGC nunca preencheria (corrige operador, nao cria particula). O vacuo que sobra e a **junção núcleo↔rim, colada no núcleo congelado** — NAO nos braços (que estao bons).

#### Mecanismo do vacuo central + analise preditiva de solucoes (2026-06-02)

**Causa-raiz:** o vacuo central e **evacuacao dinamica por orcamento de particulas**, NAO inicializacao. Massa cresce ~2% (`r_growth=0.02`) mas area ~25× (r 0.5→2.5); particulas ≈ constantes (Pass N off) → o shell **migra para os braços** (Marangoni) evacuando o centro; o núcleo (`rho_b≥0.8`) esta congelado (pin) e o crescimento aumenta `rho_b` no lugar SEM criar particulas. O vacuo abre na fronteira pin↔shell e **cresce no tempo** (t10 pequeno → t50 claro). Confirmado dinamico.

**Implicacao:** o vacuo e (a) dinamico (regenera) e (b) deficit de particulas. Solucoes estaticas (inicializacao) NAO previnem.

| Familia | Predicao sobre o vacuo | Veredito |
|---|---|---|
| Iniciar sem núcleo pinado / pre-densificar | atrasa, EOS equaliza, re-evacua | paliativo (dinamico) |
| **Resolucao global mais fina (dx menor)** | centro esparso ganha ~2.25× particulas → rho/rho0 sobe; **reduz muito** (talvez nao 100% — ha gap de deslocamento) | candidato robusto; custo dt~0.6×, wall~3.4× |
| Suavizar pin (drag) | fragmenta | DESCARTADO (§22) |
| Subir threshold pin (0.8→0.9) | vacuo se realoca + âncora enfraquece | negativo (K.25a) |
| Pass N split na junção | over-pack adjacente ao pin → colapso dt | DESCARTADO (#27) |
| Pass N + timesteps individuais (Rota B) | preenche (= divisao celular) | correto, caro (reescreve integrador) |
| **C3 — inserir particulas na rede dx no vacuo** | preenche ao rho0; **dt intacto** (h=h0), **sem over-pack** (dx, nao split); nao move particula | candidato mais direcionado |
| Aceitar + KGC | nao preenche; operador ja valido | base cientifica OK |

**Decisao:** pin PRESERVADO (mexer nele = fragmentacao §22 ou realocacao K.25a). Testar **C3 primeiro** (mais barato, direcionado), depois finer-res se C3 falhar.

#### Rota C3 — Inserção de partículas no vácuo (IMPLEMENTADO 2026-06-02, aguardando validacao)

Preenche o vacuo INSERINDO particulas frescas na rede dx — escapa de TODAS as armadilhas anteriores:
- **NAO e split Vacondio** (1→7 colado, over-packa #27). Insere na rede dx → densidade local ≈ rho0, sem over-pack.
- **h=h0** (nao encolhe) → **dt INTACTO** (escapa #29, o killer do Pass N).
- **NAO move particula existente** → nao congela (#31), nao depende de separar vacuo↔frontier por rho_b (#32).
- Inseridas herdam `rho_b` alto da mãe-vácuo → pinadas → congelam e preenchem estavel.

**Implementacao** ([main.py](main.py) `post_step`, gated `use_insert`): a cada `INSERT_FREQ=200` iter, detecta `void_mask = rho_b>0.5 & rho/rho0<0.7` (zona estrutural, poupa tips <0.5); p/ cada particula-vacuo gera 6 candidatos hexagonais a distancia dx; insere onde o spot esta VAZIO (nenhum existente a <`0.7·dx`, via cKDTree + dedupe); campos herdados da mãe (`rho_b,cs,c_o,c_n,noise`), `u=v=0`, `m=dx²`, `h=h0`. `rho`/`sigma_a`/`L` recomputados no proximo passo. Cap `INSERT_MAX=100`/call. KGC mantido ON (ortogonal).

**Predicao C3:** vacuo central eliminado onde inserido (rho/rho0→~1); dt mantem ~0.021; `mass_total` sobe (controlado pelo cap — trade-off de inserir = adicionar "celulas"); morfologia dos braços intacta (tips <0.5 nao tocados). **Risco:** mecanismo ad hoc (inserir, nao Vacondio); `mass_total` crescente (monitorar); possivel descontinuidade de campo se interpolacao grosseira (herda só da mãe). Reporta no log via coluna `pass_n_spawned` (Pass N off).

**Criterios de sucesso C3:** (1) painel 3 SEM aglomerado de gaps no centro; (2) dt ~0.021; (3) `mass_total` controlada (<~110 em t=50s); (4) morfologia preservada (braços finos, sem blob); (5) `a_pressure`≤4 (sem over-pack das inseridas).

**RESULTADO C3 (2026-06-02, t=47s) — SUCESSO: PRIMEIRO mecanismo a preencher o vacuo central.**
- ✅✅ **Vacuo central PREENCHIDO** — painel 1: centro vira núcleo solido (nao mais anel oco); painel 3: centro vira aglomerado VERDE (`rho/rho0≈1`), os gaps azuis que persistiam em TODAS as rotas anteriores SUMIRAM. Primeira vez no projeto.
- ✅ dt=0.0215 (custo ZERO — h=h0, escapou #29), ✅ `a_pressure`=3.0 (sem over-pack — rede dx, escapou #27), ✅ morfologia preservada (~18-20 dendritos de núcleo solido), ✅ `a_marangoni`=14 (KGC mantido), `mean_v`=0.0009, `contrast_cs`=13 (≈T2g).
- ⚠️ **`mass_total` 101→104.5 (+3.4%/47s)** — ~2× o crescimento natural do BiomassGrowth (T2g +2.5%/47s). Inserções subiram 35→95/call (perto do cap 100). O vacuo regenera com a expansao → inserções crescem → massa sobe. Trade-off de inserir (= adicionar celulas), controlado mas com tendencia crescente.

**Diagnostico do caveat:** como o centro JA aparece preenchido (frames), nao e sub-preenchimento (cap-limite) — e **over-inserção** na franja de under-density LEVE (rho/rho0 0.6-0.7) que nao sao buracos reais de kernel. Inseri-las so inflava massa.

#### Rota C3.2 — `INSERT_RHO_TRIG: 0.7 → 0.6` (IMPLEMENTADO 2026-06-02, aguardando validacao)

Alavanca unica: fillar so vacuo PROFUNDO. Violeau §3.6: `σ_a<0.85` (~15% erro) ~ `rho/rho0<0.6` (regime severo); a franja 0.6-0.7 e under-density leve (tolerable, nao e buraco real). **Predicao:** inserções caem ~30-50%, massa ~+2%/50s (t=100s ~107 vs ~110+), vacuo profundo (<0.6) continua preenchido. **Trade-off:** pode sobrar anel TENUE 0.6-0.7 no painel 3.

**RESULTADO C3.2 t=50s — OK (massa +2.1%, inserções decaindo 25→0-6, vacuo profundo cheio).** Validacao curta PASSOU.

**RESULTADO C3.2 t=100s — RUNAWAY na 2a metade (lição #34):** a auto-limitacao do t=50 era TRANSIENTE. Inserções: t=0-50 decaem 25→15 (✓), mas **t=57-100 explodem 37→87→102 e CRAVAM no cap** (100-102/call de t~67 em diante). `mass` 101→**110.2 (+9%)** (+2.5% na 1a metade, +6.4% na 2a). `contrast_cs` colapsou 360→**8.2** (vs T2g ~12-15, e ainda caindo) — motor degradando. Frame t=98: painel 3 com distribuicao **assimetrica/bagunçada** (centro verde limpo do t=46 deu lugar a aglomerado lopsided). dt=0.0197, `a_pressure`=3.0 (numericamente estavel, mas motor+massa NAO).

**Causa-raiz do runaway:** partículas inseridas (`rho_b` 0.5-0.6) NAO eram filler inerte — elas (a) **crescem** (BiomassGrowth) → maturam p/ `rho_b≥0.8` → expandem o núcleo pinado → **nucleus maturation (lição #18)** re-disparada → núcleo maior → junção maior → MAIS inserção (feedback positivo); (b) **produzem cs** → inflam cs no interior → `contrast_cs` colapsa. **Inserir como biomassa ativa realimenta o problema que gera o vacuo.**

#### Rota C3.3 — FILLER INERTE (IMPLEMENTADO 2026-06-02, aguardando validacao)

Corrige o erro de modelagem do runaway: as inseridas viram **suporte de kernel inerte**, nao biomassa ativa. Nova propriedade `is_filler` (0=real, 1=inserida). Tres gates:
- **[src/equations.py](src/equations.py) `BiomassGrowth.loop`:** filler NAO cresce (gate `is_filler<0.5`) → quebra o feedback de nucleus maturation (#18).
- **[src/equations.py](src/equations.py) `SurfactantEquation.post_loop`:** filler NAO produz cs (`qs=0` se filler) → preserva `contrast_cs` (so difunde passivamente).
- **[src/scheme.py](src/scheme.py) `CustomEulerStep`:** filler e PINADO (`u=v=0`) → congela como suporte de densidade na junção (evita que filler em rho_b~0.5-0.6, dentro do gate flagelar, se mova e reabra o vacuo).
- **[main.py](main.py):** propriedade `is_filler` (init 0); inserção C3 marca `is_filler=1`.

**Predicao C3.3:** inserções MUITO menores (so vacuo geometrico da expansao, sem o feedback de maturacao) → massa cresce devagar (sem cravar o cap); `contrast_cs` preservado (~T2g, sem o colapso a 8); vacuo profundo preenchido por filler congelado; morfologia simetrica (sem o lopsided do feedback). Filler carrega `rho_b` so p/ densidade/pin (suporte de kernel), nao para o motor. **Resta crescimento de massa geometrico** (cada filler = dx²) — bem menor que o feedback, mas monitorar; se ainda crescer demais, A3 (finer-res global) e o fallback. **Validacao pendente:** `make run` t=100s.

**Criterios C3.3:** (1) inserções NAO cravam o cap em t>57s (sem runaway); (2) `mass_total` < ~106 em t=100s; (3) `contrast_cs` ≥ ~12 sustentado (motor preservado); (4) painel 3 centro preenchido E simetrico; (5) morfologia dendritica sustentada.

**RESULTADO C3.3 t=100s — SUCESSO nos eixos criticos (motor/morfologia/vacuo central):**
- ✅✅ **`contrast_cs` PRESERVADO em platô ~14** (até subiu 13.7→14.9 no fim) vs colapso a 8.2 da C3.2 — o filler-nao-produz-cs consertou a inflacao. Motor ≈ T2g. (criterio #3 PASS)
- ✅ **centro preenchido E SIMETRICO** (frame t=97) — sumiu o lopsided da C3.2 (maturacao desligada removeu a assimetria). (criterios #4,#5 PASS)
- ✅ dt=0.0211 (zero custo), `a_pressure`=2.99, `mean_v`=0.0009 (saudavel).
- ⚠️ **massa +6.7% (101→107.8)** e inserções ainda sobem na 2a metade (27→102, cap em t~84) — mas e crescimento GEOMETRICO **bounded pelo cap** (linear, NAO exponencial como o feedback da C3.2). Mass so +1.8 acima do T2g. (criterios #1,#2 marginais)

**Causa do residual geometrico:** o filler congelado forma uma CASCA frozen crescente; o rim real expande; o vacuo **migra para a borda externa da casca** → novas inserções na fronteira que cresce (mesma cinematica do pin, movida para fora). Linear-limitado, nao patologico. **C3.3 resolve o vacuo CENTRAL de forma sustentavel.**

#### Rota C3.4 — TRIGGER CANONICO `σ_a < 0.85` (IMPLEMENTADO 2026-06-02, aguardando validacao)

Substitui o proxy `rho/rho0 < 0.6` (C3.2) pela metrica DIRETA da partição da unidade `σ_a < 0.85` ([main.py](main.py) `void_mask = rho_b>0.5 & sigma_a<0.85`). Motivo (auditoria de rigor — §10): o trigger de deficit de kernel deve usar `σ_a` (Violeau §3.6, Liu §3.3.3 — ~15% erro), nao `rho/rho0` que e enviesado por massa (lição #28). `σ_a` ja era calculado por `KernelSum` (so usado no Pass N ate agora). O gate estrutural `rho_b>0.5` protege contra falso-positivo de free-surface (Liu §6.5 — rim diluido tem σ_a<1 por kernel truncado, nao por vacuo real).

**Valor:** `0.85` e o limiar canonico (vs o `0.6` em rho/rho0, que correspondia a σ_a≈0.6 ~40% erro — ESTRITO demais, escolhido so p/ conter massa). Em massa uniforme `ρ_a/ρ0 ≈ σ_a` (lição #28), entao C3.4 e MAIS PERMISSIVO que C3.2: pega vacuos MODERADOS (σ_a 0.6-0.85) que o 0.6 ignorava.

**Predicao:** mais particulas-vacuo detectadas → mais inserções → `mass_total` provavelmente cresce MAIS que +6.7% da C3.3. **Trade-off rigor↔massa explicito:** o controle de massa agora e responsabilidade dos parametros de ENGENHARIA (`INSERT_MAX`, `INSERT_FREQ`), NAO de deturpar o threshold fisico. Se a massa for inaceitavel, a discussao e sobre o cap, nao sobre afrouxar o critério de Violeau.

**RESULTADO C3.4 (2026-06-10, t=97s) — KEEP: melhoria de rigor sem custo, resultado ≈ C3.3.**
- ✅ dt=0.0203 (4800 iter/97s — custo ZERO), `a_pressure`=3.01 platô (sem over-pack), `contrast_cs` platô **~14** (NAO colapsa — = C3.3, vs 8 da C3.2), `a_marangoni` 12-16 (KGC on), `mean_v` 0.0007-0.0009 (tip-only vivo).
- ✅ **Vacuo CENTRAL preenchido** — painel 3 (frames 12/25): centro verde solido (`rho/rho0≈1`), como na C3.3.
- ⚠️ `mass_total` 101.1 → **108.7 (+7.5%)** — só **+0.8% acima** da C3.3 (+6.7%). Predicao confirmada (0.85 mais permissivo → mais massa) MAS ganho pequeno.
- **DESCOBERTA (cap binding):** `pass_n_spawned` crava **100-102/call de t≈54s em diante** (cap `INSERT_MAX=100`). Ou seja, `σ_a<0.85` detecta MAIS vacuo estrutural do que o cap permite preencher → **a escolha do threshold (proxy 0.6 vs canonico 0.85) quase NAO afeta o regime permanente; quem controla a massa e `INSERT_MAX`/`INSERT_FREQ`**. Isso **valida empiricamente** o claim C3.4: o controle de massa e engenharia, nao deturpar o critério fisico de Violeau. Migrar para o trigger canonico foi de graca (mesmo resultado validado da C3.3) e mais defensavel em revisao ("déficit de partição da unidade σ_a<0.85, Violeau §3.6").
- Vacuo da **frontier ativa** (anel + dendritos, `rho_b<0.5`) intocado — gate `INSERT_RHO_B_MIN=0.5` exclui de proposito (escopo do #2, lição #35).

**Constantes de engenharia (declaradas explicitamente, NAO derivadas):** `INSERT_FREQ=200` (cadencia), `INSERT_MAX=100` (cap de massa), `INSERT_PROX=0.7·dx` (anti-over-pack, conceito Liu §6.4 mas valor empirico), `INSERT_RHO_B_MIN=0.5` (gate estrutural↔frontier do modelo, lição #35 — posicionado acima do pico flagelar 0.4). So `σ_a<0.85` e literatura-derivado.

**VACUO REMANESCENTE (diagnosticado pelo usuario, 2026-06-02):** o vacuo central foi resolvido, mas PERSISTE (a) num **anel intermediario** entre o centro cheio e as bases dos dendritos, e (b) **dentro dos dendritos**. Ambos estao na **FRONTIER ATIVA (`rho_b<0.5`)** — a zona motil (swarmers/tips) que C3.3 EXCLUI de proposito (`INSERT_RHO_B_MIN=0.5`). Distincao fisica fundamental:
- **Vacuo estrutural (`rho_b>0.5`, junto ao pin):** frozen, nao-motor → filler inerte congelado e correto (C3.3 ✅).
- **Vacuo da frontier ativa (`rho_b<0.5`, anel+arms):** movel, crescente, E O MOTOR → filler congelado o MATARIA (#31, shifting v1); biomassa ativa realimenta (#34). NAO pode ser preenchido por C3-filler.
- A preocupacao CIENTIFICA da frontier (∇cs espurio nos arms) JA e mitigada pela **KGC** (on) — a amplificacao `a_marangoni` 8→14 e a KGC corrigindo o gradiente nos arms sub-resolvidos. Operador valido apesar do vacuo.
- Preencher FISICAMENTE o vacuo dos arms exige refinamento ATIVO que se move com o arm: Pass N (velocidade herdada) + Rota B (timesteps individuais [T9]) para o custo de dt; OU inserção ATIVA (C4: is_filler=0 + u,v herdados do arm) a h=h0; OU A3 (finer-res global) reduz todos os vacuos pelo orcamento. Ver §9 lição #35.

#### Roadmap de decisao — vacuo da FRONTIER ATIVA (anel + dendritos) — A DECIDIR (2026-06-02)

4 opcoes documentadas para decisao posterior. Todas pressupoem C3.3 mantido (vacuo central) + KGC mantido. Ordem por custo/risco crescente:

**Opcao 1 — ACEITAR (KGC ja resolve a preocupacao cientifica). Custo ZERO, ja implementado.**
A KGC (on) corrige o `∇cs` da Marangoni nos arms sub-resolvidos (a amplificacao `a_marangoni` 8→14 e exatamente isto — Bonet-Lok, Liu §3.3). Se o objetivo e "o gradiente nos braços e valido" (a base da tese), JA esta feito. O vacuo visual no painel 3 dos arms seria entao **cosmetico**, nao cientifico. **Pre-passo recomendado antes de qualquer implementacao:** instrumentar `σ_a` (ja computado por KernelSum) e/ou o erro do gradiente nos arms no log/plot, para MEDIR se o vacuo da frontier degrada a fisica ou se a KGC ja o neutralizou. Decide empiricamente se 2/3/4 sao necessarios.

**Opcao 2 — Inserção ATIVA nos arms (C4). Custo dt ZERO (h=h0), risco medio.**
Estender a inserção C3 a frontier (`rho_b<0.5`) MAS com particulas ATIVAS e velocidade-herdada: `is_filler=0`, `u,v` herdados do arm-pai (movem-se COM o arm → escapam #31), `rho_b/cs/c_n` herdados, `m=dx²`, `h=h0` (dt intacto). Crescem/produzem cs como funcao natural do arm — **sem feedback #18** porque nos arms NAO ha nucleo pinado para expandir (arm crescer e desejado). E "divisao celular no dendrito" barata. **Risco:** crescimento de massa/cs — precisa testar se NAO infla cs como C3.2 (a diferenca: nos arms a producao de cs e correta/local, nao afogamento do interior). Monitorar `contrast_cs` e `mass`. Se inflar → recuar para opcao 4.

**Opcao 3 — Pass N (split) nos arms + Rota B (timesteps individuais [T9]). Canonico, custo ALTO (reescrita do integrador).**
Refinamento Vacondio herda velocidade (move com arm) e CONSERVA massa (divide, nao adiciona → sem crescimento de massa, melhor que C3/C4 nesse aspecto). O unico bloqueio era o dt (#29, h pequeno) — resolvido por timesteps individuais/blocos hierarquicos (Springel 2005, limiter Saitoh-Makino [T9]). Exige reescrever `CustomEulerStep` + pin para integracao multi-dt. **Fix "fisicamente correto" de longo prazo** (= divisao celular real com resolucao verdadeira sub-dx), mas o mais caro de implementar.

**Opcao 4 — A3 (resolucao global mais fina). Robusto, custo dt~0.6× + wall~3.4×.**
`dx 0.054→0.036` (ou menor) → ~2.25× particulas globais → centro, anel E arms todos mais densos pelo ORCAMENTO de particulas. Sem over-pack, sem mexer no pin, sem congelar, sem feedback. Reduz TODOS os vacuos de uma vez. Incerteza: parte do vacuo e gap de deslocamento → pode reduzir sem zerar. **Fallback robusto** se 2 inflar e 3 for caro demais.

**Recomendacao:** comecar pela **Opcao 1** (instrumentar `σ_a`/erro de gradiente nos arms p/ medir se o vacuo da frontier importa cientificamente ou se KGC ja resolveu). So implementar 2/3/4 se a medicao mostrar degradacao real. Se precisar preencher: **2** (mais barata, dt intacto) → se inflar cs, **4** (robusta) → **3** como fix definitivo.

#### #2 MEDIDO — instrumentacao σ_a nos braços (2026-06-10) → OPCAO 1 JUSTIFICADA

Implementada a Opcao 1 (medir antes de agir): 3 colunas no log ([main.py](main.py) — `min_sig_arms`, `mean_sig_arms`, `frac_lowsig_arms` na frontier `rho_b∈[0.1,0.5]`, lição #35) + 4º painel `σ_a` no [plot.py](plot.py) (RdYlGn, braços<0.85 contornados). `sigma_a` adicionado a `add_output_arrays`. `KernelSum` ja computava σ_a todo passo (independente do Pass N).

**Resultado (run t=99s, HDF5 final):** o déficit de σ_a na frontier NAO esta nas pontas — esta no **anel interno da junção**. Perfil radial dos braços:

| Anel radial | n | mean σ_a | frac<0.85 |
|---|---|:---:|:---:|
| r<0.8 (junção interna) | 38 | 0.936 | **0.18** |
| r[0.8,1.5) | 17 | 0.983 | 0.12 |
| r[1.5,2.5) | 8 | 0.901 | 0.12 |
| **PONTAS r≥2.33 (outer 25%)** | 21 | **0.977** | **0.05** |

mean σ_a braços (agregado) = 0.95; 11 partículas (13%) com σ_a<0.85, concentradas em r≈0.93.

**Conclusao — vácuo da frontier NAO e problema científico:**
1. **Pontas (onde Marangoni dirige Mullins-Sekerka, base da comparacao Trinschek) bem resolvidas:** σ_a=0.977, ~2% erro, só 5% em déficit. O ∇cs que dirige o fingering NAO e espurio ali.
2. **O déficit real e o anel interno** (r<0.8, 18% em déficit) — ATRAS da frente ativa, cinematico (pin congela núcleo + rim expande, #32/#35), NAO dirige morfologia.
3. **KGC ja cobre o residual** — corrige ∇cs onde σ_a<1 (Liu §3.3); a amplificacao a_mar 8→14 e a prova. `contrast_cs` platô ~14, `a_pressure` 3.0.

**Veredito: Opcao 1 (aceitar+KGC) confirmada pelos dados.** Nao precisa C4/PassN/finer-res para validade científica. Limpar o anel interno (r<0.8, rho_b∈[0.1,0.5]) seria decisao ESTETICA (C3 estendido com cuidados #34/#35), nao científica.

---

### Pass L — Superficies rugosas (Objetivo 1) — **BLOQUEADO ate Pass M**

> ⚠️ **PRE-REQUISITO ATUALIZADO (2026-05-01):** Pass L agora requer **Pass M (substrato consumivel) implementado e validado primeiro**. So apos a simulacao reproduzir a morfologia de [reference.jpg](reference.jpg) com **dendritos finos AR≥1:5 + agar limpo entre eles** (Pass M criterio), iniciar Pass L. Ver §1 Objetivo 1, §2.2 e §10 Proibicoes.
>
> **Justificativa atualizada:** Rugosidade e uma extensao fisica do modelo (objetivo de tese final), nao uma muleta para compensar mecanismos faltantes. K.27 demonstrou que mecanismos hidrodinamicos puros (Marangoni + Flagelar + EOS) **nao produzem morfologia 100% conforme reference.jpg sem o substrato consumivel `c_n`**. A baias precisam ser suprimidas ATIVAMENTE para que o swarm reproduza o padrao natural de P. aeruginosa. Adicionar rugosidade antes de Pass M:
> - Oculta o bug real (halo radial das baias seria mascarado por canalizacao topografica espuria).
> - Impossibilita isolar a contribuicao da rugosidade vs. supressao quimica das baias.
> - Compromete a validade cientifica da tese — comparacao "swarm liso vs swarm rugoso" requer baseline liso CORRETO primeiro.

**Objetivo (quando desbloqueado pos-Pass M):** Substituir paredes planas em `particles.py` por topografia irregular para estudar o efeito da rugosidade no padrao de swarming ja estabelecido.
- **Geometria:** Senoidais com amplitude `A` e comprimento de onda `Lambda` controlados, ou perfil aleatorio com espectro de potencia definido.
- **Interacao:** Particulas solidas com no-slip (velocidade zero imposta via `MomentumEquation` `sources=["solid"]`).
- **Validacao:** Colonia ja dendritica deve apresentar modulacao da morfologia pela rugosidade (canalizacao pelos vales, ancoragem de dedos, etc.), nao formar dendritos pela primeira vez.