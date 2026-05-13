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

### 3.3 Meta de superficies rugosas

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
| Coef. Marangoni | `beta` | **1.0** | Pass K.5: 4→1 (4x reduzido para domar feedback amplificado por K.6/K.7) |
| Producao surfactante | `sigma` | **1.5** | M-B.2 atual: K.15 0.8→1.2, M-B.2 1.2→1.5; K.20 k_consume controla saturacao |
| Forca flagelar | `f0` | **3.0** | K.22: 0.5→3.0 — alvo teorico §8 (f0~γ·v_term/2=3); K.21 diagnosticou regime subcritico |
| Difusao (biofilme) | `D` (`D_int`) | 1.5e-3 | Pass I.3: gradiente afiado na interface (`L_D_int`=0.93h) |
| Difusao (agar) | `D_ext` | **0.08** | K.18 0.01→0.04; M-B.2 0.04→0.08; `L_D_ext=0.327`; habilita focalizacao Mullins-Sekerka |
| Decaimento | `lambda_` | 0.15 | Manter — confina `cs` mas permite penetracao de `L_D_ext` |
| Sumidouro enzimatico | `k_consume` | **1.0** | K.20 0.5→2.0; revisado M-B para 1.0; cs_∞(interior)≈0.10 |
| Taxa crescimento | `r_growth` | **0.02** | M-B.7: K.21 0.15 → 0.02 — elimina tip pumping (bloqueio E) |
| Modelo producao | `sigma*qs*(1.2-rho_b)*noise*tip_boost*motile_boost*c_n_factor` | — | K.6 adicionou `tip_boost`, K.7 adicionou `motile_boost`, M-B.2 adicionou `c_n_factor` |
| Visc. artificial Monaghan | `alpha_mon` | **0.12** | K.23: I.2 revert parcial — previne instabilidade de tração SPH (I.2=0.06, pre-I.2=0.15) |
| Vel. som (EOS) | `c0` | 0.8 | B ~ 0.09; tensao `tension_ratio=0.08` (M-B.9b) |
| Tensao superficial EOS | `tension_ratio` | **0.08** | M-B.9b: I.5 0.02 → 0.08 (4×) — `B_tension`=0.00275; lição §23 |
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