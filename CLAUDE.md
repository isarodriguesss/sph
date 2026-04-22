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
- **Referencia visual obrigatoria:** [reference.jpg](reference.jpg) (Michiels et al.) — colonia com ~15-20 dendritos radiais longos, `AR >= 1:5`, separados por agar limpo, nucleo central coeso com picos de surfactante nas pontas (painel C da referencia).
- **Criterio de sucesso (Pass I-K):** aproximacao monotonica da morfologia de `reference.jpg` — dendritos cada vez mais finos, longos e separados.
- **Criterio de falha:** expansao circular uniforme, anel oco, colapso para disco compacto, ou blob ameboide com muitos bumps curtos (AR ~1:2).
- Nao declarar um Pass como "bem-sucedido" baseando-se apenas em metricas escalares — a forma da colonia **comparada a `reference.jpg`** e o validador final.
- **Gatilho para Pass L (rugosidade):** so propor apos frames mostrarem dendritos de AR >= 1:5 separados por agar limpo. Ate la, refinar mecanismos hidrodinamicos (Pass I-K).

### 2.3 Predicao antes da acao

- Para cada mudanca de parametro, a IA deve apresentar uma **predicao quantitativa** do efeito esperado (ex: "reduzir `D` de 4e-3 para 2e-3 deve aumentar `|nabla cs|` em ~40% e reduzir `L_D` de 1.5h para 1.1h").
- Apos a simulacao, comparar predicao vs. resultado. Discrepancias > 2x devem ser investigadas antes do proximo passo.

### 2.4 CRITERIO BLOQUEANTE — Core Pinning (diagnosticado 2026-04-22, pos-K.15)

**Contexto biologico:** na literatura (Michiels et al., Tremblay et al., Kearns — P. aeruginosa PA14 em agar swarming), a colonia exibe **nucleo denso imovel** (matriz EPS madura) com **dendritos finos sendo a unica frente de expansao**. A velocidade radial da colonia = velocidade das pontas. Interior e baias entre dendritos sao **estaticos**.

**Falha observada na simulacao (K.15, t=0→33.8s):** apesar de motor sustentado e dendritos visiveis (AR 1:5-1:7), a **colonia inteira expande radialmente** — toda a interface (borda, baias, pontas) avança junto. Dendritos sao "bumps arrastados" na superficie de um blob inflando, nao extensoes de um nucleo ancorado. Consequencia morfologica: ausencia de tip-splitting e de estrutura fractal — pontas nao tem trajetoria propria para bifurcar.

**Metricas diagnosticas (core pinning failure):**
- **`a_pressure` (orcamento §8: <3) observada 15-30 rotina, picos >50** (~5-17× excesso). EOS com `edge_fade=1.0` no interior (`rho_b>0.5`) pressuriza o nucleo quando `rho > rho0`, empurrando radialmente.
- **`n_fast` 900-1300** — ativacao em massa; fracao ativa/total > 50%. Nao localizada nas pontas.
- **Razao drag core/borda = 2.2×** (`gamma_mature=1.5·gamma`) — insuficiente (§3.2 sugere 5-10×).

**Criterios de aceitacao para desbloquear qualquer avanco:**
1. `a_pressure` sustentada **<5** em t>20s (max eventual <15).
2. `n_fast / N_interface < 2` — ativacao predominantemente nas pontas, nao no bulk.
3. Visualmente nos frames: **baias entre dendritos devem permanecer estacionarias** entre frames sucessivos (t-spaced ≥3s); so as pontas avançam radialmente.
4. Razao efetiva `v_core / v_tip < 0.2` — nucleo pinnado relativamente as pontas.

**BLOQUEIO DURO:** enquanto core pinning failure persistir, sao **proibidas**:
- Novas mudancas em parametros de surfactante (`sigma`, `beta`, `lambda`, `D`, `k_consume`).
- Propostas de Etapa 2 (tip-splitting, Peclet-Mullins via β).
- Propostas de nova fisica (substrato consumivel, Ising, etc.).
- **Pass L (rugosidade) permanece BLOQUEADO** (segundo bloqueio ativo, alem de §1/§2.2).

**Unicas acoes permitidas:** tunar `gamma_mature` (core drag), `c0` (EOS stiffness), e `edge_fade` geometria/domínio da BiomassEOS — parametros que atuam diretamente no pinning.

---

## 3. Skill: SPH Bacterial Swarming Specialist

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

**Frente 4 — Pressao Osmotica (Equacao de van't Hoff)** *(a implementar)*
- `Pi = iCRT` — pressao proporcional a concentracao de solutos intra/extracelulares.
- `OsmoticForce` existe em[equations.py:318-363](src/equations.py#L318-L363) mas esta desabilitada. A EOS coesiva atual (ramo atrativo) funciona como substituto simplificado.
- Meta: reativar com formulacao termodinamica quando a estabilidade numerica permitir.

**Frente 5 — Motilidade Flagelar orientada por gradiente** *(✅ implementado — Pass K)*
- `F_flag = -f0 * gate(rho_b) * n̂(∇cs)` — magnitude **constante** (`f0 * gate`), direcao oposta a `∇cs`.
- **Escolha-chave:** magnitude nao depende de `|∇cs|` — por isso sobrevive a saturacao do reservatorio (quando `|∇cs|` cai mas direcao persiste). Diferente da Marangoni, que colapsa quando `cs` uniformiza.
- `motility_gate`: smoothstep em `rho_b ∈ [0.1, 0.6]` com pico em `rho_b = 0.35` — so swarmers de borda; nucleo maduro (`rho_b ≈ 1`) e agar livre (`rho_b ≈ 0`) sao imunes.
- Implementacao usa padrao loop/post_loop (como `BiomassGradient`): `loop` acumula `grad_cs_x/y` pelo gradiente SPH simetrico, `post_loop` normaliza o vetor total e aplica a forca com sinal `-f0` (mesma direcao da Marangoni — para agar fresco).
- `f0 = 2.0` atual; orcamento teorico `f0 ~ 0.5*gamma*v_term = 3`.
- **Invariante critico:** `FlagellarForce` precisa executar **apos** `SurfactantEquation` para usar `cs` atualizado. Ver ordem em [scheme.py](src/scheme.py).

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
| `BiomassGrowth` | Crescimento logistico `drho_b/dt = r * rho_b * (1 - rho_b/rho_max)` | [equations.py:4-18](src/equations.py#L4-L18) |
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
| Viscosidade | `mu` | 0.012 | Pass I.2: reduzida de 0.025 para filamentos finos |
| Arrasto (base) | `gamma` | 60.0 | `a_drag` = `gamma*v_term` = 6 em v=0.1 |
| Arrasto (nucleo) | `gamma_mature` (scheme) | `1.5*gamma` | Pass I.4: razao core/edge 2.5x (era 0.3*gamma) |
| Coef. Marangoni | `beta` | **1.0** | Pass K.5: 4→1 (4x reduzido para domar feedback amplificado por K.6/K.7) |
| Producao surfactante | `sigma` | **0.8** | Pass K.9: 0.6→0.8 (ponto medio entre explosao K.7 e one-tip K.8) |
| Forca flagelar | `f0` | **0.5** | Pass K.5: 2.0→0.5 (4x reduzido junto com beta) |
| Difusao (biofilme) | `D` (`D_int`) | 1.5e-3 | Pass I.3: gradiente afiado na interface (`L_D_int`=0.93h) |
| Difusao (agar) | `D_ext` | 0.01 | Pass I.7, **validado K.10**: reduzir para 0.005 mata o campo long-range de Mullins-Sekerka (pior morfologia) |
| Decaimento | `lambda_` | 0.15 | Manter — confina `cs` mas permite penetracao de `L_D_ext` |
| Taxa crescimento | `r_growth` | 0.4 | Pass K: estende janela de swarmers de borda |
| Modelo producao | `sigma*qs*(1.2-rho_b)*noise*tip_boost*motile_boost` | — | K.6 adicionou `tip_boost = 1+3*grad_rho_b_mag`; K.7 adicionou `motile_boost = 1+50*min(|v|,0.1)` |
| Visc. artificial Monaghan | `alpha_mon` | 0.06 | Pass I.2: reduzida de 0.15 para gradientes afiados |
| Vel. som (EOS) | `c0` | 0.8 | B ~ 0.09; tensao `tension_ratio=0.02` (Pass I.5) |
| Smoothing kernel | `h_factor` | 1.8*dx | ~35 vizinhos por particula |
| Timestep | `dt` | 5e-5 | Adaptivo, CFL=0.4 |
| Grade | `x_dim, y_dim` | 150x150 | Pass I.8: resolucao aumentada |
| Dominio | `x/y_min/max` | [-3, 3]^2 | 6x6 centrado na origem |
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

**Proximos diagnosticos obrigatorios:**
1. **Pass K.14 (Etapa 1 — aplicado 2026-04-22):** sumidouro `-k_consume*rho_b*cs` em [src/equations.py:147](src/equations.py#L147), `k_consume=0.5`. Meta: estabilizar `mean_cs ∈ [0.25, 0.35]` e `contrast_cs ≥ 80` sustentado em t=60-100s. Executar `make run` e comparar log + frames vs K.13.
2. **Se K.14 estabilizar motor:** iniciar Etapa 2 (tip-splitting). Testar em ordem: (a) β:1.0→1.5 capilaridade, (b) ruido estocastico 0.01→0.05 na noise de producao, (c) se insuficiente, arquitetura nova — campo substrato `c_n` consumivel por ρ_b + chemotaxis positivo para `+∇c_n` na `FlagellarForce`.
3. **Se K.14 nao estabilizar:** k_consume insuficiente — testar 0.8 ou reformular lambda_int mais alto diretamente.
4. **Pass L (rugosidade) permanece BLOQUEADO** ate frames mostrarem morfologia fractal com tip-splitting matching `reference.jpg` de P. aeruginosa (Michiels et al.), nao apenas dendritos radiais simples.

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
- **Nao sugerir nenhuma mudanca em parametros de surfactante ou nova fisica enquanto core pinning failure (§2.4) persistir.** Criterios de aceitacao (todos): `a_pressure < 5` sustentado, `n_fast / N_interface < 2`, baias estacionarias entre frames t-spaced ≥3s, `v_core/v_tip < 0.2`. Unicas alavancas permitidas ate passar: `gamma_mature`, `c0`, `edge_fade` geometria. Motor sustentado sem core pinning nao produz tip-splitting — dendritos sao bumps arrastados (K.15).

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

### Pass L — Superficies rugosas (Objetivo 1) — **BLOQUEADO**

> ⚠️ **PRE-REQUISITO:** Pass L so pode ser iniciado apos a simulacao reproduzir a morfologia de [reference.jpg](reference.jpg) (dendritos de `AR >= 1:5`, separados por agar limpo, nucleo coeso) via mecanismos hidrodinamicos puros. Ver §1 Objetivo 1, §2.2 e §10 Proibicoes.
>
> **Justificativa:** Rugosidade e uma extensao fisica do modelo (objetivo de tese), nao uma muleta para compensar motor insuficiente. Se Marangoni + Flagelar + EOS nao produzem dendritos finos em paredes planas, a causa-raiz esta em um desses mecanismos. Adicionar rugosidade antes de resolver o problema hidrodinamico:
> - Oculta o bug real (dendritos "emergem" mas por razao errada).
> - Impossibilita isolar a contribuicao da rugosidade no mecanismo de fingering.
> - Compromete a validade cientifica da tese.

**Objetivo (quando desbloqueado):** Substituir paredes planas em `particles.py` por topografia irregular para estudar o efeito da rugosidade no padrao de swarming ja estabelecido.
- **Geometria:** Senoidais com amplitude `A` e comprimento de onda `Lambda` controlados, ou perfil aleatorio com espectro de potencia definido.
- **Interacao:** Particulas solidas com no-slip (velocidade zero imposta via `MomentumEquation` `sources=["solid"]`).
- **Validacao:** Colonia ja dendritica deve apresentar modulacao da morfologia pela rugosidade (canalizacao pelos vales, ancoragem de dedos, etc.), nao formar dendritos pela primeira vez.