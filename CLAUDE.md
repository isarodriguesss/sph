# CLAUDE.md

Este arquivo governa todas as interações do Claude Code com este repositório. As regras aqui descritas sao **obrigatorias** e tem precedencia sobre qualquer comportamento padrao.

---

## 1. Projeto e Referencial Teorico

**Tese de Doutorado:** *"Simulacao do Swarming de Pseudomonas aeruginosa em superficies rugosas utilizando Hidrodinamica de Particulas Suavizadas (SPH)"* (Rodrigues, 2025).

Simulacao 2D em SPH (via [PySPH](https://pysph.readthedocs.io/)) da dinamica de biomassa bacteriana: crescimento logistico, forcas de Marangoni induzidas por gradientes de surfactante (ramnolipideo), fluxo viscoso, arrasto diferenciado e coesao via EOS. O objetivo final e reproduzir o padrao de *swarming* dendritico de *P. aeruginosa* (fingering de Mullins-Sekerka) observado experimentalmente (Michiels et al.).

### Objetivos de longo prazo (guiam decisoes arquiteturais)

1. **Superficies rugosas** *(em aberto — prioridade apos Pass K)* — Partículas de contorno estaticas com geometria irregular que interagem mecanicamente com o fluido e alteram os campos de difusao/escoamento. Toda decisao de refatoracao deve preservar a capacidade de substituir as paredes planas atuais por topografias arbitrarias.
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
- **Criterio de sucesso:** formacao de dendritos/gavinhas com perimetro irregular (fingering). **Criterio de falha:** expansao circular uniforme, anel oco, ou colapso para disco compacto.
- Nao declarar um Pass como "bem-sucedido" baseando-se apenas em metricas escalares — a forma da colonia e o validador final.

### 2.3 Predicao antes da acao

- Para cada mudanca de parametro, a IA deve apresentar uma **predicao quantitativa** do efeito esperado (ex: "reduzir `D` de 4e-3 para 2e-3 deve aumentar `|nabla cs|` em ~40% e reduzir `L_D` de 1.5h para 1.1h").
- Apos a simulacao, comparar predicao vs. resultado. Discrepancias > 2x devem ser investigadas antes do proximo passo.

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
| Coef. Marangoni | `beta` | 4.0 | Estavel — `a_mar` pico 148 com sigma=2.0 |
| Producao surfactante | `sigma` | 2.0 | Pass I.7: +67% para compensar drenagem por `D_ext` |
| Difusao (biofilme) | `D` (`D_int`) | 1.5e-3 | Pass I.3: gradiente afiado na interface (`L_D_int`=0.93h) |
| Difusao (agar) | `D_ext` | 0.01 | Pass I.7: campo de longo alcance (`L_D_ext`=2.4h) |
| Decaimento | `lambda_` | 0.15 | Manter — confina `cs` mas permite penetracao de `L_D_ext` |
| Taxa crescimento | `r_growth` | 0.8 | Pass G: crescimento lento estende tau_sat 0.8s -> 2.5s |
| Modelo producao | `qs * (1.2 - rho_b)` | `[rho_b^2/(rho_b^2+0.01)] * (1.2 - rho_b)` | Pass H aplicado: piso de 17% no nucleo maduro |
| Visc. artificial Monaghan | `alpha_mon` | 0.06 | Pass I.2: reduzida de 0.15 para gradientes afiados |
| Vel. som (EOS) | `c0` | 0.8 | B ~ 0.09; tensao `tension_ratio=0.02` (Pass I.5) |
| Smoothing kernel | `h_factor` | 1.8*dx | ~35 vizinhos por particula |
| Timestep | `dt` | 5e-5 | Adaptivo, CFL=0.4 |
| Grade | `x_dim, y_dim` | 100x100 | |
| Dominio | `x/y_min/max` | [-3, 3]^2 | 6x6 centrado na origem |
| Perturbacao inicial (rho_b) | `azimuthal_perturb` | `0.5*cos(8*theta)` | Pass I.1: N=8 unificado (era N=16) |
| Perturbacao inicial (noise) | `noise` | `1.0 + 0.4*sin(8*theta) + 0.03*rand()` | Pass I.1: N=8, ruido reduzido |

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

- **Pass K — Motilidade flagelar implementada; morfologia ainda lobular:** Com Pass J+K, motor sustentado ate t=24.5s confirmado (`log.csv`). `a_mar` pico 182, `mean_v` estavel 0.0075. Porem, frames 030/047 mostram colonia lobular com ~15-20 bumps curtos (AR ~1:2), nao dendritos finos (meta: 8-12 dedos com AR >= 1:5). No painel `cs` o anel da borda e quase uniforme — falta seleção competitiva entre pontas. Hipoteses:
  - `f0 = 2.0` pode ser insuficiente — `FlagellarForce` redundante com Marangoni.
  - Seleção competitiva requer canalização geométrica (Pass L, superfícies rugosas).
  - Gate de Marangoni (`grad_rho_b_mag ∈ [0.05, 0.6]`) pode ser largo demais — ativa em quase todo o perímetro.

**Proximos diagnosticos obrigatorios:**
1. Adicionar `a_flag` no log e `au_flag` no `plot.py` (painel extra) para validar que `FlagellarForce` esta ativa na borda.
2. Rodar ate t=100s para confirmar motor sustentado alem da janela atual.
3. Se `a_flag ≈ 2` e morfologia ainda blob → aumentar `f0: 2.0 → 3.5` ou pular para Pass L.

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
- **Nao** declarar sucesso de um Pass sem evidencia visual da morfologia.
- **Nao** introduzir dependencias externas sem justificativa.
- **Nao** refatorar a ordem das equacoes em `scheme.py` sem verificar invariantes.
- **Nao** remover equacoes desabilitadas (como `OsmoticForce`) que sao parte do roadmap.

---

## 11. Roadmap: Passes I-L (Correcao Morfologica e Extensoes Fisicas)

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

**Resultado (`log.csv` ate t=24.5s, 2026-04-17):** Motor sustentado confirmado — `mean_v` estavel 0.0075, `n_fast` 150-300, `max_cs` cresce ate 8.2, `a_mar` pico 182 (+23% vs I.7). Morfologia: colonia lobular com ~15-20 bumps curtos (AR ~1:2) — transicao para dendritos finos nao confirmada. Gap: seleção competitiva de dedos longos continua aberta.

**Proximos diagnosticos pendentes:**
1. Instrumentacao de log/plot — adicionar `a_flag` ao `LOG_HEADER`/writer em `main.py` (variavel `fluid.au_flag` ja existe); adicionar painel `au_flag` em `plot.py`.
2. Rodar ate t=100s para confirmar motor indefinido.
3. Se morfologia continuar lobular apos (1)+(2): testar `f0: 2.0 → 3.5` ou avancar para Pass L.

### Pass L — Superficies rugosas (Objetivo 1)

**Objetivo:** Substituir paredes planas em `particles.py` por topografia irregular.
- **Geometria:** Senoidais com amplitude `A` e comprimento de onda `Lambda` controlados, ou perfil aleatorio com espectro de potencia definido.
- **Interacao:** Particulas solidas com no-slip (velocidade zero imposta via `MomentumEquation` `sources=["solid"]`).
- **Validacao:** Colonia deve canalizar pelos vales da rugosidade.
- **Hipotese:** A canalizacao geometrica pelas paredes pode ser o mecanismo faltante para selecionar 8-12 dedos em vez de ~15-20 bumps — complementa os mecanismos hidrodinamicos (Marangoni+Flagelar) com restricao topologica.