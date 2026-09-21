# CRITÉRIOS PRÉ-REGISTRADOS — série de rugosidade R1-R4 (Pass L)

> Escrito em **2026-09-18, ANTES de qualquer rodada da série**. §2.3 do CLAUDE.md.
>
> Fica em `docs/` e não em `runs/<R>/CRITERIOS.md` porque `tools/archive_run.sh` faz `rm -rf`
> no destino e já apagou o CRITERIOS.md pré-registrado do P2R22. Copiar para cada run
> **depois** de arquivar.

## 1. A pergunta

**A rugosidade atrapalha o swarm?** É uma comparação, não um mapa de fases. Por isso o padrão é
fixo (pilares), a alavanca é a **dose** (`φ`, via `Λ`, com `A = 10 dx` constante), e a evidência é
a **monotonicidade ao longo de R1→R3**, não a diferença liso × rugoso isolada.

Com uma realização por dose, **uma resposta não-monotônica não autoriza conclusão** — exigiria
réplicas. Isso está pré-registrado para não ser decidido depois de ver o resultado.

## 2. Referência — P2R23, medido em 2026-09-18

| grandeza | valor |
|---|---:|
| `R99`(t=50) | **4.067** |
| **`dR/dt` em t ∈ [25,50]** | **0.0756** |
| largura EDT em 0.6R / ciclo | 5.4 dx / 0.61 |
| tortuosidade | 1.00-1.02 |
| dedos / AR | 15 / 9.9 |
| `a_mar_bio_p95` (mediana t>30) | 0.73 |
| `void_15` | 0.02% |
| `a_press_med` / picos > 4 | 3.22 / 31% |
| **`frac(c_s > 0.45)` nas vivas** | **0.007** |
| `frac(c_s > 0.45)` no corpo | 0.164 |
| **`max c_s` no ágar a < 2h do corpo** | **0.300** |
| frac. desse ágar acima de `COL_CS_MIN` = 0.3 | **0.0000** |

## 3. Predição (eixo 5 da análise preditiva)

Dois canais de estorvo, com **limiar de dose diferente** — é isso que torna a predição
falsificável, e não só "a rugosidade atrapalha".

> **CORRIGIDA em 2026-09-18 pelo pré-passo (seção 8).** Eu havia previsto R1 "transparente"
> raciocinando pela fresta estática (18 dx = 3.3 larguras de braço). Medido: **100% dos braços
> interceptam pelo menos um pilar nas TRÊS doses**, com 1.50 / 2.38 / 4.02 contatos por braço.
> A fresta diz que o braço *pode* passar entre dois pilares, não que passe — ao percorrer 53 dx
> radiais numa rede 2D ele encontra pilares necessariamente. O zero da escada de dose é o
> **run liso**, não R1.

| | `Λ` | **contatos/braço** | `φ` | **τ** | **`dR/dt`** | dedos |
|---|---:|---:|---:|---|---|---|
| liso | — | **0** | 0 | 1.00-1.02 | 0.0756 | 15 |
| **R1** | 28 dx | **1.50** | 11.6% | 1.05-1.15 | 0.064-0.072 (−5 a −15%) | 12-15 |
| **R2** | 20 dx | **2.38** | 22.7% | 1.15-1.30 | 0.053-0.064 (−15 a −30%) | 9-13 |
| **R3** | 16 dx | **4.02** | 35.4% | > 1.30 | 0.038-0.053 (−30 a −50%) | 6-10 |

As faixas de `dR/dt` e de dedos são estimativas ancoradas na contagem de contatos, não medições;
o que está de fato pré-registrado é a **ordenação monotônica** e o desacoplamento dos dois canais
abaixo. **O resultado científico é a distribuição de desfechos dos contatos** (desvio / divisão /
parada), não um alvo numérico — qualquer das três respostas é informativa.

1. **Mecânico** — o braço é desviado antes de ser barrado: **τ responde em R2 com `dR/dt` ainda
   dentro do ruído; `dR/dt` só cai em R3.** Se os dois caírem juntos já em R1, a predição erra e
   a causa provável é atrito (ver R4), não bloqueio.
2. **Químico** — os pilares são fluxo zero para `c_s`; ele acumula a montante e o gradiente à
   frente do líder enfraquece. Apareceria como queda de `a_mar_bio_p95` **já em R1-R2**, antes de
   qualquer efeito mecânico. Se aparecer, é resultado próprio: nenhum trabalho levantado na
   revisão acopla obstáculo ao campo de surfactante.
3. **Perda de dedo** tem mecanismo específico (lição #104): líder que **para** deixa de ser líder,
   e cada dedo é o rastro de um líder (#98). Um pilar que intercepta um líder elimina um dedo
   inteiro — não é afinamento, é extinção.

## 4. Métricas — decididas ANTES, com população e estatística explícitas (§2.5, lição #46)

**Todas em média ± desvio na janela t ∈ [35, 50]**, contra o ruído entre as duas realizações do
P2 como piso. Nenhuma contagem ou massa **absoluta**: a rodada rugosa nasce com −11.6% a −35.4%
de fluido, então só grandezas **por área** ou mascaradas.

- **Primário (responde à pergunta):** `R99`(t=50) e `dR/dt` em t ∈ [25,50].
- **Secundário (como atrapalha):** área da colônia / disco, largura EDT e ciclo em 0.6R e 0.8R,
  contagem de braços, **tortuosidade**.
- **Mecanístico (por que):** fração de líderes que passa a menos de uma largura de braço (5.4 dx)
  de um pilar; quantos deixam de ser líder após o contato (rastreio de identidade da lição #104);
  desfecho classificado em **desvio / divisão / parada**.

## 5. Vigilância do confundidor químico — o risco principal

O `c_s` máximo no ágar vizinho ao corpo do P2R23 é **exatamente 0.300**, encostado no gate
`COL_CS_MIN = 0.3`, com fração acima do gate **0.0000**. Não há folga nenhuma: **qualquer
acúmulo local de `c_s` abre o gate da colonização**, a base é recrutada, o núcleo cresce e a
colônia vira disco — pelas lições #90, #94-H, #97 e #101, por razão **química**, não por bloqueio.

Em R3, `φ = 35%` derruba `D_eff` a ~55% por tortuosidade de meio efetivo. Como o `c_s` do mid-arm
já está a 0.98 do teto (#73), a subida se dá no ágar das baias — exatamente onde o gate está no
limite.

**Reportar obrigatoriamente em toda dose, junto do resultado primário:**

| sentinela | referência P2R23 | leitura |
|---|---:|---|
| `frac(c_s > 0.45)` nas vivas | 0.007 | subiu? o campo está saturando |
| frac. do ágar a < 2h do corpo acima de `COL_CS_MIN` | 0.0000 | **> 0.02 = gate aberto** |
| `R_min` / núcleo sólido | 0.39 R | cresceu? é disco, não bloqueio |

**Regra de interpretação, pré-registrada:** se `dR/dt` cair **e** o gate abrir **e** o núcleo
crescer, a conclusão **não pode** ser "a rugosidade atrapalha o swarm" — é "a rugosidade abriu o
gate da colonização". Nesse caso a série tem de ser repetida com `COL_CS_MIN` reescalado, e a
regra da lição #78-b se aplica: o invariante ao reescalar um limiar absoluto é a **fração da
população que ele seleciona**, nunca a razão ao pico.

## 6. Critérios de ABORTO (matar a rodada, não esperar t=50)

1. `max_cs` > 0.6 na iteração 200 — instabilidade difusiva explícita (lição #96). `D` não muda
   nesta série, então não deve ocorrer; custa zero verificar.
2. `n_pen` > 0 sustentado por mais de um quadro — a força de contato não segura, e a geometria
   perdeu sentido.
3. `dt` médio < 0.5× o do liso no mesmo `t` — colapso de `dt`.
4. `max_v` > 1.0 ou partícula ejetada do domínio.
5. Massa com runaway (razão 2ª/1ª metade > 1.3, guardrail do §2.5).

> **`n_pen` é estrito — vigiar `n_pen_sup` junto (medido em 2026-09-18, ANTES de R1).** A
> definição pré-registrada de `n_pen` é `d(centro) < R_p − 0.5 dx`, e ela lê **0 em todos os
> seis quadros do V1** — a rodada que correu **sem** a guarda de deposição e cujo vazamento
> motivou a guarda. As 6 partículas que o plano §5 registrou como "dentro dos pilares" (5 filler)
> estão a `d/R_p` entre **0.971 e 0.998**: nos vãos do anel de superfície, não no corpo sólido.
> A distância mínima fluido–partícula-de-pilar no V1 é **0.799 dx**, acima do `INSERT_PROX` de
> 0.7 dx.
>
> A definição pré-registrada **fica como está** — mexer nela depois de ver o dado seria escolher
> a régua pelo resultado. Foi acrescentada a coluna `n_pen_sup` (`d(centro) < R_p`, = 6 no V1)
> como **indicador precoce**: é ela que sobe primeiro se material começar a atravessar. Leitura:
> `n_pen` é o aborto; `n_pen_sup` subindo entre quadros é o aviso de que o aborto vem.

## 7. Como o resultado será lido

- **Atrapalha:** `dR/dt` cai monotonicamente com `φ`, τ sobe, **e** as sentinelas químicas da
  seção 5 ficam na referência. A magnitude do bloqueio geométrico é **R3 − R4**, não R3 − liso.
- **Não atrapalha:** `dR/dt` dentro do ruído em **todas** as doses, **com** a fração de
  interceptação de líderes demonstrando que houve contato. Sem essa fração medida, um negativo
  não é publicável — seria a lição #102 de novo.
- **Inconclusivo:** resposta não-monotônica, ou sentinela química fora da referência.

## 8. Pré-passo — EXECUTADO em 2026-09-18

**(A) Empírico, V1 (Λ = 28 dx).** Distância dos líderes (`~filler & rho_b>=0.1 & r >= 0.7·R99`,
a mesma seleção do `_alarga_rastro`) à *superfície* do pilar:

| t | R99 | líderes | d p50 (dx) | d mín (dx) | **< 1 braço** |
|---:|---:|---:|---:|---:|---:|
| 8.2 | 0.958 | 63 | 9.9 | 4.2 | 3% |
| 12.4 | 1.176 | 33 | 7.3 | 1.3 | **42%** |
| 16.6 | 1.403 | 24 | 4.5 | 1.0 | **62%** |
| 20.0 | 1.682 | 23 | 5.9 | 1.3 | **48%** |

Com a colônia mal tendo passado da zona de exclusão (1.2), metade dos líderes já está em alcance
de contato.

**(B) Geométrico, Monte Carlo sobre a rede real** (4000 raios radiais de r=1.2 a 4.067; braço de
5.4 dx toca pilar quando o centro fica a < R_p + w/2). A rede é *ordenada*, então a estimativa de
Poisson por livre caminho médio não se aplica — daí o Monte Carlo:

| Λ | φ | pilares | fresta/braço | ≥ 1 contato | **contatos/braço** |
|---:|---:|---:|---:|---:|---:|
| 28 dx | 11.6% | 92 | 3.3 | **100%** | **1.50** |
| 20 dx | 22.7% | 180 | 1.9 | **100%** | **2.38** |
| 16 dx | 35.4% | 256 | 1.1 | **100%** | **4.02** |

**Conclusão:** o risco de falso negativo por transparência está **afastado por medição** — nenhuma
das doses é transparente, e a escada de dose vale 2.7× em número de contatos. R1 é um ponto
informativo, não uma âncora vazia. O zero é o run liso.

Script: `scratchpad/interceptacao.py` (análise de geometria sobre run existente; não toca o solver).

## 9. Estado da infraestrutura (2026-09-18)

Bloqueantes fechados e verificados: máscara de pilar em `void_fraction`, em `fig_tese.campo` (e
por herança em `rank_runs`, `lit_rank`, `diag_disco`, `cmp_colonia`) e em `ciclo_largura`; guardas
geométricas em `_insere_vacuo`, no `livre()` do wake e no `ParticleShift`. Com `use_pilares=False`
o código é **bit-a-bit idêntico** ao anterior (`runs/swarm/_verif_mask_novo` × `_verif_mask_ctrl`,
10 threads, iterações 0/200/400/600).

**Instrumentação — FEITA em 2026-09-18:** colunas `n_pen`, `n_pen_sup`, `n_contato`, `a_rep_max`,
`a_rep_med` no FIM do `LOG_HEADER` (47 → 52 colunas), mais os arrays `ax_rep`/`ay_rep` acumulados por `ForcaContornoPilar` e **descontados** do resíduo `a_press_max`/`a_press_med` (plano §4: senão a coluna de pressão viraria a coluna de contato). `n_pen` também é impresso a cada quadro no console, com alerta, para o critério de aborto 2 ser verificável em tempo de execução e não só no CSV.
