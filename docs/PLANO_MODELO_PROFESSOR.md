# Plano — adaptações do modelo do professor (Cesgrav 2026) ao nosso modelo

> Escrito em 2026-10-28, branch `feat/rugosidade`. Fonte: *"Smoothed Particle Hydrodynamics models of
> bacterial colony growth in random porous media: percolation geometry, chiral active walkers, and a
> Fresnel–Brewster lattice rule"* (Cesgrav, LCSI/UPE, 13/10/2026). A análise do artigo está na nota do
> Obsidian *"Modelo do professor (1 partícula = 1 célula) × nosso modelo"*; aqui fica só o que entra no
> código, em que ordem e com que critério.
>
> Lições citadas como **#N** (CLAUDE.md §9) e seções de `docs/CRITERIOS_RUGOSIDADE.md` como **§N**.

## 0. Princípio

O modelo dele e o nosso respondem a perguntas em escalas diferentes (microcolônia × placa). **Não
importamos a arquitetura** (uma partícula por célula, EOS só repulsiva, ramificação por falta de
nutriente). Importamos **quatro peças** que são válidas na nossa escala e uma ideia de colaboração:

| # | peça do artigo | onde entra no nosso modelo | mexe no liso? |
|---|---|---|---|
| A | modelo booleano de discos sobrepostos (§4) | geometria do padrão aleatório do Pass L | não |
| B | distância com sinal s(x) até o sólido (§2.7, §4.3) | ferramenta de medida e base de C/D | não |
| C | mobilidade de Kozeny–Carman k(x) perto da parede (§4.3) | controle R4 (atrito do substrato) | não (k = 1 sem pilar) |
| D | parede por mola linear + rejeição (§2.7) | alternativa ao anel de Monaghan & Kajtar | não |
| E | mobilidade dependente de nutriente A(c) (Eq. 8) | ~~substituto do pin químico~~ descartada (§3) | — |
| F | relação de dispersão σ(k) = V\|k\| − Γk² (§3.3) | texto da tese (formaliza a #107) | não |
| G | homogeneização: lei efetiva medida no modelo dele | atrito do substrato no `LinearDrag` | não |

As peças A–D e F–G não alteram o modelo liso: o pareamento rugoso × liso com o P2R23 (mesma semente,
mesmo número de threads) continua valendo. A peça E foi descartada sem rodada (§3).

---

## 1. Fase 0 — offline, sem rodada

### F0.1 Campo s(x) analítico (peça B)
- Pilares circulares: `s(x) = min_i (|x − c_i| − R_s)`. Para o booleano (A) é o mesmo mínimo.
- **Definir a superfície R_s.** Hoje há três raios (anel 4.6 dx, densidade até 4.9 dx, fluido a partir
  de 5.4 dx). Medir nos quadros do R1″ a distribuição da distância das partículas de fluido ao centro do
  pilar e fixar R_s onde a densidade cai a metade. Sem isso, k(x) e a mola não têm referência.
- Uso imediato: reescrever as medidas de contato (P(parada|contato), §20.6/§21.6) em função de s, em
  vez de "distância ao centro menos raio".

### F0.2 Geometrias booleanas (peça A)
- Sortear centros por Poisson com densidade λ tal que `φ = exp(−λπR²)` iguale as porosidades da série
  em rede: **R1 0.88, R2 0.77, R3 0.65**. Semente fixa (`PILAR_SEED`).
- **Filtrar gargalos** mais estreitos que o suporte do kernel (2h = 3.6 dx): fundir os dois discos ou
  remover um. Medir quantos são removidos e quanto φ muda; φ final tem de ficar a ±0.01 do alvo.
- **Podar o anel de contato** dentro dos discos vizinhos (partícula de anel com s < 0 em relação a
  outro disco sai). Conferir lisura da parede composta como no R1″ (alvo: desvio ≤ 0.1%).
- Conferir que nenhum disco cai na exclusão do inóculo (r < 1.2).

### F0.3 Escala do Kozeny–Carman (peça C)
- No artigo `d_p` é o diâmetro da célula; aqui a escala equivalente é **dx** (espaçamento de partícula).
  Com `k = min(1, max(0.2, s/dx))` o atrito age só na primeira camada junto à parede.
- Medir no R1″ quantas partículas vivas e de filler ficam a s < 1 dx (por quadro, janela t ∈ [35, 50]):
  é a população que a alavanca toca. Se for < 1% das vivas, prever efeito nulo e reconsiderar a escala
  (alternativa: `s/(2h)`, que alcança 3.6 dx).

### F0.4 Estabilidade do atrito alto (peça C)
- A lição #85 mediu que o `dt` do baseline é limitado pelo arrasto (2/γ_eff ≈ 0.013). Dividir γ por
  k = 0.2 baixaria o `dt` global 5×, porque o mínimo é sobre todas as
  partículas. **Integrar só a parte extra do arrasto de forma implícita**,
  `v ← v / (1 + dt·γ_extra)`, com `γ_extra = γ_eff(1/k − 1)`; onde k = 1 o fator é exatamente 1 e o
  liso fica bit-a-bit. É a forma incondicionalmente estável já usada no M-B.9a, tentativa 2 (#22).
- Conferir antes, num teste de 200 iterações, que o `dt` do PySPH de fato inclui o arrasto (se não
  incluir, o esquema explícito com γ/k explode, como no M-B.9a, tentativa 1).

### F0.5 Distribuição de `c_n` das pinadas (peça E) — dispensada, ver E1 no §3

### F0.6 Texto (peça F)
- Escrever na tese a relação de dispersão e por que os dois termos faltam no nosso modelo
  (desestabilizador preso pelo teto de `c_s`; estabilizador ausente porque a pressão está 219× abaixo
  da Marangoni, #85, #107). Sem código.

---

## 2. Fase 1 — dentro da série de rugosidade (não muda o liso)

Ordem: termina-se a série em rede (R3, Λ = 16) **antes** de qualquer item abaixo, para não misturar
alavancas na mesma comparação.

### R4 — atrito de Kozeny–Carman (peça C)
- **Uma alavanca:** geometria do R1″ + `k(x)` no `LinearDrag`, parte extra implícita (F0.4).
- **Comparação:** R4 × R1″ (efeito do atrito) e R4 × P2R23 (rugosidade + atrito × liso).
- **Fundamentação:** a parede tem condição de não deslizamento (Violeau 2012, cap. sobre contornos;
  Liu §6.5 para o kernel truncado junto ao sólido); k(x) é a forma de lubrificação do artigo. Comparar
  na discussão com o atrito proporcional à área de contato de [T15].
- **Pré-registrar antes de rodar** (formato das §20/§21): P(parada|contato) absoluto; R99 médio em
  t ∈ [35, 48] interpolado; dR/dt; motor `a_mar_bio_p95`; dedos, AR, largura e núcleo pela régua com
  os pilares descontados; tudo contra o piso de ruído do par de realizações do P2.
- **Teste de efeito antes** (#102): 600 iterações com e sem a alavanca; se o estado for bit-a-bit
  igual, a alavanca não tem população e o R4 não roda.

### R5–R7 — série booleana (peça A)
- Mesmas porosidades de R1–R3, geometria de F0.2. Pergunta: **a mesma porosidade em arranjo aleatório
  freia mais que em rede?** (a rede tem corredores retos; o booleano tem gargalos e becos).
- Duas realizações por porosidade se o custo permitir (~18 min cada); o artigo usa uma só, e isso foi
  uma das críticas que levamos a ele.
- Mesmos critérios do R4; acrescentar fração de pilares tocados e P(parada|contato) por tamanho de
  gargalo.

### Opcional — contorno por mola + rejeição (peça D)
- Só se o custo dos pilares virar gargalo (o R2 tem 13 260 partículas de densidade e 20 880 de anel).
- Implementação: força `a_W = k_w (s0 − s)₊ n̂` com `n̂ = ∇s` analítico, e projeção que anula a
  componente da velocidade que entra no sólido quando s < 0 (a mesma ideia da `ProjecaoShiftParede`,
  aplicada à velocidade). **Manter a grade de densidade**: a nossa `BiomassEOS` tem ramo atrativo, e
  sem suporte de kernel completo a parede atrairia o fluido (Liu §6.5; é por isso que o artigo, com EOS
  só repulsiva, dispensa a grade e nós não).
- **Critério de equivalência:** refazer o R1″ com a mola; aceitar se R99, dR/dt e P(parada|contato)
  ficarem dentro do piso de ruído. Só então trocar, e só numa série nova.
- **Ressalva:** a rejeição quebra a conservação de momento; aceitável porque a dinâmica já é dominada
  por arrasto e pin (mesmo argumento usado para a KGC), mas tem de ser dito.

---

## 3. Fase 2 — série nova sobre o modelo liso (depois do Pass L)

### E1 — mobilidade A(c) no lugar do pin químico (peça E) — DESCARTADA (2026-10-28, sem rodada)
- **Não substitui o pin químico, pela mesma razão da molhabilidade.** O pin químico segura 55–64
  vivas, todas no núcleo (mediana r/R99 = 0.07); a função dele é ancorar o núcleo. Para A(c) fazer
  isso, teria de dar mobilidade ≈ ε no núcleo e ≈ 1 no resto.
- **No regime rico em nutriente isso é impossível.** Com `k_src` = 0.3 o `c_n` nunca cai abaixo de
  0.375, e o do núcleo fica em 0.41 (#82). Uma Michaelis–Menten não separa 0.41 de 1.0:

  | `K_c` | A(0.41) núcleo | A(1.0) fora | razão núcleo/fora |
  |---:|---:|---:|---:|
  | 0.1 | 0.81 | 0.91 | 0.89 |
  | 1 | 0.33 | 0.53 | 0.62 |
  | 3 | 0.16 | 0.29 | 0.57 |
  | 7.4 | 0.10 | 0.16 | 0.61 |
  | 20 | 0.07 | 0.10 | 0.73 |

  O contraste máximo é ~1.75×, contra 0 / 1 do pin. Para o núcleo chegar a 0.1 (`K_c` = 7.4) a colônia
  inteira fica com mobilidade ≤ 0.16: é um arrasto ~6× maior em tudo (a Saída B da #86, que perdeu o
  relevo), não uma âncora. Com `K_c` pequeno a mobilidade fica ~0.8–0.9 em todo lugar e o núcleo solta,
  que é o P19 (#83).
- **A forma do artigo funciona lá porque o regime é outro:** limitado por nutriente, com `c` indo a
  zero no interior. Isso é biofilme ([T2]); no nosso swarm, rico em nutriente, não acontece.

### E2 — kernel Wendland C2
- O artigo usa Wendland; nós já testamos (P2W, #89): resolve a cascata de pareamento e leva o run a
  t = 100 sem mudar a morfologia. Como a janela numérica do P2R23 é t ≲ 60 (#88), adotar o Wendland
  ao abrir uma série nova, só se a usuária quiser estender a janela.

---

## 4. Colaboração — homogeneização (peça G)

- O modelo dele, na escala do poro, mede uma lei efetiva: mobilidade da colônia em função da
  porosidade φ e do tamanho de poro (análoga ao D_eff(φ) da Tabela 2 dele).
- Essa lei entra no nosso `LinearDrag` como atrito do substrato γ_sub(φ), sem pilares explícitos. É
  um **R4 alternativo**, de base física independente, e põe o trabalho dele na tese sem trocar a
  escala dela.
- Pedido a fazer a ele: rodar 3–4 porosidades (0.88, 0.77, 0.65 e uma perto de φ_c ≈ 0.324) e medir
  a velocidade da frente normalizada pela da colônia livre, **no mesmo t** (a Tabela 1 compara N em
  tempos diferentes).

---

## 5. O que não entra

| do artigo | por que não |
|---|---|
| uma partícula = uma célula | ~10⁹ partículas na escala da placa |
| nascimento de partículas com massa fixa | é a nossa mitose (#81, #85, #86): povoa, não empurra a ponta; no máximo, como fonte da fase passiva (série O de `PLANO_FASE_PASSIVA.md`) |
| EOS só repulsiva | quebra a coesão que segura o pescoço (#15, #23) |
| mobilidade A(c) no lugar do pin | `c_n` ≥ 0.375 no regime rico em nutriente: contraste núcleo/fora ≤ 1.75× (§3) |
| ramificação por falta de nutriente | regime de biofilme, não de swarming de PA14 ([T2]) |
| nutriente em grade euleriana | o fluxo zero pelo Laplaciano de Brookshaw já foi medido no R1 (sombra 1.66 contra 1.69 previsto) |
| quiralidade, Fresnel–Brewster, n fendas | sem campo de orientação; dendrites de PA14 radiais; o próprio artigo mostra o sinal quiral apagado em φ ≈ 0.7 |
| expoentes de percolação | exigem muitas realizações perto do limiar |

---

## 6. Ordem e custo

| etapa | rodadas | custo aproximado | depende de |
|---|---:|---|---|
| F0.1–F0.6 | 0 | horas de análise | — |
| R3 (já planejado) | 1 | ~18 min | — |
| teste de efeito do R4 + R4 | 1 (+2 min) | ~20 min | F0.1, F0.3, F0.4, R3 |
| R5–R7 booleanos | 3 a 6 | 1–2 h | F0.2, R3 |
| mola + rejeição (opcional) | 1 | ~15 min | F0.1 |
| E2 (opcional) | 1 | ~30 min | fim do Pass L |

Cada rodada: autorização antes, trava `.sph_rodando`, 10 threads (ou o mesmo número do par), janela
t = 50, pré-registro escrito antes em `docs/CRITERIOS_RUGOSIDADE.md`, e resultado registrado depois,
inclusive as predições que errarem.
