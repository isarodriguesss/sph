# Validação do modelo — testes, referências e o que cada um significa

Rodada de referência: **`runs/C4_t100`** (baseline C4, domínio `[-7,7]` com 261², t=100 s).
Reproduzir com:

```bash
python tools/validate_model.py runs/C4_t100
```

Cada teste abaixo cita a referência da literatura, o valor esperado e o medido.
Onde o modelo diverge, o texto diz por quanto e em que regime — não há aprovação cega.

---

## 1. Lei de expansão — o teste mais forte

### O que é

Como o raio da colônia cresce no tempo: `R(t) ~ t^α`. O expoente α é uma
**assinatura de mecanismo**, não um parâmetro ajustável — ele decorre de qual
processo limita a expansão.

### Referências

| fonte | regime | α previsto |
|---|---|---:|
| **[T2]** Srinivasan, Kaplan & Mahadevan 2019, *eLife* 8, e42697 | swarming *nutrient-rich*: capilaridade dominada, steady-state com velocidade constante `V = C₁·g₀·H·Ca^(1/3)` | **1.0** |
| **[T3]** Giverso, Verani & Ciarletta 2016, *Biomech Model Mechanobiol* 15, 643 | dedos difusão-limitados | **0.45** |

### Medido

| janela | α | R² | n |
|---|---:|---:|---:|
| estabelecido (t = 33–55 s) | **1.03** | 0.997 | 7 |
| pós-esgotamento (t = 55–99 s) | **0.39** | 0.971 | 14 |

### O que significa

O modelo reproduz **os dois regimes e a transição entre eles**. Enquanto há
nutriente, expande a velocidade constante como Srinivasan prevê para swarming;
quando o nutriente esgota, migra para o expoente difusão-limitado de Giverso.

Isso é uma validação **estrutural**, não um ajuste: nenhum parâmetro foi calibrado
para produzir esses expoentes. Eles emergem do balanço Marangoni + arrasto +
crescimento gateado por `c_n`.

### Ressalva

O `t` da simulação **não está calibrado em segundos biológicos**. A comparação
válida é do expoente, que é adimensional. Comparar velocidades absolutas exigiria
ancorar a escala de tempo numa taxa experimental.

---

## 2. Morfologia dendrítica

### O que é

Número de dedos, razão de aspecto e ancoragem do núcleo — a forma que a
instabilidade de Mullins–Sekerka seleciona.

### Referências

- **`reference.jpg`** (Michiels et al., *P. aeruginosa* PA14 em ágar swarming):
  15–20 dendritos radiais, separados por ágar limpo, `AR ≥ 1:5`.
- **[T1]** Trinschek, John & Thiele 2018, *Soft Matter* 14, 4464, painel (b)
  *Fingering*: 7–9 dedos finos partindo de núcleo compacto.

### Medido

| | valor | alvo |
|---|---:|---|
| dendritos (t = 55 s) | **14–18** | 15–20 (PA14) |
| AR = (R − r_núcleo)/largura | **5.7** | ≥ 5 (§2.2) |
| núcleo pinado | r = 0.27 | ancorado |

### O que significa

A contagem **cai de ~45 para ~18** ao longo do run e entra na faixa do PA14. Os 45
iniciais são perturbação de ruído; a competição entre dedos elimina os perdedores
até o número observado experimentalmente. **Isso é a seleção de Mullins–Sekerka
acontecendo** — o modelo não foi semeado com 18 dedos, ele converge para 18.

O núcleo permanece em r = 0.27 (o raio da gaussiana inicial) — confirma que o hard
pin K.17 ancora a matriz madura, como na biologia do PA14.

---

## 3. Halo de surfactante

### O que é

O campo de surfactante deve **extrapolar** a biomassa — o ramnolipídeo difunde no
ágar à frente das bactérias.

### Referência

**[T1]** Trinschek 2018, painel (b): o campo Γ estende-se além da biomassa num
halo radial. É um dos critérios de sucesso de §2.2.

### Medido

| | valor |
|---|---:|
| raio da biomassa (ρ_b > 0.05, p99) | 5.28 |
| raio do surfactante (c_s > 0.02, p99) | **5.83** |
| halo | **+0.55 (1.10×)** |

### O que significa

O halo existe. Ele apareceu quando `lambda_eff` foi corrigido para decaimento
**lento no ágar** (0.5λ) e rápido no biofilme (2λ) — a passagem B3. Antes disso o
código tinha o oposto e o `c_s` morria na borda.

Vale registrar o mecanismo, porque é contraintuitivo: sob a saturação de T1
(`produção ∝ 1 − c_s/c_s_max`), o `c_s` **interior** é quase insensível a λ, mas o
do ágar — onde a produção é ~zero — é linear em 1/λ. Por isso quadruplicar o tempo
de vida no ágar cria o halo sem esvaziar o interior.

---

## 4. Consistência SPH — o operador é confiável?

### O que é

Toda força do modelo é uma soma SPH sobre vizinhos. Se a vizinhança for pobre, o
gradiente calculado é espúrio e a comparação com Trinschek perde sentido.

### Referências

- **[T6]** Liu & Liu 2003, §3.3: ~**20** vizinhos para gradiente confiável,
  ~**35** para Laplaciano.
- **[T7]** Violeau 2012, §3.4/3.6: partição da unidade
  `σ_a = Σ_j V_j W_aj → 1`; **`σ_a < 0.85` ≈ 15% de erro** nos operadores.

### Medido

Esperado numa rede uniforme: `π(2h)²/dx² = 41` vizinhos.

| região | vizinhos | ≥20 | ≥35 |
|---|---:|---:|---:|
| núcleo (ρ_b > 0.8) | 37.0 | 100% | 91% |
| braços (ρ_b 0.1–0.5) | 46.4 | 100% | 95% |
| ágar interior | 40.2 | 100% | 100% |

`σ_a` na colônia: média **0.916**, **6.3%** abaixo de 0.85, p05 = 0.842.

### O que significa

O critério de gradiente de Liu é satisfeito em **100%** da colônia, e o de
Laplaciano em 91–100%. O operador de Marangoni é confiável onde ele atua.

Os 6.3% abaixo do limiar de Violeau concentram-se na junção núcleo–braço, atrás da
frente ativa. A **KGC** (Bonet & Lok 1999) corrige o gradiente mesmo sob suporte
incompleto — validada independentemente por [`test_kgc_consistency.py`](../test_kgc_consistency.py),
que reproduz o gradiente de um campo linear a ~1e-16 sob suporte truncado.

---

## 5. Cobertura espacial (Critério §2.5 C1)

### O que é

Fração da **área** da colônia sem nenhuma partícula próxima. Diferente de `σ_a`,
que só existe onde há partícula — um buraco geométrico é invisível a qualquer média
feita "onde há partícula".

### Referência

**[T7]** Violeau §3.4 — a representação SPH de um fluido exige cobertura espacial
contínua. Vácuo físico na zona de expansão é erro de consistência fundamental.

### Medido

| limiar | área sem vizinho |
|---|---:|
| > 0.7 dx | 5.67% |
| > 1.0 dx | 0.82% |
| **> 1.5 dx** | **0.12%** (33 dx²) |

### O que significa

0.12% da área, ou ~33 partículas de área numa colônia de ~23 000 dx². Contra
**1.28%** (301 dx²) do baseline anterior S4 — uma redução de 9×. O vácuo deixou de
ser bloqueio de roadmap.

---

## 6. Conservação e estabilidade numérica

### O que é

A massa cresce porque inserção adiciona células e `BiomassGrowth` cresce biomassa —
isso é físico. O que **não** pode acontecer é crescer sem limite (runaway) ou a
pressão explodir (over-packing).

### Referências

- **§8** do CLAUDE.md — orçamento de acelerações: `a_pressure` < 3–4.
- **Lição #34** — runaway de massa por realimentação da inserção.
- **[T6]** Liu §6.4 — instabilidade de tração sob pressão negativa.

### Medido

| | valor |
|---|---:|
| massa | 197.5 → 212.9 (**+7.8%**) |
| aceleração de massa (2ª/1ª metade) | **0.77** (desacelerando) |
| `a_pressure` mediana | **3.00** |
| picos > 4 | 10% das amostras |

### O que significa

A massa **convergiu** — a razão 0.77 indica que a taxa está caindo, o oposto de
runaway. A pressão fica no orçamento de §8, com picos ocasionais e isolados.

---

## Resumo

| # | teste | resultado |
|---|---|---|
| 1 | lei de expansão | α = 1.03 (nutrient-rich) e 0.39 (difusão-limitado) — **dois regimes** |
| 2 | morfologia | 14–18 dendritos, AR 5.7, núcleo ancorado |
| 3 | halo | presente, 1.10× |
| 4 | consistência SPH | Liu 100% (gradiente); Violeau 6.3% em déficit |
| 5 | cobertura | vácuo 0.12% (9× melhor que o baseline anterior) |
| 6 | conservação | massa convergida, pressão no orçamento |

### Limites conhecidos

1. **Janela útil t ≲ 55 s.** `c_n` inicia em 1.0 sem fonte; a colônia consome todo o
   reservatório e congela (`mean_v` cai 22×).
2. **Escala de tempo não calibrada** — expoentes são comparáveis, velocidades não.
3. **Réplica única.** Com `SEED = None`, 3–5 réplicas dariam a dispersão de α e da
   contagem de dendritos.
