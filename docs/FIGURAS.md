# Guia das figuras — o que cada painel mostra e como ler

Três geradores, cada um respondendo a uma pergunta diferente.

| script | pergunta | saída |
|---|---|---|
| [`tools/plot_expansion.py`](../tools/plot_expansion.py) | **como** a colônia cresce? | `runs/expansao.png` |
| [`tools/plot_fields.py`](../tools/plot_fields.py) | **por que** ela cresce assim? | `runs/campos.png` |
| [`tools/compare_frames.py`](../tools/compare_frames.py) | qual rodada é melhor? | `runs/comparativo.png` |

---

## 1. `expansao.png` — dinâmica da expansão

```bash
python tools/plot_expansion.py runs/C4_t100 --window 5 55
```

A `--window` marca a faixa em que a rodada é confiável (antes do nutriente esgotar).
Ela aparece sombreada em (b) e como linha tracejada vermelha nos demais.

### (a) Lei de expansão vs literatura — *o painel principal*

`R(t)` em log-log. Numa escala log-log, uma lei de potência `R ~ t^α` vira **reta
de inclinação α** — por isso o eixo é log dos dois lados.

- **pontos pretos** — o modelo
- **linhas cheias** — ajustes em duas janelas, com α e R² na legenda
- **linhas tracejadas** — as inclinações de referência da literatura, ancoradas no
  primeiro ponto para servirem de régua visual

**Como ler:** se os pontos correm paralelos à tracejada verde, o modelo está no
regime de Srinivasan (velocidade constante); se correm paralelos à vermelha, está
no de Giverso (difusão-limitado). Aqui eles fazem as duas coisas, em sequência.

### (b) A transição de regime é o nutriente

`R(t)` em escala linear (preto, eixo esquerdo) com `c_n` nos braços sobreposto
(azul, eixo direito).

**Como ler:** a linha pontilhada azul horizontal é `c_n = 0.4`, o piso do gate de
crescimento. O ponto em que a curva azul a cruza é onde o crescimento se desliga —
e é exatamente onde a curva preta muda de inclinação. **Este painel é a explicação
causal do painel (a).**

### (c) Kymograph — cada faixa clara é um dendrito

Eixo x = ângulo (−180° a 180°), eixo y = tempo, cor = raio da frente naquele setor.

**Como ler:** uma coluna clara que persiste de baixo para cima é **um dendrito
sobrevivendo**. Uma que aparece e some é um dedo que perdeu a competição. O
alargamento das faixas com o tempo é o engrossamento dos braços.

É a melhor vista para ver **competição entre dedos**, que é o mecanismo de
Mullins–Sekerka em ação.

### (d) Contorno da frente — cor = tempo

O perfil `R(θ)` desenhado em coordenadas cartesianas, em 12 instantes, colorido
pelo tempo (colorbar à direita).

**Como ler:** contornos concêntricos e regularmente espaçados = expansão uniforme.
Contornos que se acumulam em algumas direções e não em outras = dedos avançando
enquanto baías estagnam. As "pétalas" são os dendritos.

> Versões anteriores usavam preenchimento com legenda de instantes; ficava ilegível
> por sobreposição. Só contorno com colorbar é mais limpo.

### (e) Seleção competitiva de dendritos

Contagem de picos no perfil angular, ao longo do tempo, com as faixas de referência
sombreadas: verde = 15–20 (PA14, `reference.jpg`), roxo = 7–9 ([T1] Trinschek).

**Como ler:** a curva descendo mostra dedos sendo eliminados. Entrar e permanecer
na faixa verde é o resultado desejado. **Este painel mostra que 18 não foi
imposto — foi selecionado.**

### (f) Velocidade da frente

`dR/dt` no tempo, com a mediana da janela útil como linha tracejada.

**Como ler:** um patamar horizontal é o steady-state de [T2] — velocidade constante.
A queda após a linha vermelha é o congelamento por falta de nutriente.

---

## 2. `campos.png` — os três escalares, como campo contínuo

```bash
python tools/plot_fields.py runs/C4_t100 --times 15 40 55 99 [--n 420]
```

Linhas = instantes, colunas = os três campos. Escalas de cor **fixas** entre
instantes, para que a comparação temporal seja honesta.

### Por que campo contínuo e não pontos

Os campos são reconstruídos por **interpolação SPH Shepard-normalizada**
(Price 2007, `splash`, PASA 24:159):

```
A(x) = Σ_j V_j W(|x−x_j|, h) A_j  /  Σ_j V_j W(|x−x_j|, h)
```

A normalização pelo denominador é o que torna o campo contínuo mesmo onde a
amostragem é irregular — é a mesma razão pela qual `σ_a` mede consistência.

Isso importa conceitualmente: **renderizar por pontos mostra as partículas;
renderizar o campo mostra o fluido**, que é o objeto físico. As partículas são o
método numérico, não o fenômeno. Uma versão anterior desta figura usava scatter e
ficava abstrata — a colônia parecia um aglomerado de pontos em vez de um tecido.

> Isto é pós-processamento puro. Nunca adicionar partículas ao solver para
> "preencher" visualmente — ver lições #34 e #36.

### Coluna 1 — biomassa `ρ_b` (verde)

A **estrutura** da colônia.

- contorno externo em `ρ_b = 0.1` — a borda da colônia
- contorno interno em `ρ_b = 0.8` — o **limiar do hard pin** (K.17): dentro dele o
  núcleo maduro é imóvel

**Como ler:** o verde escuro central é o núcleo pinado; os braços em verde médio são
swarmers móveis. Comparando as linhas, o núcleo **não cresce** — ele é a gaussiana
inicial, e toda a expansão vem do transporte dos braços.

### Coluna 2 — surfactante `c_s` (amarelo→vermelho)

O **motor**. A força é `F = −β·∇c_s`, então o que move a colônia é o **gradiente**,
não o valor — um campo alto e uniforme não empurra nada.

O contorno **verde** sobreposto é a borda da biomassa (`ρ_b = 0.1`), desenhado
justamente para que se possa ver o surfactante **além** dela.

**Como ler:** procure a região laranja que ultrapassa o contorno verde — é o
**halo** de [T1] Trinschek, painel (b). Os pontos vermelho-escuros nas pontas são
picos locais de produção, onde o gradiente é mais forte.

### Coluna 3 — nutriente `c_n` (azul)

O **limitador**. Reservatório finito: inicia em 1.0 em todo o domínio e não tem
fonte.

O contorno **vermelho** marca `c_n = 0.4`, o piso do gate de crescimento — dentro
dele o crescimento está desligado.

**Como ler:** azul escuro = ágar virgem, branco = esgotado. Acompanhe o contorno
vermelho entre as linhas: ele nasce como um círculo pequeno no núcleo e vai
engolindo a colônia inteira. Em t=99 s ele contém tudo — e é por isso que a
colônia parou.

### A leitura conjunta

> `ρ_b` estrutura a colônia · `c_s` a move (via `∇c_s`) · `c_n` a limita.

Em t=99 s a coluna 2 ainda mostra surfactante e a coluna 1 ainda mostra braços,
mas a coluna 3 está toda dentro do contorno vermelho. O motor existe, a estrutura
existe — falta combustível.

## 3. `comparativo.png` — qual rodada é melhor

```bash
python tools/compare_frames.py runs/C1 runs/C2 runs/C4
```

Linhas = rodadas, colunas = 4 diagnósticos. Primeira linha traz as duas referências
obrigatórias de §2.2 (`reference.jpg` e `reference_result.png`) como âncora visual.
Escalas de cor **fixas entre rodadas** — sem isso a comparação é inválida.

| coluna | mostra | usar para |
|---|---|---|
| swarm (`rho`) | empacotamento SPH | morfologia: continuidade dos braços, baías limpas |
| biomassa vs inseridas | verde = real, vermelho = inserida | quanto da colônia é biomassa vs suporte |
| `σ_a` | consistência do kernel, `<0.85` circulado | onde o operador é confiável |
| `c_s` | surfactante (só partículas químicas) | motor e halo |

**Como ler:** a coluna `σ_a` é a mais decisiva — a **fração** circulada de azul
(no título) é o critério C2 de §2.5. Julgar pela média em vez da fração já levou a
ranking errado (lição #46).

---

## Ordem sugerida de leitura

1. **`campos.png`** — entender o mecanismo (o que cada campo faz)
2. **`expansao.png` (b) e (a)** — ver a consequência dinâmica e a comparação com a literatura
3. **`expansao.png` (c) e (e)** — ver a seleção competitiva
4. **`comparativo.png`** — decidir entre configurações
