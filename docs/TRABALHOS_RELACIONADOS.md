# Trabalhos relacionados — pares metodológicos e pares de ESCALA

> Levantado em 2026-09-18. Pergunta: quais trabalhos publicados se parecem com o nosso, e quais
> operam na **mesma escala** (nossa partícula SPH ≈ 0.4-0.6 mm de aresta, colônia de 3-5 cm).
> Complementa os [T1]-[T14] do §3.0 do CLAUDE.md; a numeração [T15]+ é proposta, ainda não
> incorporada lá.

## 1. A nossa escala, em número

Ancorando `dx` pela largura de braço de PA14 ([T11] Deng et al.: 2-5 mm ≈ 5.6-8.7 `dx`, lição #100):

| | valor |
|---|---|
| `dx` | **0.4-0.6 mm** |
| área por partícula (`dx²`) | 0.16-0.36 mm² |
| `R99` em t=50 | 83.6 `dx` ≈ **33-50 mm** (uma placa de 9 cm inteira) |
| domínio `[-7,7]` | ~100-155 mm de lado |

**Quantas bactérias cabem numa partícula nossa?** Tomando a densidade que Ben-Jacob usa para uma
placa ao fim do crescimento (10⁹-10¹⁰ células numa placa de 9 cm ⇒ 1.6×10⁵-1.6×10⁶ células/mm²):

> **~3×10⁴ a 6×10⁵ bactérias por partícula SPH.**

Isso é a mesma ordem de grandeza do *communicating walker* de Ben-Jacob (**10⁴-10⁵ células por
walker**) e **10²-10³ vezes mais grosseiro** que o agente de Du et al. (≈10² bactérias). É o número
que localiza o nosso trabalho na literatura, e vale declará-lo na tese: **a partícula SPH não é uma
célula, é um retalho de colônia** — a mesma decisão de granularidade que a escola de Ben-Jacob tomou
em 1994, só que com hidrodinâmica de verdade em vez de passeio aleatório.

## 2. A matriz que organiza tudo

|  | **escala de CÉLULA (µm)** | **escala de COLÔNIA (mm-cm)** |
|---|---|---|
| **método de PARTÍCULAS** | NUFEB/IbM, DEM, DPD, SPH de adesão (2025), SPH pore-scale (Tartakovsky), SPH de biofilme (Soleimani) | **communicating walkers** (Ben-Jacob 1994) · agentes acoplados a thin-film (Du et al. 2011) · **← nosso trabalho** |
| **método CONTÍNUO** | — | Trinschek [T1], Srinivasan [T2], Giverso [T3], modelo generalizado 2025, Kawasaki-Matsushita, campo de fase |

**A célula que o nosso trabalho ocupa — partículas + escala de colônia + hidrodinâmica SPH com
Marangoni — não tem nenhum artigo publicado que eu tenha encontrado.** Os dois vizinhos mais próximos
(Ben-Jacob e Du et al.) usam agentes estocásticos ou elásticos, não SPH. É a lacuna que a tese preenche,
e é declarável como tal.

## 3. Pares de ESCALA — os que você pediu

### [T15] Ben-Jacob, Shochet, Tenenbaum, Cohen, Czirók & Vicsek (1994), *Nature* 368:46-49
*Generic modelling of cooperative growth patterns in bacterial colonies.* **O precedente conceitual
mais direto.** Cada "walker" representa **10⁴-10⁵ bactérias** (porque uma placa tem 10⁹-10¹⁰ células e
não se simula célula a célula), com 10⁴-10⁶ walkers por experimento numérico. Walkers consomem
nutriente, reproduzem, caminham com viés quimiotático e são confinados por um envelope. Reproduz
padrões ramificados em função de nutriente e dureza do ágar. **É a justificativa histórica para
partícula ≠ célula.** Extensão com lubrificação: Kozlovsky, Cohen, Golding & Ben-Jacob (1999),
*Phys. Rev. E* 59:7025 — *Lubricating bacteria model for branching growth of bacterial colonies*,
que é o antecessor direto da ideia de fluido secretado dirigindo a expansão.

### [T16] Du, Xu, Shrout & Alber (2011), *Math. Models Methods Appl. Sci.* 21(S1):939-954
*Multiscale modeling of Pseudomonas aeruginosa swarming*, doi:10.1142/S0218202511005428.
**O par mais próximo em FENÔMENO: mesma espécie, mesmo motor, mesma morfologia.** Acopla três
submodelos — (i) equação de filme fino por lubrificação com tensão superficial modulada por
ramnolipídeo e **tensões de Marangoni**; (ii) convecção-difusão-reação para nutriente, sinal de quorum
sensing e ramnolipídeo; (iii) modelo discreto estocástico off-lattice de bactérias como corpos
elásticos de três nós (0.3-0.8 × 1.0-1.2 µm), com flagelo, crescimento e divisão. ~10⁵ agentes, cada um
≈10² bactérias; malha 200×200 num domínio 2π×2π; camada líquida inicial de 100 µm. **Reproduz tendrils
na borda do enxame** e mostra que o mutante *rhlAB* (sem ramnolipídeo) não espalha.

**Ressalva de escala, e ela é favorável a nós:** 10⁵ agentes × 10² bactérias = **10⁷ células**, ou seja
~1% de uma placa. Du et al. resolvem **a borda do enxame em resolução quase-celular**; nós resolvemos
**a colônia inteira em resolução de retalho**. Não competem — são níveis de descrição diferentes do
mesmo sistema, e isso é exatamente como posicionar os dois na revisão.

### [T17] Modelo generalizado de swarming em substrato poroso (PMID 39655366, publicado jan/2025)
*A generalized model for predicting different morphologies of bacterial swarming on a porous solid
surface.* Modelo bifásico (fase celular + fase aquosa) em aproximação de filme fino, com produção de
surfactante, arrasto entre fases, **influxo osmótico**, tensões de Marangoni, pressão de disjunção
para molhabilidade e filme precursor para a singularidade da linha de contato. **Prediz sete
morfologias: arrested, circular, modulated, branching, droplet, fingering e dendrite** — é o sucessor
direto do Trinschek [T1], com o influxo osmótico de Srinivasan [T2] embutido. Mostra que maior
molhabilidade acelera a expansão e maior tensão superficial redistribui biomassa radialmente.
**Autores e periódico não confirmados** — não consegui a página do registro; só o PMID e o resumo.
É o alvo de leitura prioritário: cobre exatamente o eixo `cs_max`/wettability que o §12 do CLAUDE.md
identificou como o bloqueio aberto do projeto.

### [T18] Kawasaki, Mochizuki, Matsushita, Umeda & Shigesada (1997), *J. Theor. Biol.*
*Modeling spatio-temporal patterns generated by Bacillus subtilis.* Reação-difusão contínua na escala
de colônia (mm-cm), com difusividade bacteriana dependente de nutriente e de dureza do ágar. Reproduz
o diagrama morfológico de cinco classes de Matsushita (DLA-like, Eden-like, anéis concêntricos, disco,
dense-branching). **É a referência canônica de "morfologia de colônia em placa por modelo contínuo"**
e o par histórico do nosso diagrama de fases rugoso × liso.

### [T19] Porter, Trenado-Yuste, Martinez-Calvo, Su, Wingreen, Datta & Huang (2025), *Nat. Rev. Phys.* 7:535
*On the growth and form of bacterial colonies.* Revisão recente de modelos biofísicos de agregados
bacterianos. **Uso recomendado: enquadramento do capítulo de revisão da tese** — é a fonte que
justifica, com autoridade de 2025, tratar a colônia como agregado multicelular denso e não como soma
de células isoladas, que é a premissa da nossa granularidade.

## 4. Pares de MÉTODO — SPH aplicado a bactérias (todos em escala fina)

Servem para fundamentar **a escolha do SPH**, não para comparar morfologia. Todos operam µm a mm.

| trabalho | o que faz | escala |
|---|---|---|
| **[T20]** Soleimani, Wriggers, Rath et al. (2016), *Comput. Mech.* 58:619-633, doi:10.1007/s00466-016-1308-9 | Biofilme 3D em SPH contínuo: crescimento por reação biológica + difusão de nutriente, deformação e erosão da interface por escoamento (FSI). Validação experimental. | biofilme, µm-mm |
| **[T21]** *Simulation of bacterial adhesion on a rough surface based on SPH*, **Comput. Part. Mech.** (2025), doi:10.1007/s40571-025-00924-1 | SPH (PySPH) + Lennard-Jones para adesão de uma célula rígida tipo *S. aureus* sobre topografia medida. | célula única, µm |
| **[T22]** Tartakovsky et al. (2009), *J. Porous Media* 12(5) | SPH lagrangiano multi-escala: escoamento poro-escala, transporte reativo e crescimento de biomassa por Monod duplo; adesão/destacamento por forças de par. | poro, µm |
| **[T23]** *A conservative SPH method for surfactant dynamics* | Formulação conservativa de transporte de surfactante em SPH (bulk + interface, adsorção/dessorção de Langmuir) com tensão superficial acoplada à concentração local e **forças de Marangoni**. | interface, genérico |
| **[T24]** *An alternative SPH formulation to simulate chemotaxis in porous media*, PMC5388734 | Quimiotaxia bacteriana em SPH. | poro, µm |

**[T20] é relevante por um motivo extra:** o §12 do nosso CLAUDE.md já cita **Soleimani 2017 §3.2.6**
como fonte do esquema de refinamento Vacondio/Feldman do Pass N. É o mesmo grupo — ou seja, o
projeto já usa a metodologia SPH dele sem ter registrado que ele a aplicou a biofilme. Vale fechar
essa citação.

**[T23] é o par direto da nossa `SurfactantEquation` + `MarangoniForce`** e o único trabalho que
encontrei fazendo surfactante com Marangoni em SPH de forma conservativa. Merece leitura antes de
qualquer mudança futura na forma do termo de produção/sumidouro.

## 5. Pares de partículas em escala de célula (para citar como contraste)

- **NUFEB** — Li, Jayathilake et al. (2019), *PLoS Comput. Biol.* 15:e1007125. IbM massivamente
  paralelo sobre LAMMPS, 10⁷ indivíduos, biofilme 3D com deformação e destacamento.
- **Jayathilake et al. (2017), *PLoS ONE* 12:e0181965** — IbM mecanicista de comunidades microbianas.
- **DPD de crescimento de biofilme** — arXiv:1808.00910.
- Todos resolvem **a célula**; a placa inteira lhes é inacessível por custo. É o argumento quantitativo
  para a nossa escolha de granularidade.

## 6. Como declarar a posição na tese

Três frases defensáveis, todas apoiadas no que está acima:

1. **Granularidade:** a partícula SPH representa ~10⁴-10⁵ bactérias, a mesma granularidade do
   *communicating walker* de Ben-Jacob (1994) — precedente estabelecido para modelar colônia sem
   resolver célula.
2. **Fenômeno:** o mecanismo (ramnolipídeo → gradiente de tensão superficial → Marangoni → tendrils
   em PA14) é o de Du et al. (2011) e Trinschek et al. (2018); a diferença é que aqui ele é resolvido
   por um método **lagrangiano de partículas na escala da colônia inteira**, e não por filme fino
   euleriano nem por agentes celulares numa fatia da borda.
3. **Lacuna:** não há, na literatura levantada, SPH aplicado a swarming em escala de colônia — o SPH
   bacteriano existente está em adesão de célula única, biofilme sob escoamento e transporte reativo
   em meio poroso. A rugosidade em escala de braço (Pass L) é extensão natural dessa lacuna.

## Fontes
- Ben-Jacob, Shochet, Tenenbaum, Cohen, Czirók & Vicsek, *Nature* 368:46-49 (1994). PMID 8107881.
- Kozlovsky, Cohen, Golding & Ben-Jacob, *Phys. Rev. E* 59:7025 (1999).
- Du, Xu, Shrout & Alber, *Math. Models Methods Appl. Sci.* 21(S1):939-954 (2011), doi:10.1142/S0218202511005428, PMC3182104.
- *A generalized model for predicting different morphologies of bacterial swarming on a porous solid surface*, PMID 39655366 (2025) — autores/periódico não confirmados.
- Kawasaki, Mochizuki, Matsushita, Umeda & Shigesada, *J. Theor. Biol.* (1997). PMID 9379672.
- Porter, Trenado-Yuste, Martinez-Calvo, Su, Wingreen, Datta & Huang, *Nat. Rev. Phys.* 7:535 (2025), doi:10.1038/s42254-025-00849-x.
- Soleimani, Wriggers, Rath et al., *Comput. Mech.* 58:619-633 (2016), doi:10.1007/s00466-016-1308-9.
- *Simulation of bacterial adhesion on a rough surface based on SPH*, *Comput. Part. Mech.* (2025), doi:10.1007/s40571-025-00924-1.
- Tartakovsky et al., *Pore-scale model for reactive transport and biomass growth*, *J. Porous Media* 12(5) (2009).
- *A conservative SPH method for surfactant dynamics* — referência completa a confirmar.
- Li, Jayathilake et al., *PLoS Comput. Biol.* 15:e1007125 (2019) — NUFEB. PMID 31830032.
