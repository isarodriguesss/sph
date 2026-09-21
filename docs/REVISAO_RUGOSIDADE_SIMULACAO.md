# Revisão — como a literatura SIMULA superfície rugosa em biofilme e adesão bacteriana

> Levantada em 2026-09-18 para orientar o Pass L. Pergunta: a metodologia é 2D ou 3D, como os
> resultados são mostrados, e quais padrões de rugosidade são usados.

## 1. O quadro geral

| trabalho | método | dimensão | plano | padrão de rugosidade | como mostra |
|---|---|---|---|---|---|
| Rodriguez, Einarsson & Carpio, *Phys. Rev. E* 86, 061914 (2012) — biofilme em superfície rugosa | autômato celular estocástico | **2D** (versão 3D discutida, não usada) | **corte vertical**: duto retangular com fundo rugoso, escoamento ao longo de x | **degraus periódicos**: altura de pico ε = 2, comprimento λ = 2, distância entre picos δ = 5, em unidades de célula — "a rugosidade tem a mesma ordem de grandeza da bactéria" | sequências de instantâneos 2D de perfil |
| Jayathilake et al., CFD-DEM de twitching (bioRxiv 648915) | CFD + elementos discretos | **3D** | domínio 50 × 20 × 20 µm³ | **sulcos (grooves)**: paredes de 6 × 5 µm, largura do sulco 19 e 46 µm (2× e 4× o tamanho da bactéria) | renders 3D e trajetórias |
| *Simulation of bacterial adhesion on a rough surface based on SPH*, **Comput. Part. Mech.** (2025) | **SPH (PySPH!)** + potencial de Lennard-Jones para bactéria–superfície, WCSPH para fluido | 2D | corte vertical, bactéria rígida em escoamento | **topografia medida** por *shear force microscopy* — rugosidade real, não sintética | perfis e trajetórias da célula |
| *Modelling the combined effect of surface roughness and topography on bacterial attachment* (2021) | XDLVO + convecção-difusão-reação | 2D/3D conforme o caso | superfície vista de perfil | rugosidade **reconstruída de estatísticas de AFM** (Ra, Rq) somada a topografia projetada | mapas de energia de interação |
| *Microbial Response to Micrometer-Scale Multiaxial Wrinkled Surfaces*, ACS AMI (2022) | experimental (referência de escala) | — | — | **rugas senoidais**: λ ≈ 2 µm, amplitude 0.2 µm, Ra ≈ 0.07 µm | AFM 2D e 3D |
| Modelos individuais (iDynoMiCS e afins) | agentes + campos contínuos | 2D e 3D | corte vertical com cisalhamento no topo | substrato liso; rugosidade **emerge** da própria colônia, não é imposta | mapas de biomassa |

## 2. As quatro conclusões que interessam ao nosso desenho

**(a) 2D é a norma, e o motivo declarado é custo.** Carpio et al. dizem isso literalmente: "para reduzir
o custo computacional, focamos numa redução bidimensional", com a versão 3D discutida mas não usada.
O 3D aparece onde a física exige — escoamento cisalhante em torno de uma célula (CFD-DEM) — e mesmo lá
o domínio é minúsculo, 50 × 20 × 20 µm.

**(b) O plano 2D quase sempre é o CORTE VERTICAL, porque o fenômeno é biofilme submerso sob
escoamento.** O que se quer ver é o filme crescendo para cima e sendo erodido pelo fluxo. **O nosso
fenômeno é outro:** swarming se espalha NO plano do ágar, e as duas referências do projeto
(`reference.jpg` e Trinschek) são imagens de topo. Para nós, a redução 2D correta é a vista de topo —
o corte vertical mostraria a espessura do filme, que é justamente o que o nosso modelo não tem.

**(c) Os padrões usados são quatro, e todos são periódicos ou medidos:**
1. **degraus/cristas periódicas** com altura, comprimento e espaçamento (Carpio);
2. **sulcos** de largura fixa, calibrados em múltiplos do tamanho da célula — 2× e 4× (CFD-DEM);
3. **pilares** em rede, e poços de confinamento (literatura experimental de microtopografia);
4. **topografia medida** por AFM ou ShFM, reconstruída por estatísticas (Ra, Rq) — é o que o artigo de
   SPH de 2025 faz.

A regra comum a todos: **a escala da rugosidade é ancorada no tamanho da célula**, entre 1× e 4×.

**(d) Há um artigo de SPH com PySPH fazendo exatamente adesão em superfície rugosa (2025)** — mesma
biblioteca que usamos. Ele modela a bactéria como corpo rígido e a interação com a parede por
Lennard-Jones. É a referência metodológica mais próxima do nosso contato de pilar, e vale ler inteiro
antes de fechar a força de contato (hoje usamos Monaghan & Kajtar com kernel curto, ver
`docs/PLANO_RUGOSIDADE.md` seção 4).

## 3. O problema de escala do nosso caso — e ele é grande

A literatura ancora a rugosidade no tamanho da **célula** (µm). O nosso modelo não resolve célula: a
partícula SPH tem `dx` ≈ 0.054 em unidades do domínio, e ancorando pela largura de braço de PA14
([T11]: 2-5 mm) isso dá **`dx` ≈ 0.4-0.6 mm**. Ou seja:

- o pilar de 10 `dx` do nosso plano mede **~5 mm** — rugosidade macroscópica, da escala do braço;
- a rugosidade micrométrica da literatura (λ ~ 2 µm, Ra ~ 0.07 µm) está **três ordens de grandeza
  abaixo do nosso `dx`** e é irrepresentável aqui sem mudar de escala;
- o análogo correto na literatura não é AFM nem nanoestrutura: é **confinamento geométrico** (poços de
  31-90 µm que organizam enxames) e **microtopografia como barreira à motilidade de superfície**, que é
  o efeito descrito para *P. aeruginosa* em pilares.

**Consequência:** o que o nosso Pass L pode dizer é "obstáculo da escala do braço modula a morfologia
dendrítica", e não "rugosidade de superfície altera adesão". Isso precisa estar explícito na tese;
caso contrário a comparação com a literatura de adesão é uma troca de escala silenciosa.

## 3b. E em estudos de BIOFILME, qual padrão é adotado?

Pergunta levantada em 2026-09-18. A resposta muda conforme o experimento seja de **adesão/anti-fouling**
(padrão imposto pelo pesquisador) ou de **crescimento de biofilme sob escoamento** (rugosidade como
condição de contorno).

| família | padrão | escala típica | usado para |
|---|---|---|---|
| riblets em losango (**Sharklet AF**) | barras retangulares entrelaçadas | elementos de ~2 µm de largura e espaçamento, ~3 µm de altura | inibir adesão inicial; é o padrão sintético mais replicado |
| **pele de tubarão replicada** | dentículos dérmicos reais (riblets + nanoprotuberâncias nos sulcos) | dezenas a centenas de µm | anti-fouling "biomimético" |
| **pilares / poços** em rede | cilindros ou cavidades periódicas | 1-10 µm (adesão); 30-90 µm (confinamento de enxame) | barreira à motilidade de superfície, confinamento |
| **degraus/cristas periódicas** | perfil serrilhado no fundo do duto | da ordem do tamanho da célula | crescimento de biofilme sob cisalhamento (Carpio) |
| **topografia medida** (AFM/ShFM) | reconstrução por Ra/Rq | sub-µm | adesão, DLVO/XDLVO e o SPH de 2025 |

**Os dois artigos consultados pela usuária caem nas duas primeiras famílias.** Textos completos lidos
em 2026-09-18; os números abaixo são dos próprios artigos.

### Chien, Chen, Tsai & Lee (2020), *Colloids Surf. B* 186:110738 — pele de tubarão replicada

Pele de um **mako-de-barbatana-curta** macho de 140 cm (Pacífico, Taiwan), em quatro regiões: flanco
abdominal (A1), cauda abdominal (A2), nadadeira peitoral (F1) e nadadeira caudal (F2). Réplicas em
PMMA/PS/epóxi sobre PET, moldadas em PDMS.

| geometria do dentículo | abdômen | nadadeira |
|---|---|---|
| comprimento | 166-179 µm | 147-166 µm |
| largura | 86-99 µm | 65-70 µm |
| altura da crista | 10.5-13.5 µm | 6.2-8.8 µm |
| densidade | 119-125 /mm² | 84-131 /mm² |
| **Ra da réplica** | **6.08 µm (A2)** | **3.57 µm (F1)** |

Controle liso: **Ra = 0.166 µm**. Molhabilidade **anisotrópica** — ângulo de contato paralelo aos sulcos
105.1° (A2) e 100.8° (F1), perpendicular 115.8° e 105.6°, contra 73.8° no liso. Largura das cristas
42-55 µm.

**O resultado é mais interessante do que "rugoso inibe":** com *E. coli* (ATCC 23501), o padrão
**aumentou** a adesão no dia 1 (DO 0.45 e 0.39 contra 0.21 no liso) e mesmo assim **bloqueou o biofilme**
no dia 14 (DO 0.51 e 0.47 contra **3.22** no liso). Com *S. aureus* (ATCC 21351) não houve aumento inicial
(atribuído ao ar aprisionado na superfície hidrofóbica) e o dia 14 deu ~0.2 contra 2.07. Ou seja: a
topografia não impede a célula de pousar — **impede a microcolônia de se expandir**, atuando como
barreira física entre as reentrâncias e interferindo no quorum sensing.

### May et al. (2014), *Clin. Transl. Med.* 3:8 — Sharklet em tubo endotraqueal

Sharklet transferido por fotolitografia para wafer de silício e **gravado a 3 µm de profundidade** (deep
RIE), replicado em PDMS (Silastic T-2) em filmes de 0.4 mm, discos de 12 mm. **As larguras/espaçamentos
de ~2 µm não estão neste artigo** — ele remete à referência de fabricação; a fonte primária é
Schumacher/Chung (2007).

Organismos: MRSA, *K. pneumoniae*, *A. baumannii*, *E. coli* e ***P. aeruginosa*** em três cepas — dois
isolados clínicos (ATCC 9027, ATCC 10197) e a hiper-formadora de biofilme **PA14ΔbifA**, ou seja a mesma
linhagem PA14 do nosso projeto.

| ensaio | redução |
|---|---|
| colonização (5 espécies, CFU) | **95.6-99.9%** (p < 0.05) |
| biofilme MRSA, 4 d em TSB | 67% (p = 0.12, **não significativo**) |
| biofilme PA14ΔbifA, 24 h em meio mínimo | **52%** (p = 0.05) |
| biofilme *P. aeruginosa* ATCC 9027, meio com mucina | **58%** (p = 0.009) |

O meio "clínico" é mínimo + **2 mg/mL de mucina suína** + 400 µg/mL de oxacilina, 24 h estático a 37 °C —
montado porque a ATCC 9027 não formava biofilme no meio simples.

**Achado com consequência direta para a escolha do nosso padrão:** o artigo afirma (citando a literatura
de topografias ordenadas) que o Sharklet é **o mais eficaz entre as geometrias testadas — pilares,
canais e outras** — para inibir bio-adesão. É uma afirmação de adesão em microescala, não de swarming,
mas é a única comparação direta pilar × riblet que apareceu nesta revisão, e ela favorece o riblet.

**Leitura para o Pass L:** os dois são de **adesão/inibição em microescala**, não de morfologia de swarm,
e ambos estão três ordens de grandeza abaixo do nosso `dx` (§3) — a rugosidade da pele de tubarão (Ra 3.6
a 6.1 µm) vale ~0.01 `dx`. O que se aproveita é (a) a **geometria**: riblet/sulco domina nos dois e vence
pilares na única comparação direta, o que reforça a recomendação 3 abaixo; e (b) o **mecanismo**: em
Chien et al. a topografia não barra a célula individual, barra a **expansão da microcolônia** — que é
exatamente a classe de efeito que o nosso modelo (obstáculo da escala do braço modulando a frente
dendrítica) pode representar, ainda que em outra escala.

## 4. Recomendação para o projeto

1. **Manter a vista de topo** (item b). O corte vertical mudaria o fenômeno e invalidaria toda a
   métrica morfológica; o 3D é migração de esquema.
2. **Manter os pilares periódicos** como padrão principal — é o item 3 da lista, e a rede triangular
   com Λ = 28 `dx` e A = 10 `dx` já está ancorada na largura do braço, que é o análogo do "tamanho da
   célula" da literatura no NOSSO nível de descrição.
3. **Acrescentar depois um segundo padrão**, para não depender de uma geometria só: **sulcos**
   (faixas paralelas) na direção radial e transversal, que é o padrão com efeito mais documentado —
   guiam ou bloqueiam a motilidade conforme a largura.
4. **Terceiro padrão, se houver tempo:** campo aleatório com Ra/Rq controlados, à moda do que se faz
   reconstruindo AFM — no nosso caso, posições de pilar sorteadas com densidade fixa em vez de rede.
5. **Ler o artigo de SPH de 2025** antes de fechar a força de contato.

## Fontes
- Rodriguez, Einarsson & Carpio, *Biofilm growth on rugose surfaces*, Phys. Rev. E 86, 061914 (2012) — arXiv:2401.07135.
- Jayathilake, Li, Zuliani, Curtis & Chen, *Modelling bacterial twitching in fluid flows: a CFD-DEM approach*, bioRxiv 648915.
- *Simulation of bacterial adhesion on a rough surface based on smoothed particle hydrodynamics (SPH)*, Computational Particle Mechanics (2025), doi:10.1007/s40571-025-00924-1.
- *Modelling the combined effect of surface roughness and topography on bacterial attachment*, J. Mater. Sci. Technol. (2021), doi:10.1016/j.jmst.2020.12.009.
- *Microbial Response to Micrometer-Scale Multiaxial Wrinkled Surfaces*, ACS Appl. Mater. Interfaces (2022), PMC9284519.
- Chien, Chen, Tsai & Lee, *Inhibition of biofilm formation by rough shark skin-patterned surfaces*, Colloids Surf. B: Biointerfaces 186 (2020) 110738, doi:10.1016/j.colsurfb.2019.110738, PMID 31869602. **Texto completo lido.**
- May et al., *Micro-patterned surfaces reduce bacterial colonization and biofilm formation in vitro*, Clin. Transl. Med. 3:8 (2014), doi:10.1186/2001-1326-3-8, PMC3996152. **Texto completo lido.**
- Schumacher, Chung, Brennan et al. (2007) — dimensões canônicas do Sharklet AF (largura/espaçamento ~2 µm). **Não lida**; é a fonte primária do padrão, citada por May et al. para a fabricação.
- *Confinement discerns swarmers from planktonic bacteria*, eLife (2021), e64176.
- *Surface Hardness Impairment of Quorum Sensing and Swarming for P. aeruginosa*, PLOS One (2011), e20888.
