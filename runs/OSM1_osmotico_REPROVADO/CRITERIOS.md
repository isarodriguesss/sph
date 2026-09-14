# OSM1 — influxo osmotico (absorcao de agar) no lugar da deposicao

Pre-registrado em 2026-09-11, ANTES de lancar.

## Config
- Base: E11 (`git fe49aba`), reproducao conferida bit-a-bit nas iteracoes 0 e 200
  (`runs/OSM_repro` contra `runs/E11_t100/log.csv`, inclusive a sequencia do dt).
- OSM1: `use_wake=False`, `use_insert=False`, `use_osm=True`, `OSM_K=32`, `OSM_PHI=0.4`,
  `OSM_FREQ=20`, t=50.
- OSM0 (controle): o MESMO codigo com `use_osm=False`. Isola deposicao de absorcao.

## Mecanismo
Agar vizinho de VIVA (`is_filler<0.5`, `rho_b>0.1`) absorvido com probabilidade
`1-exp(-OSM_K*D*dt)`, `D = sum_j V_j rho_b_j W_ij` sobre vivas. Absorvido: `phi_osm=0.4`,
`rho_b` continua 0 -> quimica de agar (nao cresce, nao produz cs, nao consome c_n,
colonizavel); a `BiomassEOS` le `max(rho_b, phi_osm)` -> pressao/coesao (fade 0.84).

Medido no E11 (t=20-50): D p50 ~0.004 e p90 ~0.025-0.045 na frente; com K=32 a taxa na
ponta e ~1/s e a mediana ~0.1/s. Alcance maximo 2h de uma viva (nao ha flood-fill).

## Predicoes (§2.3)
1. **Dinamica OSM1 ~ OSM0.** A EOS do absorvido e >=15x mais fraca que a Marangoni (a_p
   ~0.01 no anel comprimido contra a_mar 0.2-3; licao #85). R99, dR/dt, vivas, amplitude
   e dedos devem ficar a <=10% entre OSM1 e OSM0.
2. **C5 do OSM1 alto POR CONSTRUCAO** — o absorvido e o rastro das vivas dilatado ~2h, e
   as vivas saem do inoculo por trajetorias continuas. Por isso C5 SOZINHO nao conta como
   evidencia; so conta junto dos criterios de forma.
3. **Vazio (C1) pior que o E11** nos dois (sem wake, o lugar que a viva desocupa nao e
   preenchido e o agar nao cicatriza, licao #40); OSM1 ~ OSM0.
4. OSM0 com mais dedos e menos R99 que o E11 (mesmo sinal do I4 contra o I3, licao #85-F).

## Criterios de aprovacao (OSM1 contra E11; forma em R99 IGUAL; C5 em media t in [35,50])
- C5b (colonia = rho_b>=0.1 ou filler ou phi_osm>0, ligacao 1.05 dx) >= 50% e acima do
  E11 por mais de 2 desvios; reportar tambem a coluna +limbo.
- Amplitude >= 0.8x e dedos >= 0.8x o E11.
- Agar limpo na baia >= o do E11.
- R99(t=50) >= 0.8x o E11 (3.95 -> >= 3.16).
- Massa sem runaway; a_pressure mediana <= 4.
- Leitura visual (§11) da usuaria.

## Atribuicao
Se OSM1 passar mas OSM0 tiver a mesma forma, o que a absorcao acrescenta e so a
CLASSIFICACAO (continuidade do conjunto), sem efeito dinamico — isso sera dito assim.
