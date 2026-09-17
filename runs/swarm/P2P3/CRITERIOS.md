# P2P3 — P2 + gate flagelar [0.1, 0.8] (P3)

Pre-registrado em 2026-09-11, antes de lancar. Controle: `runs/P2_fillerdonor` (bit-a-bit
com o codigo atual em `FLAG_GATE` = [0.2, 0.6], verificado em `runs/_verif_P2`). Alavanca
unica: `FLAG_GATE_LO/HI` 0.2/0.6 -> 0.1/0.8. t=50.

## Hipotese
O P3 sobre o E5 deu o melhor relevo da classe dendritica (amplitude 0.412 vs 0.345, 40
dedos, fundo das baias 1.08) e piorou a continuidade (C5a 0.24 -> 0.21), porque forca
constante numa banda mais larga cria cavaleiros (p90/p50 de deslocamento 5.0 -> 6.8). O P2
e o que liga os bracos ao nucleo. A pergunta e se as duas coisas somam.

## Predicao
- No P2 o gate alargado engaja **4.3x** mais vivas em t=34-50 (388 vs 93; 1162 vs 271),
  contra 2-3x no E11 — o P2 tem 4x mais vivas, a maioria recrutas em rho_b 0.1-0.2, onde
  o smoothstep do gate ainda e pequeno (0.05 em 0.15). O efeito deve ser MAIOR que no P3.
- Relevo sobe para perto do P3 (amplitude 0.34-0.40 em R99=2.4); dedos 38-40.
- Risco principal: os cavaleiros saem do corpo e voltam a formar ilhas — bracos soltos
  sobem de 1.3 para 2-4.

## Criterios (contra o P2, janela t in [35, 50]; forma em R99 = 2.4)
- **Bracos soltos** (figura renderizada, > 20 dx^2) <= 2.0 e area fora do corpo <= 3%.
- **Relevo**: amplitude >= 0.34 **ou** fundo das baias (Rmin em t=50) <= 1.5.
- **C5a** >= 0.40.
- Agar limpo nas baias >= 75%; a_pressure mediana <= 3 e picos > 4 <= 10%.
- Classe morfologica dendritica na leitura visual (§11).
- Reportar tambem a razao p90/p50 de deslocamento (cavaleiros).
