# Handoff — Implementação do Pass M-A (Osmólito c_o)

## Contexto do Projeto

Simulação SPH 2D de *P. aeruginosa* swarming (tese Isadora Rodrigues, UFAL/NEES).
Objetivo: reproduzir morfologia dendrítica de `reference.jpg` (15-20 dendritos, AR ≥ 1:5).

**Estado atual (branch `develop`):**
- K.17: hard pinning (`u=v=0` em `ρ_b ≥ 0.8`) ✅
- K.20: `k_consume=2.0` — Bloqueio B (saturação química cs) resolvido ✅
- K.24: `r_growth=0.05` — K.17 Trap resolvido ✅
- K.26: `f0=1.5`, `threshold=0.8` — 8 dendritos persistentes ✅
- K.27: `f0=3.0` — dendritos com AR ~1:6, **mas halo radial de baia persiste** ❌

**Problema diagnosticado (K.27):** O modelo SPH é um sistema de massa fechada. Todas as partículas da borda recebem empurrão radial uniforme (Marangoni + Flagelar). O diferencial motile_boost tip/baia (~2.3×) é insuficiente para suprimir baias. **Calibração paramétrica esgotada.**

---

## Tarefa: Implementar Pass M-A (Campo de Osmólito `c_o`)

**Mecanismo físico alvo (Srinivasan 2019, Bru 2023):**
Bactérias secretam osmólitos → acumulam no agar das baias (rodeadas por colônia dos dois lados) → agar virgem nas pontas tem c_o ≈ 0 → gradiente |∇c_o| grande nas pontas, pequeno nas baias → influxo de fluido do agar assimétrico → pontas incham, baias estagnám → dendritos finos com agar limpo.

---

## Arquivos a Modificar

### 1. `src/particles.py` — Adicionar propriedade `c_o`

Dentro do bloco `if pa.name == "fluid":`, após as propriedades existentes, adicionar:
```python
pa.add_property("c_o")        # osmólito produzido pelas bactérias
pa.add_property("a_c_o")      # acumulador (derivada)
pa.c_o[:] = 0.0               # agar começa sem osmólito
pa.add_output_arrays(["c_o"]) # para visualização
```

### 2. `src/equations.py` — Nova equação `OsmolyteProduction`

Adicionar após `class OsmoticForce` (linha ~363), antes de `class FlagellarForce`:

```python
class OsmolyteProduction(Equation):
    """
    Campo de osmólito c_o produzido pelas bactérias.
    Satura nas baias (agar confinado entre dois dendritos).
    Fresco nas pontas (agar virgem).
    Influxo V0 ∝ |∇c_o| adiciona massa às pontas, não às baias.
    Ref: Srinivasan et al. 2019 (eLife), Bru et al. 2023 (Biophys. Rev.)
    """

    def __init__(self, dest, sources, D_o=1e-3, k_o=0.5, lambda_o=0.05, Q0=5e-4):
        self.D_o = D_o
        self.k_o = k_o
        self.lambda_o = lambda_o
        self.Q0 = Q0
        super(OsmolyteProduction, self).__init__(dest, sources)

    def initialize(self, d_idx, d_a_c_o, d_grad_co_x, d_grad_co_y):
        d_a_c_o[d_idx] = 0.0
        d_grad_co_x[d_idx] = 0.0
        d_grad_co_y[d_idx] = 0.0

    def loop(
        self,
        d_idx, s_idx,
        d_c_o, s_c_o,
        d_a_c_o,
        d_grad_co_x, d_grad_co_y,
        d_rho, s_rho,
        s_m,
        DWIJ, WIJ,
        XIJ, RIJ,
        d_h,
    ):
        # Gradiente SPH simétrico de c_o
        co_ij = s_c_o[s_idx] - d_c_o[d_idx]
        Vj = s_m[s_idx] / s_rho[s_idx]

        d_grad_co_x[d_idx] += Vj * co_ij * DWIJ[0]
        d_grad_co_y[d_idx] += Vj * co_ij * DWIJ[1]

        # Difusão Brookshaw
        if RIJ > 1e-12:
            eij_dot_dwij = (XIJ[0] * DWIJ[0] + XIJ[1] * DWIJ[1]) / (RIJ * RIJ)
            d_a_c_o[d_idx] += (
                2.0 * self.D_o * Vj * co_ij * eij_dot_dwij
            )

    def post_loop(
        self,
        d_idx,
        d_a_c_o,
        d_c_o,
        d_rho_b_grown,
        d_grad_co_x, d_grad_co_y,
        d_au, d_av,
        d_m,
    ):
        rho_b = d_rho_b_grown[d_idx]

        # Produção: bactérias secretam osmólito (proporcional à densidade)
        qs = rho_b * rho_b / (rho_b * rho_b + 0.01)
        production = self.k_o * qs * (1.0 - d_c_o[d_idx])  # satura em c_o=1

        # Decaimento lento (osmólitos persistem)
        decay = self.lambda_o * d_c_o[d_idx]

        d_a_c_o[d_idx] += production - decay

        # Influxo osmótico: V0 ∝ |∇c_o|, só em swarmers de borda
        # gate: smoothstep em rho_b ∈ [0.1, 0.6]
        if rho_b < 0.1:
            gate = 0.0
        elif rho_b > 0.6:
            gate = 0.0
        else:
            t = (rho_b - 0.1) / (0.6 - 0.1)
            gate = t * t * (3.0 - 2.0 * t)
            # pico em rho_b=0.35
            t2 = (rho_b - 0.35) / 0.25
            gate *= max(0.0, 1.0 - t2 * t2)

        grad_mag = (
            d_grad_co_x[d_idx] * d_grad_co_x[d_idx]
            + d_grad_co_y[d_idx] * d_grad_co_y[d_idx]
        ) ** 0.5

        # Influxo aumenta massa da partícula
        # (equivalente a V0 = Q0 * |∇c_o| * gate entrando na partícula)
        d_m[d_idx] += self.Q0 * gate * grad_mag
```

Adicionar também `grad_co_x` e `grad_co_y` às propriedades em `src/particles.py`:
```python
pa.add_property("grad_co_x")
pa.add_property("grad_co_y")
```

### 3. `src/scheme.py` — Integrar `c_o` no `CustomEulerStep` e importar/usar equação

**Em `CustomEulerStep.stage1`**, adicionar parâmetros e integração:
```python
# Adicionar à assinatura do método:
d_c_o,
d_a_c_o,

# Adicionar ao corpo, junto com rho_b_grown e cs:
d_c_o[d_idx] += dt * d_a_c_o[d_idx]
d_c_o[d_idx] = max(0.0, min(d_c_o[d_idx], 1.0))
```

**No import**, adicionar `OsmolyteProduction` ao import de equations:
```python
from .equations import (
    ...
    OsmolyteProduction,
)
```

**No `MyBiomassScheme.__init__`**, adicionar parâmetros:
```python
D_o=1e-3,
k_o=0.5,
lambda_o=0.05,
Q0=5e-4,
```
e `self.D_o = D_o` etc.

**No `get_equations()`**, adicionar `OsmolyteProduction` ao `equations_main` **antes** de `FlagellarForce`:
```python
OsmolyteProduction(
    dest="fluid",
    sources=["fluid"],
    D_o=self.D_o,
    k_o=self.k_o,
    lambda_o=self.lambda_o,
    Q0=self.Q0,
),
```

### 4. `main.py` — Adicionar parâmetros e passar ao scheme

```python
# Parâmetros Pass M-A
D_o = 1e-3       # difusão lenta (osmólito alto peso molecular)
k_o = 0.5        # taxa de produção
lambda_o = 0.05  # decaimento lento (osmólitos persistem)
Q0 = 5e-4        # acoplamento influxo-massa
```

No `create_scheme()`, passar os novos parâmetros ao `MyBiomassScheme`.

---

## Critérios de Sucesso Pass M-A

1. Campo `c_o` no interior da colônia satura em ~0.8-1.0 em t > 30s
2. Campo `c_o` no agar das baias maior que no agar das pontas (visível nos frames)
3. `n_fast` concentrado nas pontas (não distribuído em todo o rim)
4. Halo radial de baia desaparece nos frames
5. Dendritos finos com agar limpo entre eles — AR ≥ 1:5
6. `mean_cs` plateau em [0.15, 0.45], `contrast_cs` > 100 sustentado (manter Pass J/K.20 intactos)

---

## Invariantes Críticos a Preservar

- `BiomassGradient` deve executar ANTES de `MarangoniForce`
- `SurfactantEquation` deve executar ANTES de `FlagellarForce`
- `OsmolyteProduction` deve executar ANTES de `FlagellarForce` (usa `c_o` atualizado)
- Hard pinning (`u=v=0` em `ρ_b ≥ 0.8`) em `CustomEulerStep` — NÃO REMOVER
- `lambda_ext_ratio=5.0` em `SurfactantEquation` — NÃO ALTERAR (Pass J)
- `k_consume=2.0` — NÃO ALTERAR (Bloqueio B)

---

## Parâmetros Físicos Atuais (não alterar nesta implementação)

```python
mu = 0.020
gamma = 60.0
beta = 1.0
sigma = 1.2
D = 1.5e-3
D_ext = 0.08
lambda_ = 0.15
r_growth = 0.05
c0 = 0.35
alpha_mon = 0.12
f0 = 3.0  # FlagellarForce em scheme.py
```

---

## Para Executar

```bash
conda activate pysph_env
make run
```

Frames gerados em `main_output/movie/`. Log em `log.csv`.
