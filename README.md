# A Unified Workflow for Sensitivity Based Kinetic Analysis in Microkinetic Models

Code for the publication: **[A Unified Workflow for Sensitivity-Based Kinetic Analysis in Microkinetic Models](https://doi.org/10.1021/acs.iecr.6c01462)**

Julia workflows for kinetic-network simulation and local sensitivity analysis used in the paper cases (`case1`, `case2`, `case3`). Each case folder contains a reaction-network input file, usually `data.csv`, and a copy of `workflow.jl`.

The script reads an elementary reaction table, builds a stoichiometric kinetic model, solves either transient or steady-state profiles, and computes local sensitivity descriptors for a selected target species.

## Contents

```text
Code/
├── case1/
│   ├── data.csv
│   └── workflow.jl
├── case2/
│   ├── data.csv
│   └── workflow.jl
└── case3/
    ├── data.csv
    └── workflow.jl
```

## Main capabilities

Each `workflow.jl` runs the following pipeline:

1. Reads the reaction table from CSV.
2. Auto-detects the reaction column.
3. Parses reversible and irreversible elementary steps.
4. Builds a species catalog classified as:
   - free site: `*`
   - adsorbed species: labels ending in `*`
   - gas species: all other species
5. Builds the stoichiometric matrix for non-gas species.
6. Prompts for gas pressures.
7. Solves either transient profiles or steady-state profiles.
8. Computes, for a selected target species:
   - Campbell DRC
   - one-sided forward DRC
   - one-sided reverse DRC
   - apparent activation energy, if Arrhenius columns are available
   - one-sided barrier-weighted `E_app`
   - state-enthalpy weighted `E_app`, when possible
   - apparent reaction order with respect to each gas species

The target rate is evaluated as the net stoichiometric formation rate of the selected target species from the direction-expanded reaction network.

## Requirements

Recommended Julia version depends on the local package environment. The script uses the following packages:

- `CSV`
- `DataFrames`
- `PrettyTables`
- `Printf`
- `DifferentialEquations`
- `DiffEqSensitivity`
- `ForwardDiff`
- `Plots`
- `Markdown`
- `NLsolve`
- `SteadyStateDiffEq`
- `LinearAlgebra`

`LinearAlgebra`, `Printf`, and `Markdown` are Julia standard-library packages. The others can be installed from the Julia REPL:

```julia
using Pkg
Pkg.add([
    "CSV",
    "DataFrames",
    "PrettyTables",
    "DifferentialEquations",
    "DiffEqSensitivity",
    "ForwardDiff",
    "Plots",
    "NLsolve",
    "SteadyStateDiffEq"
])
```

The script calls `pyplot()` through `Plots.jl`. If the plotting backend is not already available in your environment, install `PyPlot` as well:

```julia
using Pkg
Pkg.add("PyPlot")
```

## How to run

From inside a case directory:

```bash
julia workflow.jl
```

With an explicit input CSV:

```bash
julia workflow.jl data.csv
```

With an explicit input CSV and target species:

```bash
julia workflow.jl data.csv CO2
```

Arguments:

- `ARGS[1]`: input CSV path. Default: `data.csv`.
- `ARGS[2]`: target species for DRC, `E_app`, and reaction-order analysis. If omitted, the script prompts for it.

## Input CSV requirements

The input file must contain one column with elementary reaction strings. The script first looks for common names such as `Elementary_Step` or `Step`; otherwise it auto-detects the column containing reaction arrows.

Supported reaction arrows include:

```text
⇌, ↔, <->, <=>, ⇄, ⟷, ->, →, ←
```

The script normalizes reversible arrows to `⇌` for display.

### Kinetic input modes

The script auto-detects whether the CSV contains direct rate constants or Arrhenius parameters.

Supported pair types:

1. `(Keq, kf)`
   - forward rate constant and equilibrium constant
   - reverse rate is computed as `kr = kf / Keq`

2. `(kf, kr)`
   - forward and reverse rate constants supplied directly
   - equilibrium constant is computed as `Keq = kf / kr`

Supported direct-rate columns include names similar to:

```text
k, kf, k_f, kforward, kForward
kr, k_r, kreverse, kReverse
K, Keq, K_eq, KEquilibrium, Kconst, Keqconst
```

Supported Arrhenius columns include names similar to:

```text
Af, A_f, Aforward, AForward, A_fwd, A_F
Eaf, Ea_f, EaForward, Ea_forward, Eafwd, Ea_F
Ar, A_r, Areverse, AReverse, A_rev, A_R
Ear, Ea_r, EaReverse, Ea_reverse, Earv, Ea_R
```

For Arrhenius input, the script computes:

```text
k(T) = A T^b exp[-Ea / (R T)]
```

The activation-energy unit can be `kJ/mol` or `J/mol`.

### Thermodynamic consistency check

If `K`, `kf`, and `kr` are all available for a reversible row, the script checks whether:

```text
K ≈ kf / kr
```

If the mismatch exceeds the log tolerance, the script stops with an error. The default tolerance is:

```text
DRC_THERMO_CHECK_LOGTOL=1e-6
```

This is a local consistency check for the supplied kinetic inputs. It does not force a globally consistent thermodynamic cycle across the entire network.

## Interactive flow

The script asks whether to run steady-state mode:

```text
Go to steady state? (y/n)
```

### Transient mode

In transient mode, the script prompts for:

- initial coverages for all non-gas species
- ODE time span
- profile save grid
- ODE solver
- `maxiters`
- analysis time grid
- target species, if not provided through `ARGS[2]` or `DRC_TARGET`

### Steady-state mode

In steady-state mode, the script starts from a clean-surface guess when `*` is present. It uses dynamic steady-state termination and prompts for:

- `steady_state_maxiters`
- `steady_state_profile_points`
- `steady_state_accept_threshold`

If steady-state convergence fails, the script switches to transient mode and asks for initial coverages manually.

## ODE solver settings

The default solver is `Rodas5P`. Other solver names accepted by the script include:

```text
Kvaerno5, Tsit5, Vern9, RK4, Rodas5, Rodas5P, TRBDF2, Rodas4P
```

The default ODE maximum iteration count is `200000`, unless overridden.

The ODE solve uses tight tolerances:

```text
abstol = 1e-12
reltol = 1e-12
dtmin  = 1e-25
```

## Environment variables

Most prompts can be controlled using environment variables. Useful variables include:

| Variable | Purpose |
|---|---|
| `DRC_STEADY_STATE` | Use steady-state mode if set to `y`, `yes`, `true`, or `1`. |
| `DRC_SOURCE` | Force kinetic source mode: `direct` or `arrhenius`. |
| `DRC_PAIR` | Force kinetic pair type: `(keq,kf)` or `(kf,kr)`. |
| `DRC_TARGET` | Target species for analysis. |
| `DRC_EA_UNIT` | Activation-energy unit: `kJ/mol` or `J/mol`. |
| `DRC_B_ARR` | Arrhenius temperature exponent `b`. |
| `DRC_ALG` | ODE solver name. |
| `DRC_MAXITERS` | Default ODE `maxiters`. |
| `DRC_THERMO_CHECK_LOGTOL` | Log tolerance for `K ≈ kf/kr` consistency check. |
| `DRC_ENTHALPY_REFERENCES` | Reference states for state-enthalpy decomposition, for example `*=0,O2=0,C3H6=0`. |

Prompt-specific variables are also generated from prompt names by converting spaces to underscores and using uppercase. For example:

```text
analysis_start      -> ANALYSIS_START
analysis_end        -> ANALYSIS_END
analysis_points     -> ANALYSIS_POINTS
steady_state_maxiters -> STEADY_STATE_MAXITERS
```

## Output files

All outputs are written to the current working directory.

### Profiles

For steady-state mode:

```text
Profile_ss.csv
Profiles_ss.png
```

For transient mode:

```text
Profile_ts.csv
Profiles_ts.png
```

The profile CSV contains `Time_Steps` and one column for each non-gas species.

### DRC outputs

For target species `<target>`:

```text
DRC_Campbell_<target>.csv
DRC_Campbell_<target>.png
DRC_onesided_forward_<target>.csv
DRC_onesided_forward_<target>.png
DRC_onesided_reverse_<target>.csv
DRC_onesided_reverse_<target>.png
```

The Campbell DRC perturbs each reversible elementary step through the paired forward-rate/equilibrium-constant representation. The one-sided DRC files report independent sensitivity to the forward and reverse directional rate constants.

### Apparent activation energy outputs

If Arrhenius columns are available, the script writes:

```text
Eapp_<target>.csv
Eapp_<target>.png
```

The CSV contains:

```text
Time
E_app_workflow
E_app_onesided_barrier_weighted
E_app_state_enthalpy_weighted
E_app_state_minus_workflow
E_app_onesided_minus_workflow
max_state_enthalpy_residual
```

The plot shows:

- workflow-computed `E_app` as a solid line
- one-sided barrier-weighted `E_app` as scatter points
- state-enthalpy weighted `E_app` as scatter points, when the state reconstruction succeeds

The old direct Campbell DRC multiplied by activation barriers is intentionally not used for the final `E_app` decomposition.

### State-enthalpy decomposition outputs

When Arrhenius data and reversible reactions are available, the script attempts a Campbell-style state-enthalpy decomposition. It writes:

```text
State_enthalpies_<target>.csv
State_enthalpy_consistency_<target>.csv
State_DRC_stable_<target>.csv
State_DRC_TS_<target>.csv
```

Reference enthalpies must be supplied interactively or through `DRC_ENTHALPY_REFERENCES`. Example:

```bash
export DRC_ENTHALPY_REFERENCES="*=0,O2=0,C3H6=0"
```

If the enthalpy reconstruction is underdetermined, ill-conditioned, or thermodynamically inconsistent, the script prints warnings and continues with a least-squares or minimum-norm reconstruction. In that case, the result should be interpreted as a directed state-energy regrouping rather than a strictly consistent thermodynamic diagram.

### Apparent reaction-order outputs

For each gas species `X` and selected target `<target>`:

```text
nX_<X>_on_<target>.csv
nX_<X>_on_<target>.png
```

The CSV contains:

```text
Time
n_X
```

The reaction order is computed by differentiating the target rate with respect to the logarithm of the gas pressure.

## Running multiple targets

After completing analysis for one target, the script asks:

```text
Run another target? Press ESC then Enter to exit; Enter to continue
```

Press Enter to analyze another target species without restarting the script. Type `n`, `no`, `q`, `quit`, or `exit`, or press Escape then Enter, to stop.

## Notes

- For difficult stiff systems, increase `DRC_MAXITERS` or choose a stiff solver such as `Rodas5P`, `Rodas5`, `TRBDF2`, or `Kvaerno5`.
- When using Arrhenius data, keep the activation-energy unit consistent with the CSV file.
- When using state-enthalpy decomposition, choose reference states that match the thermodynamic convention used in the manuscript.

# Cite

Jay Shukla, Qin Wu; A Unified Workflow for Sensitivity-Based Kinetic Analysis in Microkinetic Models. Ind. Eng. Chem. Res. 15 July 2026; 65 (27): 14161–14176., (https://doi.org/10.1021/acs.iecr.6c01462)


