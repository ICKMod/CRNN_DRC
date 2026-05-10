# CRNN_DRC
Supplementary information/code for the publication at DOI: 

Julia workflows for kinetic network simulation and analysis used in the paper cases (`case1`, `case2`, `case3`).

## Contents

```text
Code/
├── case1/
│   ├── data.csv
│   └── workflow.jl
├── case2/
│   ├── data.csv
│   └── workflow.jl
├── case3/
│   ├── data.csv
│   └── workflow.jl
└── case3_reduced/
    ├── data.csv
    └── workflow.jl
```

Each `workflow.jl` is interactive and runs the full pipeline:
- profile generation (`Profile_ss.csv` or `Profile_ts.csv`)
- DRC analysis
- apparent activation energy (`E_app`)
- apparent reaction orders (`n_X`)

## Requirements

- Julia `1.6.7` (recommended)
- Packages used by scripts:
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

Install once in Julia REPL:

```julia
using Pkg
Pkg.add([
  "CSV","DataFrames","PrettyTables","DifferentialEquations",
  "DiffEqSensitivity","ForwardDiff","Plots","NLsolve","SteadyStateDiffEq"
])
```

## How to run

From inside any case directory:

```bash
cd case1
julia workflow.jl
```

or with explicit files:

```bash
julia workflow.jl data.csv CO2
```

Arguments:
- `ARGS[1]`: input CSV path (default: `data.csv`)
- `ARGS[2]`: target species for DRC/E_app/n_X (if omitted, script prompts)

## Interactive flow

1. Reads kinetics from `data.csv`.
2. Prompts for mode:
   - **Steady-state mode** (`ss` behavior):
     - computes steady state from initial guess
     - prompts for `steady_state_maxiters`
     - prompts for `steady_state_profile_points`
     - prompts for `steady_state_accept_resid`
   - **Transient mode** (`ts` behavior):
     - prompts for time span/save grid/solver maxiters
3. Prompts analysis grid (`analysis_start`, `analysis_end`, `analysis_points`).
4. Runs DRC, E_app, and n_X for the selected target.

## Output files

All outputs are written in the working case folder.

### Profiles
- `Profile_ss.csv` or `Profile_ts.csv`
- `Profiles_ss.png` or `Profiles_ts.png`

### DRC
- `DRC_Campbell_<target>.csv`
- `DRC_Campbell_<target>.png`
- `DRC_onesided_forward_<target>.csv`
- `DRC_onesided_forward_<target>.png`
- `DRC_onesided_reverse_<target>.csv`
- `DRC_onesided_reverse_<target>.png`

### Apparent activation energy
- `Eapp_<target>.csv`
- `Eapp_<target>.png`

`Eapp_<target>.png` includes:
- workflow `E_app` scatter
- Campbell-weighted line: `Σ(DRC_Campbell * Eaf)`
- one-sided weighted line: `Σ(DRC_onesided_forward * Eaf) + Σ(DRC_onesided_reverse * Ear)`

### Apparent reaction orders
For each gas species `X`:
- `nX_<X>_on_<target>.csv`
- `nX_<X>_on_<target>.png`

## Notes

- Scripts are designed for direct case folder execution so relative paths resolve correctly.
- If your terminal does not support interactive prompts, pass arguments and/or environment variables where available.
- For reproducible figures, keep the same `analysis_points`, tolerances, and solver settings across cases.
