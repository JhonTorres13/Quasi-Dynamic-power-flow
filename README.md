# ⚡ Quasi-Dynamic Load Flow for Distribution Systems

## 🧪 General Description

This simulation environment allows the analysis of the electrical behavior of three-phase distribution systems using load flow algorithms and quasi-dynamic simulations. The models enable:

- Representation of real-world networks, including urban European systems and German rural grids.
- Evaluation of various numerical methods, including:
  - Newton's method in the real domain.
  - Newton's method in the complex domain using Wirtinger calculus and a fixed-point method to solve the non-holomorphic complex linear system (Proposed Method).
  - Backward/Forward sweep algorithm (for radial systems).
  - Direct fixed-point method.
- Hourly behavior simulation of the grid (24-hour horizon) with variable load profiles.
- Comparison of accuracy, stability, and computational performance among the methods.

This set of tools is aimed at research and analysis of modern and active electrical distribution systems.

---

## 📁 Repository Content

| File                                      | Description                                                                 |
|------------------------------------------|-----------------------------------------------------------------------------|
| `FEEDER900.xlsx`                         | Data of the 906-node European three-phase system (complex urban network).  |
| `FEEDER14.xlsx`                          | Data of the 14-node German rural medium-voltage meshed system.             |
| `load_feeder.jl`                         | Function to load Excel files with network parameters.                      |
| `select_scenario.jl`                     | Scenario selector based on the load profile and simulation setup.          |
| `load_flow_newton.jl`                    | Load flow solution using Newton's method in the real domain.              |
| `load_flow_scl_fixed_point.jl`           | Load flow using Newton's method in the complex domain with Wirtinger calculus and fixed-point solver (proposed method). |
| `load_flow_sweep.jl`                     | Load flow solution using the Backward/Forward Sweep algorithm for radial networks. |
| `load_flow_ybus.jl`                      | Load flow using the direct fixed-point method.                             |
| `Quasi_dynamic_simulation_load_flow_scl_fixed_point.jl` | Quasi-dynamic simulation using Newton's method in the complex domain with Wirtinger calculus and fixed-point method. |
| `Quasi_dynamic_simulation_load_flow_sweep.jl`          | Quasi-dynamic simulation using the Backward/Forward Sweep method.          |
| `Quasi_dynamic_simulation_load_flow_ybus.jl`           | Quasi-dynamic simulation using the direct fixed-point method.              |

---

## ⚙️ Requirements

- [Julia](https://julialang.org/) version 1.6 or higher  
- Packages:
  - `LinearAlgebra`
  - `Printf`
  - `XLSX`
  - `Plots` (optional for visualization)

To install the required packages:

```julia
import Pkg
Pkg.add("MAT")
Pkg.add("XLSX")
Pkg.add("Plots")
```


1. 🚀 Execution:
Load feeder data:

```julia
include("load_feeder.jl")
feeder = load_feeder("FEEDER900.xlsx")  # o "FEEDER14.xlsx"
```

2. Select the load scenario:

```julia
include("select_scenario.jl")
cargas = select_scenario(feeder, tiempo)
```

3. Run the desired load flow method:

```julia
load_flow_scl_fixed_point.jl
# load_flow_newton.jl
# load_flow_sweep.jl
# load_flow_ybus.jl
```
4. Run the 24-hour quasi-dynamic simulation:

```julia
Quasi_dynamic_simulation_load_flow_scl_fixed_point.jl
# Quasi_dynamic_simulation_load_flow_sweep.jl
# Quasi_dynamic_simulation_load_flow_ybus.jl
```

## 🧑‍🔬 Author
Jhon A Torres
Electrical Engineer - M.Sc. in Electrical Engineering (in progress).
