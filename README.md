# ⚡ Quasi-Dynamic Load Flow for Distribution Systems

Este repositorio contiene scripts en Julia para el análisis de flujo de carga trifásico en sistemas de distribución eléctricos, con enfoque en simulaciones cuasi-dinámicas. Se incluyen herramientas para el modelado, solución iterativa y visualización de resultados en sistemas reales como redes europeas y rurales.

## 📁 Contenido del repositorio

| Archivo                               | Descripción                                                                 |
|--------------------------------------|-----------------------------------------------------------------------------|
| `FEEDER900.xlsx`                     | Sistema trifásico europeo de 906 nodos (media tensión, complejo y urbano). |
| `FEEDER14.xlsx`                      | Sistema rural alemán de media tensión con 14 nodos.                        |
| `load_feeder.jl`                     | Función para cargar los archivos Excel con los parámetros de red.         |
| `select_scenario.jl`                 | Selección de escenarios de carga en función del perfil y tiempo simulado. |
| `load_flow_newton.jl`                | Algoritmo iterativo basado en el método de Newton en el dominio complejo. |
| `load_flow_scl_fixed_point.jl`       | Solución de flujo con punto fijo (Self-Consistent Load Flow).             |
| `load_flow_sweep.jl`                 | Algoritmo tipo "Backward/Forward Sweep" para sistemas radiales.           |
| `load_flow_ybus.jl`                  | Solución de flujo usando la matriz de admitancia nodal (Ybus).            |
| `Quasi_dynamic_simulation_load_flow_scl_fixed_point.jl` | Simulación cuasi-dinámica con punto fijo.              |
| `Quasi_dynamic_simulation_load_flow_sweep.jl`          | Simulación cuasi-dinámica con algoritmo de barrido.      |
| `Quasi_dynamic_simulation_load_flow_ybus.jl`           | Simulación cuasi-dinámica usando Ybus.                   |

## ⚙️ Requisitos

- [Julia](https://julialang.org/) 1.6 o superior  
- Paquetes:
  - `LinearAlgebra`
  - `Printf`
  - `XLSX`
  - `Plots` (opcional para visualizaciones)

Puedes instalar paquetes ejecutando:

```julia
import Pkg
Pkg.add("MAT")
Pkg.add("XLSX")
Pkg.add("Plots")
