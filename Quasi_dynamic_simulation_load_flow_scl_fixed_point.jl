#|---------------------------------------------|
#| Power flow in power distribution networks   |
#| Active and reactive power in the substation |
#| Quasi dynamic simulation                    |
#| By:  Jhon Torres                            |
#|      jhon.torres1@utp.edu.co                |
#|---------------------------------------------|


"""Quasi-dynamic analysis of three-phase load flow on a 906-node European distribution feeder using complex linearization with the fixed-point method.
"""

# Import libraries and files
using LinearAlgebra
using Printf
using MAT
using Profile
include("load_feeder.jl")
include("select_scenario.jl")
include("funcion_graficas.jl")

file_name = "FEEDER900.xlsx"        # Feeder data file
feeder = load_feeder(file_name);    # Load feeder data from Excel

@time begin                         # Start timer to measure execution time

# Initialization of network parameters
num_n = feeder.num_n                # Number of nodes
ybus = feeder.ybus                  # Ybus admittance matrix
ynn = feeder.ynn                    # Ynn submatrix for non-slack nodes
znn = feeder.znn                    # Inverse 
yns = feeder.yns                    # Yns submatrix (non-slack to slack)
n_slack = feeder.n_slack            # Slack node indices
n_other = feeder.n_other            # Non-slack node indices
vn = feeder.vn_initial              # Initial voltages for non-slack nodes
vs = feeder.vs_initial              # Initial voltages for slack nodes
p_base = feeder.p_base              # Power base (MVA)
num_e = feeder.num_e                # Number of scenarios (typically 24h)



# Voltage initialization
global v = zeros(ComplexF64, 3*num_n, 10)  # Voltage matrix: rows = nodes, columns = iterations
v[n_other, 1] = vn                      # Set initial voltages for non-slack nodes
v[n_slack, 1] = vs                      # Set initial voltages for slack nodes

s_sub = ones(num_e, 3)                      # Matrix to store substation power per scenario

# Loop over all scenarios (e.g., 24 hours)
for tm = 1:num_e
    global loads = select_scenario(feeder, tm)  # Load scenario for current hour
    global err = 100                            # Initialize convergence error
    global iter = 1                             # Initialize iteration counter

    # Fixed-point iteration loop
    while err > 1e-9
        local s = v[:,iter].*conj(ybus*v[:, iter])  # Calculate nodal power injections

        local B = yns * vs + ynn * v[n_other, iter]       # Compute JB matrix
        local C = -conj(s[n_other] + loads[n_other])      # Net injected power at non-slack nodes
        local b = 1. ./ adjoint(v[n_other, iter])         # Inverse of diagonal voltage matrix
        local Inv_A = znn .* b                            # Approximate inverse of matrix A

        # Fixed-point proposal
        global dv = zeros(ComplexF64, 3 * (num_n - 1), 1) # Initialize voltage increment
        for k = 1:15
            dv_new = Inv_A * (C - (B .* conj(dv)))        # Update voltage increment
            if norm(dv_new - dv) < 1e-6                   # Check convergence
                break
            end
            global dv = dv_new
        end

        # Update voltages for next iteration
        v[n_other, iter + 1] = v[n_other, iter] + dv
        v[n_slack, iter + 1] = v[n_slack, iter]
        global err = norm(dv)                             # Update error
        global iter = iter + 1                            # Increase iteration count

        if iter > 9
            println("Complex linearization. After 10 iterations the error is :")
            break
        end
    end

    # Final voltage solution for the current scenario
    local v_n = ones(ComplexF64, 3 * num_n, 1)
    v_n[n_slack] = v[n_slack, iter]
    v_n[n_other] = v[n_other, iter]

    local s_n = v_n .* conj(ybus * v_n)                   # Calculate nodal powers
    s_sub[tm, 1] = tm                                     # Store time
    s_sub[tm, 2] = real(sum(real(s_n[n_slack]) * p_base * 1000))  # Real power in kW
    s_sub[tm, 3] = real(sum(imag(s_n[n_slack]) * p_base * 1000))  # Reactive power in kVAr

    # Reinitialize voltage matrix for next scenario
    global v = zeros(ComplexF64, 2718, 10)
    v[n_slack, 1] = v_n[n_slack]                          # Start from last solution
    v[n_other, 1] = v_n[n_other]
end

end

# Display substation real and reactive power for first 5 scenarios
display(s_sub[1:5, 2])
display(s_sub[1:5, 3])




















