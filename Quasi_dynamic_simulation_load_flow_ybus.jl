#|---------------------------------------------|
#| Power flow in power distribution networks   |
#| Active and reactive power in the substation |
#| Quasi dynamic simulation                    |
#| By:  Jhon Torres                            |
#|      jhon.torres1@utp.edu.co                |
#|---------------------------------------------|


"""Three-phase load flow using the fixed-point method directly"""

# Libraries
using LinearAlgebra 
using Printf
using MAT
using Profile

# Include custom functions
include("load_feeder.jl")           # Function to load feeder data
include("select_scenario.jl")       # Function to select time-based load scenario

# Load feeder data from Excel file
nombre_ar = "FEEDER900.xlsx"        # Specify feeder data file
feeder = load_feeder(nombre_ar);    # Load feeder data into struct

# Extract network parameters from feeder data
num_e = feeder.num_e                # Number of load scenarios
num_n = feeder.num_n                # Number of nodes
ybus = feeder.ybus                  # Full system admittance matrix
ynn = feeder.ynn                    # Reduced admittance matrix (non-slack nodes)
znn = feeder.znn                    # Inverse of Ynn (optional use)
yns = feeder.yns                    # Coupling admittance between slack and other nodes
n_slack = feeder.n_slack            # Indices of slack bus nodes
n_other = feeder.n_other            # Indices of all other (non-slack) nodes
vn = feeder.vn_initial              # Initial voltage guess at non-slack nodes
vs = feeder.vs_initial              # Fixed voltages at slack nodes
p_base = feeder.p_base              # Power base (used to convert to kW and kVAr)

# Start timing the algorithm
@time begin

    s_sub = ones(5, 3)                                     # Matrix to store results for each scenario
    global vnn = zeros(ComplexF64, 2715, 10)                   # Voltage history matrix for iteration
    vnn[:, 1] = vn                                             # Set initial voltage guess

    for tm = 1:5                                            # Loop over all time scenarios
        local cargas = select_scenario(feeder, tm)              # Get complex load for current time
        global err = 100;                                       # Initial error
        global iter = 1;                                        # Iteration counter

        while err > 1E-8                                        # Iterative loop (Newton-like)
            kk = norm(vnn[:, iter])                            # Store norm of current voltage estimate
            # Update voltage at non-slack nodes using complex power injection model
            vnn[:, iter + 1] = ynn \ (conj(-cargas[n_other] ./ vnn[:, iter]) - yns * vs)
            global err = kk - norm(vnn[:, iter + 1])           # Check voltage change (convergence criterion)
            global iter += 1;

            if iter > 10
                println("Complex linearization. After 10 iterations the error is:")
                break;
            end
        end

        # Reconstruct full voltage vector (slack + non-slack)
        local v_n = ones(ComplexF64, 3 * num_n, 1)
        v_n[n_slack] = vs
        v_n[n_other] = vnn[:, iter]

        
        local s_n = v_n .* conj(ybus * v_n)
        s_sub[tm, 1] = tm                                               # Time index
        s_sub[tm, 2] = sum(real(s_n[n_slack]) * p_base * 1000)          # Real power in kW
        s_sub[tm, 3] = sum(imag(s_n[n_slack]) * p_base * 1000)          # Reactive power in kVAr

        # Reinitialize voltage guess for next scenario
        global vnn = zeros(ComplexF64, 2715, 10)
        vnn[:, 1] = v_n[n_other]
    end
end

# Display results (first 5 scenarios)
display(s_sub[1:5, 2])   # Real power in kW
display(s_sub[1:5, 3])   # Reactive power in kVAr











