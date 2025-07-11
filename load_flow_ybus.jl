"""Flujo de carga trifasico utilizando el método de punto fijo directamente """

using LinearAlgebra
using Printf
using MAT
using LinearAlgebra
using Profile
include("load_feeder.jl")          
include("select_scenario.jl")      

tm = 566                           # Select scenario index

nombre_ar = "FEEDER900.xlsx"       # Specify feeder data file
feeder = load_feeder(nombre_ar);   # Load feeder data from Excel
cargas = select_scenario(feeder, tm) # Get scenario-specific complex loads


# Extract network data from the feeder struct
num_n = feeder.num_n;              # Number of nodes
ybus = feeder.ybus;                # System admittance matrix
ynn = feeder.ynn;                  # Reduced admittance matrix (non-slack nodes)
znn = feeder.znn;                  # Inverse of Ynn
yns = feeder.yns;                  # Coupling admittance between slack and other nodes
n_slack = feeder.n_slack;          # Slack node indices
n_other = feeder.n_other;          # Non-slack node indices
vn = feeder.vn_initial             # Initial voltage guess for non-slack nodes
vs = feeder.vs_initial             # Slack bus voltages
p_base = feeder.p_base             # Base power (for per-unit normalization)
sn = -cargas[n_other];             # Complex power injection at non-slack nodes


@time begin                        # Start timing the core algorithm

err = 100;                         # Initialize error
conv = zeros(10,1);                # Store error at each iteration
iter = 1;                          # Iteration counter

# Fixed-point power flow iteration
while err > 1E-8
    kk = norm(vn)                                       # Store previous voltage norm
    global vn = znn * (conj(sn ./ vn) - yns * vs);      # Fixed-point update of voltages
    global err = kk - norm(vn);                         # Compute convergence error
    conv[iter,1] = err;                                 # Store error
    global iter = iter + 1;                             # Increment iteration count
    if iter > 10
        println("Complex linearization. After 10 iterations the error is :");
        break;
    end
end

# Build the full voltage vector including slack and other nodes
v_n = ones(ComplexF64, 3*num_n, 1);
v_n[n_slack] = vs;               # Assign slack voltages
v_n[n_other] = vn;               # Assign updated voltages for non-slack nodes

# Compute power at each node
s_n = v_n.*conj(ybus*v_n);      # Nodal complex powers
p_loss = real(sum(s_n));        # Total real power loss in the system
errorres = conv;                # Convergence history
iteraciones = iter;             # Final iteration count
scenario = 566;                 # Scenario index

end


println(iteraciones)            
println(p_loss)                


