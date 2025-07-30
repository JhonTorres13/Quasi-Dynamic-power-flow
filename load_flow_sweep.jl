#|------------------------------------------------|
#| Power flow in power distribution networks      |
#| Load flow backwar\forward sweep algorithm      |
#| By:  Jhon Torres                               |
#|      jhon.torres1@utp.edu.co                   |
#|------------------------------------------------|

# Libraries and file import
using LinearAlgebra
using Printf
using MAT
using LinearAlgebra
using Profile
include("load_feeder.jl")
include("select_scenario.jl")

tm = 566                             # Select scenario
nombre_ar = "FEEDER900.xlsx"         # Load feeder data
feeder = load_feeder(nombre_ar);     # Load data from Excel
cargas = select_scenario(feeder, tm) # Select scenario load profile


@time begin
num_n = feeder.num_n;                # Number of nodes
num_l = feeder.num_l;                # Number of lines       
vs = feeder.vs_initial;              # Initial voltages (slack)
profiles = feeder.profiles           # Load profiles (per unit)
loads = feeder.loads                 # Load data (nodes, phases, PF, profile index)
graph = feeder.graph                 # Graph of node connections
z_line = feeder.z_line               # Line impedance matrices
ybus = feeder.ybus;                  # Ybus matrix of the system
p_base = feeder.p_base               # Base power
num_d = feeder.num_d                 # Number of loads
num_e = feeder.num_e                 # Number of scenarios
n_slack = feeder.n_slack             # Slack node indices

# Voltage initialization at all nodes
global vn = ones(ComplexF64, 3, num_n); # Initialize 3-phase voltage matrix
vn[:, 1] = vs                          # Assign slack voltages to node 1

for k = 2:num_n
    vn[:, k] = vs                     # Initialize all other nodes with slack voltage
end

err = 100;                            # Initial error
conv = zeros(10, 1);                  # Error per iteration
iter = 1;                             # Iteration counter

# Fixed-point backward-forward sweep method
while err > 1E-9
    kk = norm(vn)                     # Store previous norm of voltage matrix

    # Current injection at load nodes
    global i_node = zeros(ComplexF64, 3, num_n); # Initialize current matrix
    for k = 1:num_d
        n1 = Int(loads[k, 1]);        # Node index
        ph = Int(loads[k, 2]);        # Phase (1, 2, or 3)
        pf = loads[k, 3];             # Power factor
        prof = Int(loads[k, 4]);      # Profile index
        p = profiles[tm, prof];       # Active power (per unit)
        q = p * sqrt((1 / pf)^2 - 1); # Reactive power from PF
        global i_node[ph, n1] = conj((p + im * q) / vn[ph, n1]); # Current injection
    end

    # Backward sweep: sum downstream currents
    for k = num_l:-1:1
        n1 = graph[k, 1];
        n2 = graph[k, 2];
        global i_node[:, n1] = i_node[:, n1] + i_node[:, n2];
    end

    # Forward sweep: update node voltages
    for k = 1:num_l
        n1 = graph[k, 1];
        n2 = graph[k, 2];
        global vn[:, n2] = vn[:, n1] - z_line[:, :, k] * i_node[:, n2];        
    end

    # Compute convergence error
    global err = kk - norm(vn);           # Error as difference in voltage norm
    global conv[iter, 1] = err;           # Store error value
    global iter = iter + 1;               # Increase iteration count

    if iter > 10
        println("Complex linearization. After 10 iterations the error is: ", err);       
        break
    end
end

# Flatten vn matrix to vector for compatibility with Ybus model
a, b = size(vn');
Xv = reshape(vn', a * b, 1);              # Convert to column vector
v_n = conj(Xv);                           # Complex conjugate to match power flow convention

v_node = v_n;                            # Node voltages
s_node = v_n .* conj(ybus * v_n);        # Node complex powers
p_loss = real(sum(s_node));              # Total real power losses
errores = conv;                          # Convergence history
iteraciones = iter;                      # Number of iterations
                 
end

println((p_loss))  # Print total real power losses










