#|---------------------------------------------|
#| Power flow in power distribution networks   |
#| Active and reactive power in the substation |
#| Quasi dynamic simulation                    |
#| By:  Jhon Torres                            |
#|      jhon.torres1@utp.edu.co                |
#|---------------------------------------------|


"""Quasi-dynamic analysis of three-phase load flow on a 906-node European distribution feeder using backwar\forward sweep algorithm"""

# Libraries and file imports
using LinearAlgebra
using Printf
using LinearAlgebra
using Profile
using Plots
include("load_feeder.jl")
include("select_scenario.jl")

nombre_ar = "FEEDER900.xlsx"         # Feeder data filename
feeder = load_feeder(nombre_ar);     # Load feeder data from Excel file


@time begin
    # Extract system data
    num_n = feeder.num_n              # Number of nodes
    num_l = feeder.num_l              # Number of lines       
    vs = feeder.vs_initial            # Initial voltages for slack bus
    profiles = feeder.profiles        # Load profiles (per unit)
    loads = feeder.loads              # Load data: [node, phase, power factor, profile index]
    graph = feeder.graph              # Connection graph (topology)
    z_line = feeder.z_line            # Line impedance matrices
    ybus = feeder.ybus                # System Ybus matrix
    p_base = feeder.p_base            # Base power of the system
    num_d = feeder.num_d              # Total number of loads
    num_e = feeder.num_e              # Total number of scenarios
    n_slack = feeder.n_slack          # Indices of slack bus nodes (3 phases)

    # Initialize voltage matrix at all nodes (3-phase system)
    global vn = ones(ComplexF64, 3, num_n);  # Voltage matrix [phase, node]
    vn[:, 1] = vs                            # Set slack voltages for node 1

    for k = 2:num_n
        vn[:, k] = vs                        # Initialize all other nodes with slack voltage
    end

    s_sub = ones(num_e, 3)                   # Matrix to store substation power per scenario

    for tm = 1:num_e                            
        global err = 100                     
        global iter = 1                     

        while err > 1E-8                     # Convergence condition
            kk = norm(vn)                    # Voltage norm for convergence check

            # Current injection at load nodes
            global i_node = zeros(ComplexF64, 3, num_n);  # Node current matrix
            for k = 1:num_d
                n1 = Int(loads[k,1]);                                   # Node index
                ph = Int(loads[k,2]);                                   # Phase index (1, 2, or 3)
                pf = loads[k,3];                                        # Power factor
                prof = Int(loads[k,4]);                                 # Load profile index
                p = profiles[tm, prof];                                 # Active power (per unit)
                q = p * sqrt((1/pf)^2 - 1);                             # Reactive power based on PF
                global i_node[ph, n1] = conj((p + im*q) / vn[ph, n1]);  # Load current
            end

            # Backward sweep: accumulate downstream currents
            for k = num_l:-1:1
                n1 = graph[k,1];            # Sending node
                n2 = graph[k,2];            # Receiving node
                global i_node[:, n1] = i_node[:, n1] + i_node[:, n2];  # Current sum
            end

            # Forward sweep: update node voltages
            for k = 1:num_l
                n1 = graph[k,1];            # Sending node
                n2 = graph[k,2];            # Receiving node
                global vn[:, n2] = vn[:, n1] - z_line[:, :, k] * i_node[:, n2];  # Voltage drop
            end

            # Error computation for convergence check
            global err = kk - norm(vn);     
            global iter = iter + 1

            if iter > 10
                println("Complex linearization. After 10 iterations the error is :");       
                break
            end
        end

        # Prepare result storage
        tensiones = vn                     # Store voltages
        local a, b = size(vn')             # Get dimensions
        local Xv = reshape(vn', a * b, 1);
        local v_n = conj(Xv);              # Conjugate of voltage vector

        local v_node = v_n;
        local s_node = v_n .* conj(ybus * v_n);  # Calculate complex power per node

        # Store results for this scenario
        s_sub[tm, 1] = tm
        s_sub[tm, 2] = sum(real(s_node[n_slack]) * p_base * 1000)  # Real power in kW
        s_sub[tm, 3] = sum(imag(s_node[n_slack]) * p_base * 1000)  # Reactive power in kVAr

        # Reset voltage matrix for next scenario
        global vn = ones(ComplexF64, 3, num_n);
        vn .= tensiones
    end
end

# Display substation real and reactive power for the first 5 scenarios
display(s_sub[1:5, 2])
display(s_sub[1:5, 3])









