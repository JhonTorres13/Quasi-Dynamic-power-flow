"""
select_scenario(feeder, tm)

This function selects a load scenario from a set of profiles.
It calculates the complex power (S = P + jQ) for each node and phase,
based on the load definitions and a given scenario index 'tm'.

# Arguments
-feeder: A feeder data structure
-tm: An integer indicating the scenario to select from the load profiles.

# Returns
-s_load: A column vector (3*num_n x 1) containing the complex power demand on each phase of each node.
"""

using LinearAlgebra
using Printf
using MAT
using LinearAlgebra
using Profile


function select_scenario(feeder, tm)
    num_n = feeder.num_n                                        # Total number of nodes
    s_load = zeros(ComplexF64, 3 * num_n, 1)
    for k = 1:feeder.num_d
        n1 = Int(feeder.loads[k, 1])
        ph = Int(feeder.loads[k, 2])
        pf = feeder.loads[k, 3]
        prof = Int(feeder.loads[k, 4])
        p = feeder.profiles[tm, prof]                        # Real power for this time step and profile
        q = p * sqrt((1 / pf)^2 - 1)                         # Calculate reactive power using the formula PF

        # Assign the complex power to the correct location in the vector                                                    
        if ph == 1                                           
            s_load[n1, 1] = p + q*im                         # Phase A at node n1
        elseif ph == 2
            s_load[n1 + num_n, 1] = p + q*im                 # Phase B at node n1
        elseif ph == 3
            s_load[n1 + 2 * num_n, 1] = p + q*im             # Phase C at node n1
        end
    end
    return s_load                                            # Returns the complete complex load vecto
end


 
