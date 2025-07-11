using LinearAlgebra
using Printf
using MAT
using LinearAlgebra
using Profile


function select_scenario(feeder, tm)
    num_n = feeder.num_n
    s_load = zeros(ComplexF64, 3 * num_n, 1)
    for k = 1:feeder.num_d
        n1 = Int(feeder.loads[k, 1])
        ph = Int(feeder.loads[k, 2])
        pf = feeder.loads[k, 3]
        prof = Int(feeder.loads[k, 4])
        p = feeder.profiles[tm, prof]
        q = p * sqrt((1 / pf)^2 - 1)
        
        if ph == 1
            s_load[n1, 1] = p + q*im
        elseif ph == 2
            s_load[n1 + num_n, 1] = p + q*im
        elseif ph == 3
            s_load[n1 + 2 * num_n, 1] = p + q*im
        end
    end
    return s_load
end
 