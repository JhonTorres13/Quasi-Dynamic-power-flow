#|-------------------------------------------------------------|
#| Power flow in power distribution networks                   |
#| Load flow using complex linear with method fixed point      |
#| By:  Jhon Torres                                            |
#|      jhon.torres1@utp.edu.co                                |
#|-------------------------------------------------------------|


#libraries and file import
using LinearAlgebra
using Printf
using LinearAlgebra
using Profile
include("load_feeder.jl")
include("select_scenario.jl")


tm=566                            #Select scenario

nombre_ar="FEEDER900.xlsx"        #Load feeder data
feeder=load_feeder(nombre_ar);    #Load data from Excel
cargas=select_scenario(feeder,tm)


@time begin                                                  #Timer
num_n = feeder.num_n;                                        #Number of nodes
ybus = feeder.ybus;                                          #Ybus matrix of the system
ynn = feeder.ynn;                                            #Ynn matrix of the system 
znn= feeder.znn;                                             #Znn matrix of the system
yns = feeder.yns;
n_slack = feeder.n_slack;                                    #Slack nodes phases A, B and C
n_other = feeder.n_other;                                    #Other nodes
vn=feeder.vn_initial                                         #Initial voltages (non-slack)
vs=feeder.vs_initial                                         #Initial voltages (slack)
p_base=feeder.p_base                                         #Base power
num_e=feeder.num_e                                           #Number of scenarios
s_load = cargas                                              #Load variation

#Voltage initialization
v = ones(ComplexF64,3*num_n,1);
v[n_slack] = vs;
v[n_other] = vn;


err = 100;                       #Error
conv = zeros(10,1);              #Error per iteration
iter = 1; 
#Inv_A=zeros(ComplexF64,3*(num_n-1),3*(num_n-1));

while err > 1E-9
    local s = v.*conj(ybus*v);                                  #Nodal power
    local B = yns*vs+ynn*vn;                                    #Compute matrix JB
    local C = -conj(s[n_other]+s_load[n_other]);             #Compute complex powers
    local b= 1. ./ adjoint(vn);                                 #Inverse of diag(vn') as a vector
    local Inv_A=znn.*b                                          #Compute Inv A
    

    ## Fixed-point method proposal
    global dv = zeros(ComplexF64,3*(num_n-1),1)                 #Initialize the solution
    for iter = 1:15

        if iter==1
            dv_new=C;
         else 
            dv_new = Inv_A*(C-B.*conj(dv));  
         end

        # Check convergence
        if norm(dv_new - dv) < 1e-6
            break
        end
        global dv = dv_new

    end
    global vn = vn+dv;                    #Update voltages
    v[n_other]=vn;                         
    global err = norm(dv);                #Compute error
    conv[iter,1] = err;                   #Store error
    global iter = iter + 1;               #Update iteration count
    if iter > 10
        println("Complex linearization. After 10 iterations the error is : ", err);
        break;
    end
end

v_node = v;                           #Voltages
s_node = v.*conj(ybus*v);             #Nodal power
p_loss = real(sum(s_node));           #Losses
errores = conv;                       #Errors
iteraciones=iter                      #Iteration count

end

println("The system losses are: ", p_loss)
println(iteraciones)




