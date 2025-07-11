using LinearAlgebra
using Printf
using MAT
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
B = imag(ybus);
G = real(ybus); 
ynn = feeder.ynn;                                            #Ynn matrix of the system 
znn= feeder.znn;                                             #Znn matrix of the system
yns = feeder.yns;
n_slack = feeder.n_slack;                                    #Slack nodes, phases A, B, and C
n_other = feeder.n_other;                                    #Other nodes
vnn=feeder.vn_initial                                        #Initial voltages (non-slack)
vs=feeder.vs_initial
p_base=feeder.p_base                                         #Base power
num_e=feeder.num_e                                           #Number of scenarios                                   
sref = -cargas[n_other,1]                                    #Load variation
pref = real(sref);
qref = imag(sref);


#Voltage initialization
v = ones(ComplexF64,3*num_n);
an = zeros(3*num_n);

an[n_other] = angle.(vnn);
an[n_slack] = angle.(vs);
v[n_slack] = abs.(vs);


err = 100;                           #Initial error
conv = zeros(10,1);                  #Error at each iteration
iter = 1; 


num_t = 3*num_n;
num_r = length(n_other);
H = zeros(num_t,num_t);                #Define Jacobian submatrices H, N, J, and L
N = zeros(num_t,num_t);
J = zeros(num_t,num_t);
L = zeros(num_t,num_t);

while err>1E-9
    local vn = v.*exp.(an*1im);
    global sn = vn.*conj(ybus*vn);
    local p = real(sn);
    local q = imag(sn);
    for k = 1:num_t                     #Build Jacobian
        for m = 1:num_t
            if k==m
                H[k,k] = -B[k,k]*v[k,1]*v[k,1]-q[k,1];            
                N[k,k] =  G[k,k]*v[k,1] + p[k,1]/v[k,1];
                J[k,k] = -G[k,k]*v[k,1]*v[k,1] + p[k,1];
                L[k,k] = -B[k,k]*v[k,1] + q[k,1]/v[k,1];    
            else
                akm = an[k,1]-an[m,1];
                N[k,m] =  v[k,1]*(G[k,m]*cos(akm)+B[k,m]*sin(akm));
                L[k,m] =  v[k,1]*(G[k,m]*sin(akm)-B[k,m]*cos(akm));
                H[k,m] =  L[k,m]*v[m,1];
                J[k,m] = -N[k,m]*v[m,1];
            end
        end
    end
    local dp = pref-p[n_other];
    local dq = qref-q[n_other];

    local Jac = [H[n_other,n_other] N[n_other,n_other] 
           J[n_other,n_other] L[n_other,n_other]];        #Assemble Jacobian
 
    local delta=[dp
           dq]
    local dx = Jac\delta;

    global an[n_other] = an[n_other]+dx[1:num_r];
    global v[n_other]  = v[n_other]+dx[num_r+1:end];
    global err = norm(dx);
    global conv[iter,1] = err;
    global iter = iter+1;
    if iter > 10
        println("Complex linearization. After 10 iterations the error is :")
        break;
    end
    
end

vn = v.*exp.(an*1im);
v_node = vn;
s_node = vn.*conj(ybus*vn);
p_loss = real(sum(s_node));
errores = conv;
iteraciones = iter;
scenario = 566;

end


println(iteraciones)
println(p_loss)

