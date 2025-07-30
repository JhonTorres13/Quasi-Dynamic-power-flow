#|--------------------------------------------|
#| Power flow in power distribution networks  |
#| Load data from excel file                  |
#| By:  Jhon Torres                           |
#|      jhon.torres1@utp.edu.co               |
#|--------------------------------------------|


#librerias y importacion de archivos
using LinearAlgebra
using Printf
using LinearAlgebra
using Profile
using XLSX
using DataFrames

# Define a data structure to store feeder information
struct Datos
    graph::Matrix{Int64}             # Connection matrix between nodes
    z_line::Array{ComplexF64, 3}     # Line impedances
    Y_abc::Array{ComplexF64, 3}      # Line admittances
    ybus::Matrix{ComplexF64}         # Bus admittance matrix
    yns::Matrix{ComplexF64}          # Submatrix Y_NS (nodes to slack)
    ynn::Matrix{ComplexF64}          # Submatrix Y_NN (nodes to nodes)
    znn::Matrix{ComplexF64}          # Inverse of Y_NN (Z_NN)
    loads::Matrix{Float64}           # Load definitions
    profiles::Matrix{Float64}        # Load profiles
    xy::Matrix{Float64}              # XY coordinates
    vs_initial::Vector{ComplexF64}   # Initial slack bus voltages
    vn_initial::Vector{ComplexF64}   # Initial voltages for other nodes
    p_base::Float64                  # Base power
    v_base::Float64                  # Base voltage
    n_slack::Vector{Int}             # Slack node indices (phases A, B, and C)
    n_other::Vector{Int}             # All other node indices
    num_n::Int                       # Number of nodes
    num_l::Int                       # Number of lines
    num_d::Int                       # Number of loads
    num_e::Int                       # Number of load scenarios
end



function load_feeder(j)
    ruta="C:/Users/jt588/codigos_programacion_JULIA/Tesis/"*j

     # Read the Excel file and load the desired sheets
    datos = XLSX.readxlsx(ruta)
    

    lines=convert(Matrix{Float64}, datos[1][2:end,1:end])    
    codes =convert(Matrix{Float64}, datos[2][2:end,1:end])  
    general = convert(Matrix{Float64}, datos[6][2:end,1:end]) 
    xy = convert(Matrix{Float64}, datos[5][2:end,1:end])  

    p_base = general[1,2] / 3        # Base power per phase
    v_base = general[1,1] / sqrt(3)  # Line-to-neutral voltage
    z_base = v_base^2 / p_base        # Base impedance


    num_n = Int(maximum(maximum(lines[:,1:2])))  # Number of nodes
    num_l = length(lines[:,1])                   # Number of lines
    ybus = zeros(ComplexF64, 3*num_n, 3*num_n)    # Initialize Ybus matrix
    z_line = zeros(ComplexF64, 3, 3, num_l)       # Initialize line impedances
    Y_abc = zeros(ComplexF64, 3, 3, num_l)        # Initialize line admittances

    # Loop through each line and compute Y_abc, Z_line, and update Ybus

    for k = 1:num_l
        n1 = Int(lines[k,1]);
        n2 = Int(lines[k,2]);
        len = lines[k,3]/1000;
        cde = Int(lines[k,4]);
        r1 = codes[cde,2]*len;
        x1 = codes[cde,3]*len;
        r0 = codes[cde,4]*len;
        x0 = codes[cde,5]*len;
        B0 = codes[cde,6]*len*10^(-6);
        B1 = codes[cde,7]*len*10^(-6);

        zs = (2*r1+r0)/3 + (2*x1+x0)/3*1im;
        zm = (r0-r1)/3 + (x0-x1)/3*1im;  

        Bs = (2*B1+B0)/3*1im;
        Bm = (B0-B1)/3*1im; 

        Y_abc[:,:,k] = [Bs Bm Bm
                        Bm Bs Bm 
                        Bm Bm Bs]/(1/z_base);
        
        z_line[:,:,k] = [zs zm zm
                        zm zs zm 
                        zm zm zs]/z_base;
        yL = inv(z_line[:,:,k])
        Y_abc_2=(Y_abc[:,:,k])/2;
        nt1 = [n1,n1+num_n,n1+2*num_n];
        nt2 = [n2,n2+num_n,n2+2*num_n]; 

        # Update Ybus matrix using line admittance
        ybus[nt1,nt1] = ybus[nt1,nt1] + yL + Y_abc_2 ;
        ybus[nt1,nt2] = ybus[nt1,nt2] - yL;
        ybus[nt2,nt1] = ybus[nt2,nt1] - yL;
        ybus[nt2,nt2] = ybus[nt2,nt2] + yL+ Y_abc_2;    
    end

    vs =([exp(0*1im) 
        exp((-2*pi/3)*1im) 
        exp((2*pi/3)*1im)])*general[1,3];


    n_slack = [1,num_n+1,2*num_n+1];
    n_other = setdiff((1:3*num_n), n_slack)
    yns = ybus[n_other,n_slack];
    ynn = ybus[n_other,n_other];

    graph = Int.(lines[:,1:2]);            # Connection graph (from-to nodes)
    znn = inv(ynn);

    loads = convert(Matrix{Float64}, datos[3][2:end,1:end])  

    profiles = convert(Matrix{Float64}, (datos[4][2:end,1:end])/p_base/1000)  # Normalized profiles

    xy = xy[:,2:3];
    vs_initial = vs;
    vn_initial = vec(kron(vs, ones(num_n - 1, 1)))  # Initial voltages for non-slack nodes
    num_d = length(loads[:,1])                      # Number of loads
    num_e = length(profiles[:,1])                   # Number of load profiles

   
    # Create and return an instance of the Datos struct with all system data
    datos_sal = Datos(graph, z_line, Y_abc , ybus, yns, ynn, znn, loads, profiles, xy,vs_initial, vn_initial, p_base, v_base, n_slack, n_other,num_n, num_l, num_d, num_e)

    return datos_sal

end






