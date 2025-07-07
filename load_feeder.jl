#librerias y importacion de archivos
using LinearAlgebra
using Printf
using MAT
using LinearAlgebra
using Profile
using XLSX
using DataFrames


# Definir el tipo de datos para almacenar los datos del alimentador
struct Datos
    graph::Matrix{Int64}  # Matriz de conexiones entre nodos
    z_line::Array{ComplexF64, 3}  # Impedancias de línea
    Y_abc::Array{ComplexF64, 3}   # Admitancias de línea
    ybus::Matrix{ComplexF64}  # Matriz de admitancia del bus
    yns::Matrix{ComplexF64}  # Admitancias de nodos 
    ynn::Matrix{ComplexF64}  # Admitancias de nodos 
    znn::Matrix{ComplexF64}  # Impedancias
    loads::Matrix{Float64}  # Cargas
    profiles::Matrix{Float64}  # Perfiles de carga
    xy::Matrix{Float64}  # Coordenadas XY
    vs_initial::Vector{ComplexF64}  # Voltajes iniciales slack
    vn_initial::Matrix{ComplexF64}  # Voltajes iniciales otros nodos
    p_base::Float64  # Potencia base
    v_base::Float64  # Voltaje base
    n_slack::Vector{Int}   # Número de nodo slack A, B y C
    n_other::Vector{Int}  # Números de otros nodos
    num_n::Int  # Número de nodos
    num_l::Int  # Número de líneas
    num_d::Int  # Número de cargas
    num_e::Int  # Número de escenarios
end


function load_feeder(j)
    ruta="C:/Users/jt588/codigos_programacion_JULIA/Tesis/"*j
    # Lee el archivo Excel y carga la hoja deseada (por ejemplo, la primera hoja)
    datos = XLSX.readxlsx(ruta)
    

    lines=convert(Matrix{Float64}, datos[1][2:end,1:end])    
    codes =convert(Matrix{Float64}, datos[2][2:end,1:end])  
    general = convert(Matrix{Float64}, datos[6][2:end,1:end]) 
    xy = convert(Matrix{Float64}, datos[5][2:end,1:end])  

    p_base = general[1,2]/3;         # potencia nominal por fase
    v_base = general[1,1]/sqrt(3);   # voltaje linea a neutro
    z_base = v_base^2/p_base;


    num_n = Int(maximum(maximum(lines[:,1:2])));
    num_l = length(lines[:,1]);
    ybus = zeros(ComplexF64,3*num_n,3*num_n);
    z_line = zeros(ComplexF64,3,3,num_l);
    Y_abc = zeros(ComplexF64,3,3,num_l);

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

    graph = Int.(lines[:,1:2]);
    znn = inv(ynn);

    loads = convert(Matrix{Float64}, datos[3][2:end,1:end])  

    profiles = convert(Matrix{Float64}, (datos[4][2:end,1:end])/p_base/1000)

    xy = xy[:,2:3];
    vs_initial = vs;
    vn_initial = kron(vs,ones(num_n-1,1));
    num_d = length(loads[:,1]);
    num_e = length(profiles[:,1]);


   
    # Crear una instancia de Feeder para almacenar los datos
    datos_sal = Datos(graph, z_line, Y_abc , ybus, yns, ynn, znn, loads, profiles, xy,vs_initial, vn_initial, p_base, v_base, n_slack, n_other,num_n, num_l, num_d, num_e)

    return datos_sal

end







