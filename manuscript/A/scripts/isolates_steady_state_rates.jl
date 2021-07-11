using DEBmicroTrait
using JLD
using Gadfly

assimilation      = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_assimilation.jld")
enzymes           = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_enzymes.jld")
maintenance       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_metabolism.jld")
protein_synthesis = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_protein_synthesis.jld")


y_DE              = assimilation["yDE"]
N_C               = assimilation["NC"]
k2p               = 180.0*60^2
N_SB              = assimilation["NSB"]

k_M               = maintenance["kM"]
y_EM              = assimilation["yEM"]

k_E               = protein_synthesis["kE"]
y_EV              = protein_synthesis["yEV"]

y_EX              = y_EV
α_X               = enzymes["alpha"]

id_monomer        = 66
r_ss_all = zeros(39)
x_ss_all = zeros(39)

for i = 1:39
    id_isolate = i
    j_ED_max         = (1.0 .- y_DE[id_monomer])*N_C[id_monomer]*k2p*N_SB[id_monomer, id_isolate]
    r_ss             = (j_ED_max - k_M[id_isolate]*y_EM[id_isolate])/(j_ED_max/k_E[id_isolate] + (1+α_X[id_isolate])*y_EV[id_isolate])
    x_ss             = α_X[id_isolate]*y_EV[id_isolate]*r_ss/y_EX[id_isolate]
    r_ss_all[i]      = r_ss
    x_ss_all[i]      = x_ss
end



plot(x=x_ss_all[x_ss_all.>0.0], y=r_ss_all[r_ss_all.>0.0], Geom.point, Scale.y_log10, Scale.x_log10)


id_monomer = 1
id_microbe = 2
j_ED_max         = (1.0 .- y_DE[id_monomer])*N_C[id_monomer]*k2p*N_SB[id_monomer, id_isolate]
