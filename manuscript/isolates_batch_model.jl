using DEBmicroTrait
using JLD
using DifferentialEquations
using Gadfly
using CSV, DataFrames

n_monomers = 84
n_isolates = 39
n_enzymes  = 1

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)

df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")
uptake          = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_assimilation.jld")
enzymes         = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_enzymes.jld")

n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)


id_monomer  = 1
id_isolate  = 1

α_X    = [enzymes["alpha"][id_isolate]]
met    = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_metabolism.jld")
k_M    = [met["kM"][id_isolate]]
k_E    = [met["kE"][id_isolate]]
y_EV   = [met["yEV"][id_isolate]]
y_EX   = y_EV
y_EM   = ones(1)
f_αX   = ones(n_enzymes)

p_met  = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX)

N_SB   = uptake["NSB"][id_monomer, id_isolate]*ones(1,1)
K_D    = uptake["KD"][id_monomer, id_isolate]*ones(1,1)
N_C    = [df_metabolites.N_C[id_monomer]]
y_DE   = [df_metabolites.y_DE[id_monomer]]

p_ass             = AssimilationC(N_SB,K_D,y_DE,N_C)

d            = 0.23.*exp.(0.88.*gmax)
γ_V_0        = d[1].*ones(1)
γ_V_1        = 1.0.*ones(1)
γ_X          = 2.5e-4.*ones(n_enzymes)
γ_D_ads      = zeros(1)
γ_X_ads      = zeros(1)
f_ED         = ones(1)
f_VD         = ones(1)
f_VP         = ones(1)
f_V          = ones(1)
f_XD         = ones(1)
f_XP         = zeros(1)
f_X          = zeros(1)



p_turn       = Turnover(γ_V_0,γ_V_1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)

p            = Params(p_set, p_met, p_ass, nothing, p_turn)

V_0          = 1e9*0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)

u0                                                                         = zeros(p_set.dim)
u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.0
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*V_0[id_isolate]
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*V_0[id_isolate]

tspan                     = (0.0,250.0)
condition(u,t,integrator) = u[1] - 1e-6
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)

prob  = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
sol   = solve(prob, alg_hints=[:stiff], callback=cb)


D    = [sol[i][1] for i in 1:length(sol.t)]
E    = [sol[i][2] for i in 1:length(sol.t)]
V    = [sol[i][3] for i in 1:length(sol.t)]
X    = [sol[i][4] for i in 1:length(sol.t)]
CO2  = [sol[i][5] for i in 1:length(sol.t)]
Bio = E.+V
plot(x=sol.t, y=D, Geom.line)

du   = zeros(p_set.dim)
BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
BGE  = @. BP/(BP + BR)

plot(x=gmax, y=d, Geom.point)

BGE



dD     = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[1] for i in 1:size(sol.t,1)]
dE     = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] for i in 1:size(sol.t,1)]
dV     = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
dX     = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[4] for i in 1:size(sol.t,1)]
dCO2   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]

mass_balance = dD .+ dE .+ dV .+ dX .+ dCO2
