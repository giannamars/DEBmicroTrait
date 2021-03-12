using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)
df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")

n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)


Genome_size = convert(Array{Float64,1}, df.Genome_size)
rrn_copies = convert(Array{Float64,1}, df.rRNA_genes)
Min_gen_time = df.Min_gen_time
gmax = log(2)./Min_gen_time
Gram_stain = convert(Array{String}, df.gram_stain)
α = DEBmicroTrait.constrain_enzyme_production(ones(39), df.z_hydrolases)

id_microbe = 1
id_monomer = 1
p = DEBmicroTrait.init_batch_model(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set)

u0 = zeros(p_set.dim)
u0[1+n_polymers:n_polymers+n_monomers] .= (25.0)/n_monomers
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] .= 1e-3
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 1e-3
du = zeros(p_set.dim)
tspan = (0.0,500.0)
condition(u,t,integrator) = u[1] - 1e-6
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

prob  = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
sol   = solve(prob, alg_hints=[:stiff], callback=cb)
#
D    = [sol[i][1] for i in 1:length(sol.t)]
E    = [sol[i][2] for i in 1:length(sol.t)]
V    = [sol[i][3] for i in 1:length(sol.t)]
X    = [sol[i][4] for i in 1:length(sol.t)]
CO2  = [sol[i][5] for i in 1:length(sol.t)]

plot(x=sol.t, y=V)


BGEs = zeros(39, 25)
rs   = zeros(39, 25)
J_Vs = zeros(39, 25)
J_Es = zeros(39, 25)

for l = 1:39
    for m = 1:25
        id_microbe = l
        id_monomer = m
        p = DEBmicroTrait.init_batch_model(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set)
        #
        prob  = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol   = solve(prob, alg_hints=[:stiff], callback=cb)
        #

        D    = [sol[i][1] for i in 1:length(sol.t)]
        E    = [sol[i][2] for i in 1:length(sol.t)]
        V    = [sol[i][3] for i in 1:length(sol.t)]
        X    = [sol[i][4] for i in 1:length(sol.t)]
        CO2  = [sol[i][5] for i in 1:length(sol.t)]

        BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE  = @. BP/(BP + BR)
        BGE_med = median(BGE[BGE.>=0.0])
        BGEs[l,m] = BGE_med

        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
        r_med = median(r)
        rs[l,m] = r_med

        J_V  = [DEBmicroTrait.biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, [V[i]])[1] for i in 1:size(sol.t,1)]
        J_V_med = median(J_V)
        J_Vs[l,m] = J_V_med
        J_E  = [DEBmicroTrait.biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, [E[i]])[1] for i in 1:size(sol.t,1)]
        J_E_med = median(J_E)
        J_Es[l,m] = J_E_med
    end
end

df_out = DataFrame()
df_out.BGE = vec(BGEs)
df_out.rate = vec(rs)
df_out.JV = vec(J_Vs)
df_out.JE = vec(J_Es)
df_out.response = repeat(df.Rhizosphere_response, 25)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/batch_model_isolates.csv", df_out)


tmp = convert(Array{Int64,1}, 1:39)
positives = tmp[df.Rhizosphere_response .== "positive"]
negatives = tmp[df.Rhizosphere_response .== "negative"]

using HypothesisTests
KruskalWallisTest(vec(BGEs[positives, :]), vec(BGEs[negatives, :]))
