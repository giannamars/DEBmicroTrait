using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly

df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")
df_pold = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/IsolateMetadata190309.csv", DataFrame, missingstring="N/A")
df = df_pold[completecases(df_pold), :]

Genome_size = convert(Array{Float64,1}, df.actual_genome_size_bp)
rrn_copies  = convert(Array{Float64,1}, df.rrN)
Min_gen_time = convert(Array{Float64,1}, df.MinGT)
Gram_stain = convert(Array{String}, df.Gram_stain)
α = zeros(size(Genome_size,1))

n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)


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

D    = [sol[i][1] for i in 1:length(sol.t)]
E    = [sol[i][2] for i in 1:length(sol.t)]
V    = [sol[i][3] for i in 1:length(sol.t)]
X    = [sol[i][4] for i in 1:length(sol.t)]
CO2  = [sol[i][5] for i in 1:length(sol.t)]

plot(x=sol.t, y=D)

BGEs = zeros(size(Genome_size,1), 25)
rs   = zeros(size(Genome_size,1), 25)
J_Vs = zeros(size(Genome_size,1), 25)
J_Es = zeros(size(Genome_size,1), 25)

for l = 1:size(Genome_size,1)
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
df_out.BGE = vec(BGEs[:,[1,11,12]])
df_out.rate = vec(rs[:,[1,11,12]])
df_out.medium = vcat(repeat(["Glucose"], 42), repeat(["Pyruvate"], 42), repeat(["Succinate"], 42))


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/batch_model_pold.csv", df_out)
