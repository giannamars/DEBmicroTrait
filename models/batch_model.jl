using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly

n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)

df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")

Genome_size = convert(Array{Float64,1}, LinRange(2e6, 9e6, 100))
Genome_size = [Genome_size[5], Genome_size[15], Genome_size[25], Genome_size[35], Genome_size[45]]

u0 = zeros(p_set.dim)
u0[1+n_polymers:n_polymers+n_monomers] .= (1e-5)/n_monomers
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] .= 1e-5
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 1e-5
tspan = (0.0,1000.0)
condition(u,t,integrator) = u[1] - 1e-8
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)
#
BGEs = zeros(5, 25)
rs   = zeros(5, 25)
J_Vs = zeros(5, 25)
J_Es = zeros(5, 25)

id_microbe = 1
id_monomer = 1
p = DEBmicroTrait.init_batch_model(id_microbe,id_monomer,df_metabolites, Genome_size, p_set)
#
prob  = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
sol   = solve(prob, alg_hints=[:stiff], callback=cb)
#

D    = [sol[i][1] for i in 1:length(sol.t)]
E    = [sol[i][2] for i in 1:length(sol.t)]
V    = [sol[i][3] for i in 1:length(sol.t)]
X    = [sol[i][4] for i in 1:length(sol.t)]
CO2  = [sol[i][5] for i in 1:length(sol.t)]

plot(x=sol.t, y=D)

#
for l = 1:5
    for m = 1:25
        id_microbe = l
        id_monomer = m
        p = DEBmicroTrait.init_batch_model(id_microbe,id_monomer,df_metabolites, Genome_size, p_set)
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

BGEs[BGEs.>1.0] .= NaN
rs


df_out = DataFrame()
BGEs_flat = vcat(BGEs[1,:], BGEs[2,:], BGEs[3,:], BGEs[4,:], BGEs[5,:])
rs_flat = vcat(rs[1,:], rs[2,:], rs[3,:], rs[4,:], rs[5,:])
JVs_flat = vcat(J_Vs[1,:], J_Vs[2,:], J_Vs[3,:], J_Vs[4,:], J_Vs[5,:])
JEs_flat = vcat(J_Es[1,:], J_Es[2,:], J_Es[3,:], J_Es[4,:], J_Es[5,:])
bugs = vcat(fill("1",25), fill("2",25), fill("3",25), fill("4",25), fill("5",25))

df_out.BGE = BGEs_flat
df_out.rate = rs_flat
df_out.JV = JVs_flat
df_out.JE = JEs_flat
df_out.size = bugs

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/batch_model_low.csv", df_out)
