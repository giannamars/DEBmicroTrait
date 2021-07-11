using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly
using Roots
using HypothesisTests
using JLD

df_isolates = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")

Genome_size = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time = df_isolates.Min_gen_time
gmax = log(2)./Min_gen_time
Gram_stain = convert(Array{String}, df_isolates.gram_stain)
α  = DEBmicroTrait.constrain_enzyme_allocation(V_cell, Min_gen_time, Gram_stain, df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)


y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C
y_EM            = ones(size(V_cell,1))
assimilation      = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_assimilation.jld")
N_SB            = assimilation["NSB"][1:83,:]
K_D             = assimilation["KD"][1:83,:]
# ρ_ps            = zeros(size(y_DE,1), size(V_cell,1))
#
# for j in 1:size(y_DE,1)
#     if df_metabolites.Ontology[j] == "Sugars"
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[1]
#         end
#     elseif df_metabolites.Ontology[j] == "Organic acids"
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[2]
#         end
#     elseif df_metabolites.Ontology[j] == "Amino acids"
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[3]
#         end
#     elseif df_metabolites.Ontology[j] == "Fatty acids"
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[4]
#         end
#     elseif df_metabolites.Ontology[j] == "Nucleotides"
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[5]
#         end
#     else
#         for i in 1:size(V_cell,1)
#             find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
#             ρ_p = Roots.find_zero(find_ρ, 1.0)
#             closure = genome_distr[:,i]./sum(genome_distr[:,i])
#             ρ_ps[j,i] = ρ_p[1].*closure[6]
#         end
#     end
# end
#
# ρ_ps[ρ_ps.==0.0] .= 1e-8
#
# N_SB = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
# Vmax = @. 180.0*60^2*N_SB.*N_C
# K_D = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, df_metabolites.diffusivity)

n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)


id_microbe = 1
id_monomer = 1
p = DEBmicroTrait.init_batch_model_old(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)

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

DEBmicroTrait.calc_mass_balance(p,sol)
#
D    = [sol[i][1] for i in 1:length(sol.t)]
E    = [sol[i][2] for i in 1:length(sol.t)]
V    = [sol[i][3] for i in 1:length(sol.t)]
X    = [sol[i][4] for i in 1:length(sol.t)]
CO2  = [sol[i][5] for i in 1:length(sol.t)]

plot(x=sol.t, y=D)


BGEs = zeros(83, 39)
rs   = zeros(83, 39)
J_Vs = zeros(83, 39)
J_Es = zeros(83, 39)

for l = 1:39
    for m = 1:83
        id_microbe = l
        id_monomer = m
        p = DEBmicroTrait.init_batch_model_old(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)
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
        BGEs[m,l] = BGE_med

        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
        r_med = median(r)
        rs[m,l] = r_med

        J_V  = [DEBmicroTrait.biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, [V[i]])[1] for i in 1:size(sol.t,1)]
        J_V_med = median(J_V)
        J_Vs[m,l] = J_V_med
        J_E  = [DEBmicroTrait.biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, [E[i]])[1] for i in 1:size(sol.t,1)]
        J_E_med = median(J_E)
        J_Es[m,l] = J_E_med
    end
end

BGEs[BGEs.>1.0].=NaN
BGEs[BGEs.<0.05].=NaN

df_out = DataFrame()
df_out.BGE = vec(BGEs).+0.3
df_out.rate = vec(rs)
df_out.JV = vec(J_Vs)
df_out.JE = vec(J_Es)
df_out.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))

response = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_out.response = vec(response)

species = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     species[:,i] .= df_isolates.Isolate[i]
 end
df_out.species = vec(species)

class = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     class[:,i] .= df_isolates.Class[i]
 end
df_out.class = vec(class)
df_out.monomer = repeat(df_metabolites.Name, size(V_cell,1))


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/batch_model_isolates_old.csv", df_out)

df_isolates.Class
BGEs

df_p = filter(x->(x.response.=="positive") , df_out)
df_n = filter(x->(x.response.=="negative") , df_out)

df_p[isnan.(df_p.BGE), :BGE] .= 0.0
df_n[isnan.(df_n.BGE), :BGE] .= 0

median(df_p.BGE)
median(df_n.rate, na.rm=true)

df_p_sugars = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_out)

kw_vmax_sugars = KruskalWallisTest(df_n_sugars.BGE, df_p_sugars.BGE)
kw_K_sugars = KruskalWallisTest(df_n_sugars.rate, df_p_sugars.rate)
