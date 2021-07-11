using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly
using Roots
using HypothesisTests
using GLM

df_isolates = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props1.csv", DataFrame, missingstring="N/A")
df_metabolites = df_metabolites[(df_metabolites.Name.=="glucose") .| (df_metabolites.Name.=="salicylic acid") .| (df_metabolites.Name.=="alanine") .| (df_metabolites.Name.=="sinapyl alcohol") .| (df_metabolites.Name.=="oxalic acid") , :]

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

ρ_ps            = zeros(size(y_DE,1), size(V_cell,1))

for j in 1:size(y_DE,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[2]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[3]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[4]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[5]
        end
    else
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[6]
        end
    end
end


ρ_ps[ρ_ps.==0.0] .= 1e-8

N_SB_0 = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax = @. 180.0*60^2*N_SB_0.*N_C
K_D_0 = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, df_metabolites.diffusivity)


########################################
n_polymers = 0
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 1
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
########################################

u0 = zeros(p_set.dim+1)
u0[1+n_polymers:n_polymers+n_monomers] .= (25.0)/n_monomers
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] .= 1e-3
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 1e-3

tspan = (0.0,500.0)
condition(u,t,integrator) = u[1] - 1e-6
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition,affect!)

BGEs = zeros(5, 39, 3)
rs   = zeros(5, 39, 3)
D_ads_end = zeros(5, 39, 3)
du = zeros(p_set.dim+1)

for l = 1:39
    for m = 1:5
        for k = 1:3
            id_microbe = l
            id_monomer = m
            id_mineral = k

            K_D_M_glucose  = [0.024, 0.141, 0.067].*100
            K_D_M_alanine  = [0.037, 0.193, 0.059].*100
            K_D_M_salicyclic = [0.055, 0.103, 0.032].*100
            K_D_M_sinapyl = [0.008,0.009,0.006].*100
            K_D_M_oxalic = [0.101, 0.164, 0.105].*100
            K_D_M = hcat(K_D_M_alanine, K_D_M_glucose, K_D_M_salicyclic, K_D_M_sinapyl, K_D_M_oxalic)'
            K_D = hcat(K_D_0, reshape(K_D_M[:,id_mineral],5,1))
            N_SB = hcat(N_SB_0, ones(5,1))

            p = DEBmicroTrait.init_batch_model_minerals(id_microbe, id_monomer, id_mineral, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)

            prob  = ODEProblem(DEBmicroTrait.batch_model_minerals!,u0,tspan,p)
            sol   = solve(prob, alg_hints=[:stiff], callback=cb)
            #

            D    = [sol[i][1] for i in 1:length(sol.t)]
            E    = [sol[i][2] for i in 1:length(sol.t)]
            V    = [sol[i][3] for i in 1:length(sol.t)]
            X    = [sol[i][4] for i in 1:length(sol.t)]
            D_ads = [sol[i][5] for i in 1:length(sol.t)]
            CO2  = [sol[i][6] for i in 1:length(sol.t)]

            BR   = [DEBmicroTrait.batch_model_minerals!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
            BP   = [DEBmicroTrait.batch_model_minerals!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model_minerals!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
            BGE  = @. BP/(BP + BR)
            BGE_med = median(BGE[BGE.>=0.0])
            BGEs[m,l,k] = BGE_med

            r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
            r_med = median(r)
            rs[m,l,k] = r_med

            D_ads_end[m,l,k] = D_ads[end]
        end
    end
end
