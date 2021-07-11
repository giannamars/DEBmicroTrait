using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly
using Roots
using HypothesisTests

df_isolates = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")
df_metabolites = df_metabolites[1:6,:]
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

N_SB = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax = @. 180.0*60^2*N_SB.*N_C
K_D = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, df_metabolites.diffusivity)

n_polymers = 0
n_monomers = 6
n_microbes = 39
n_enzymes  = 1
n_minerals = 0
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)

p = DEBmicroTrait.init_hydroponics(df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)

u0 = zeros(p_set.dim)
u0[1+n_polymers:n_polymers+n_monomers] .= (25.0)/n_monomers
u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] .= 1e-1
u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 1e-1


xt, TOC = DEBmicroTrait.exudation_io("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/TOC_data_Hopland_Avena.csv")
class_norm = DEBmicroTrait.exudation_class_io("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/avena_exudation.csv")

TOC_amino = TOC.*class_norm[1,:]
max_exudation_rate_amino(t) = max_exudation_rate(xt, TOC_amino, t)
TOC_organic = TOC.*class_norm[2,:]
max_exudation_rate_organic(t) = max_exudation_rate(xt, TOC_organic, t)
TOC_nucleotide = TOC.*class_norm[3,:]
max_exudation_rate_nucleotide(t) = max_exudation_rate(xt, TOC_nucleotide, t)
TOC_sugar = TOC.*class_norm[4,:]
max_exudation_rate_sugar(t) = max_exudation_rate(xt, TOC_sugar, t)
TOC_auxin = TOC.*class_norm[5,:]
max_exudation_rate_auxin(t) = max_exudation_rate(xt, TOC_auxin, t)
TOC_fatty = TOC.*class_norm[6,:]
max_exudation_rate_fatty(t) = max_exudation_rate(xt, TOC_fatty, t)

ts = collect(range(0,12*7*24,step=1))
tmp_ts = zeros(size(ts,1))
for i in 1:size(ts,1)
    if mod(ts[i], 24) == 0
        tmp_ts[i] = ts[i]
    end
end
diurnal_ts = filter(!iszero, tmp_ts)

cb_amino = FunctionCallingCallback((u,t,integrator)->integrator.u .+= max_exudation_rate_amino(t), func_everystep=true)


tspan = (0.0,100.0)

prob  = ODEProblem(DEBmicroTrait.hydroponics!,u0,tspan,p)
sol   = solve(prob, alg_hints=[:stiff], cb=cb_amino, tstops=diurnal_ts)

D1 = [sol[i][9] for i in 1:length(sol.t)]
plot(x=sol.t, y=D1, Geom.point)
