using DEBmicroTrait
using CSV, DataFrames
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props1.csv", DataFrame, missingstring="N/A")
d_assimilation  = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_assimilation.jld")
d_enzymes       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_enzymes.jld")
########################################

########################################
# traits
ρ_ps            = d_assimilation["rho"]
α               = d_enzymes["alpha"]
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
y_EM            = ones(size(V_cell,1))
y_EX            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C
########################################

########################################
m_E_max         = DEBmicroTrait.steady_state_reserve_density(ρ_p[1,:], V_cell, Min_gen_time, Gram_stain, y_DE[1], N_C[1])
gmax            = log(2)./Min_gen_time
V_p             = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r             = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E             = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
j_EA_m          = m_E_max.*k_E
k_M             = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
y_EV            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
r_m             = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+(1+α)*y_EV)
x_m             = @. α*y_EV/y_EX*r_m
yield           = @. r_m*(1+m_E_max)/j_EA_m

plot(x=yield, y=r_m, Geom.point, Scale.x_log10, Scale.y_log10)
