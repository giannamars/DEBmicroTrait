using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")


Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)


k_M    = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
gmax   = log(2)./Min_gen_time
V_p    = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r    = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
k_E    = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
y_EV   = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)


save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_metabolism.jld", "kM", k_M, "kE", k_E, "yEV", y_EV)
