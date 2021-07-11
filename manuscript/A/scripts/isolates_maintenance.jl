using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using HypothesisTests

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
########################################

########################################
Genome_size     = convert(Array{Float64,1}, df_isolates.Genome_size)
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
k_M    = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cell, Min_gen_time, Gram_stain)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_metabolism.jld", "kM", k_M)
df_out              = DataFrame()
df_out.k_M          = k_M
df_out.response     = df_isolates.Rhizosphere_response
########################################

########################################
# statistics
df_p       = filter(x->(x.response.=="positive"), df_out)
df_n       = filter(x->(x.response.=="negative"), df_out)
df_u       = filter(x->(x.response.=="undefined"), df_out)

k_M_kw     = KruskalWallisTest(df_p.k_M, df_n.k_M)
