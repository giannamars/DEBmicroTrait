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
Min_gen_time    = df_isolates.Min_gen_time
gmax            = log(2)./Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
########################################

########################################
γ_V_0 = DEBmicroTrait. max_specific_death_rate(Min_gen_time::Vector{Float64})
########################################

########################################
dry_mass        = 0.47*DEBmicroTrait.cell_volume_to_dry_mass(V_cell, gmax, Gram_stain)
ρ_bulk          = 1.21 # g/cm^3
Bio_0           = 1e9*1e6*ρ_bulk*dry_mass./12.011
γ_V_1           = median(Bio_0)
########################################

########################################
# I/O
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_turnover.jld", "gV0", γ_V_0, "gV1", γ_V_1)
df_out              = DataFrame()
df_out.gV0          = γ_V_0
df_out.response     = df_isolates.Rhizosphere_response
########################################

########################################
# statistics
df_p       = filter(x->(x.response.=="positive"), df_out)
df_n       = filter(x->(x.response.=="negative"), df_out)
df_u       = filter(x->(x.response.=="undefined"), df_out)

k_M_kw     = KruskalWallisTest(df_p.gV0, df_n.gV0)
