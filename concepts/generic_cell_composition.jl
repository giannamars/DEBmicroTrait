using DEBmicroTrait
using CSV, DataFrames
using Gadfly, Cairo, Fontconfig

Genome_size = convert(Array{Float64,1}, LinRange(1.25e6, 1.12e7, 100))
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain = repeat(["+"], size(Genome_size,1))

p = IsolateComposition(Genome_size, Min_gen_time, Gram_stain)
v_N = 1.47e-27
V_DNA = Genome_size.*v_N

JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/generic_cell_composition.jld", "V_c", p.Cell_volume, "V_DNA", V_DNA, "V_p", p.Protein_volume, "V_r", p.Ribosome_volume, "V_env", p.Envelope_volume)
