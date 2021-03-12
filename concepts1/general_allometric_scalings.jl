using DEBmicroTrait, Roots
using Gadfly
using DataFrames, GLM

Genome_size  = convert(Array{Float64,1}, LinRange(2e6, 1e7, 100))
V_cs         = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies   = DEBmicroTrait.genome_size_to_rRNA_copy_number((Genome_size))
gmax         = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain   = repeat(["+"], size(Genome_size,1))
