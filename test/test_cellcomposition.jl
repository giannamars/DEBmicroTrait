using DEBmicroTrait, Test

Genome_size = convert(Array{Float64,1}, df.Genome_size)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain = repeat(["+"], 100)

p = IsolateComposition(Genome_size, Min_gen_time, Gram_stain)
