using DEBmicroTrait
using CSV, DataFrames
using Gadfly, Cairo, Fontconfig

Genome_size = convert(Array{Float64,1}, LinRange(1.25e6, 1.12e7, 100))
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Protein_volume = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Ribosome_volume = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)

k_E_prediction = DEBmicroTrait.translation_power(Protein_volume, Ribosome_volume, Min_gen_time)

plot(x=V_cell, y=k_E_prediction, Geom.point, Scale.x_log10, Scale.y_log10)
