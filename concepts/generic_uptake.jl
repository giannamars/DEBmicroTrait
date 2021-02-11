using DEBmicroTrait
using Roots
using Gadfly, Cairo, Fontconfig
using JLD
using CSV, DataFrames, Statistics


Genome_size = convert(Array{Float64,1}, LinRange(2e6, 9e6, 100))
V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain = repeat(["+"], size(Genome_size,1))
α = fill(0.0, size(Genome_size,1))


ρs = zeros(size(V_cs,1))
N_SBs = zeros(size(V_cs,1))
Ks = zeros(size(V_cs,1))

for i in 1:size(V_cs,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρs[i] = ρ_p
    N_SBs[i] = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cs[i]*ones(1), ρs[i]*ones(1,1), [Min_gen_time[i]], [Gram_stain[i]])[1]
    Ks[i] = DEBmicroTrait.specific_reference_affinity(V_cs[i]*ones(1), ρs[i]*ones(1,1), 1e-10*ones(1))[1]
end

Vmax = 180.0*60^2*N_SBs

df_isolates = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)

z_sugars = median(df_isolates.z_sugars)
z_organics = median(df_isolates.z_organic_acids)
z_aminos = median(df_isolates.z_amino_acids)
z_fattys = median(df_isolates.z_fatty_acids)
z_nucleos = median(df_isolates.z_nucleotides)
z_auxins = median(df_isolates.z_auxins)
genome_distr  = reshape(repeat(vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins), 100),6,100)

ρs_closure = DEBmicroTrait.transporters_closure(ρs, genome_distr)

N_SB = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cs, ρs_closure, Min_gen_time, Gram_stain)
Vmax_closure = 180.0*60^2*N_SB
K_closure = DEBmicroTrait.specific_reference_affinity(V_cs, ρs_closure, 1e-10*ones(6))


JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/generic_uptake.jld", "Vmax", Vmax, "K", Ks, "Vc", V_cs, "Vmaxc", Vmax_closure, "Kc", K_closure)
