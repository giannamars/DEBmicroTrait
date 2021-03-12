using DEBmicroTrait, Roots
using Gadfly
using DataFrames, CSV, GLM

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)

Genome_size = convert(Array{Float64,1}, df.Genome_size)
V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = convert(Array{Float64,1}, df.rRNA_genes)
Min_gen_time = df.Min_gen_time
gmax = log(2)./Min_gen_time
Gram_stain = convert(Array{String}, df.gram_stain)
#α = DEBmicroTrait.constrain_enzyme_production(ones(39), df.z_hydrolases)
α = zeros(39)

ρs_zero      = zeros(size(V_cs,1))
for i in 1:39
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], α[i])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρs_zero[i] = ρ_p
end
plot(x=V_cs, y=ρs_zero, Scale.x_log10, Scale.y_log10)

Vmax = zeros(size(V_cs,1))
K_Ds = zeros(size(V_cs,1))
k2p  = 180.0*60^2
for i in 1:39
    N_SB = DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cs[i]], ρs_zero[i]*ones(1,1), [Min_gen_time[i]], [Gram_stain[i]])
    Vmax[i] = k2p*N_SB[1,1]
    K_D = DEBmicroTrait.specific_reference_affinity([V_cs[i]], ρs_zero[i]*ones(1,1), 1e-10*ones(1))
    K_Ds[i] = K_D[1,1]
end
plot(x=Vmax, y=K_Ds, Scale.x_log10, Scale.y_log10)

df = DataFrame()
df.Vc = V_cs
df.VDNA = Genome_size
df.rrna = rrn_copies
df.gmax = gmax
df.Vmax = Vmax
df.K = K_Ds

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/isolates_steady_state_trade_off.jld",df)
