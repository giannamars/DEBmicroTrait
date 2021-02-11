using DEBmicroTrait
using CSV, DataFrames
using Gadfly, Cairo, Fontconfig
using Roots

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)

Genome_size = convert(Array{Float64,1}, df.Genome_size)
rrn_copies = df.rRNA_genes
gmax = log(2)./df.Min_gen_time
Gram_stain = convert(Array{String}, df.gram_stain)
#α = fill(0.0, size(Genome_size,1))
α = DEBmicroTrait.constrain_enzyme_production(ones(39), df.z_hydrolases)

V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)

yield = zeros(size(V_cs,1), 13)
rate = zeros(size(V_cs,1), 13)
ρs = zeros(size(V_cs,1))
mEs = zeros(size(V_cs,1))

for i in 1:size(V_cs,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρs[i] = ρ_p
    mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
end

plot(x=V_cs, y=mEs, Scale.x_log10, Scale.y_log10)

FCR = zeros(39)
for i in 1:size(FCR,1)
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
end

plot(x=V_cs, y=FCR, Scale.x_log10)
