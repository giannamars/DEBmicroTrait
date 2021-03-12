using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using Gadfly, Cairo, Fontconfig
using JLD

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/ulas_terrestrial_genomes.csv", DataFrame)

Genome_size = convert(Array{Float64,1}, df.genomesize)
Genome_size = Genome_size[1:end .!= 61]
Genome_size = Genome_size[1:end .!= 576]
Genome_size = Genome_size[1:end .!= 1117]
Min_gen_time = convert(Array{Float64,1}, df.mingentime)
Min_gen_time = Min_gen_time[1:end .!= 61]
Min_gen_time = Min_gen_time[1:end .!= 576]
Min_gen_time = Min_gen_time[1:end .!= 1117]
gmax = log(2)./Min_gen_time
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
Gram_stain = repeat(["+"], size(Genome_size,1))
α = fill(0.0, size(Genome_size,1))
V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)


yield = zeros(size(V_cs,1), 100000)
rate = zeros(size(V_cs,1), 100000)
ρs = zeros(size(V_cs,1))
mEs = zeros(size(V_cs,1), 100000)

for i in 1:size(Genome_size,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    #print(i)
    ρs[i] = ρ_p
    mE_max = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    rtmp = convert(Array{Float64,1}, LinRange(1e-5, 1e1, 100000))
    filter!(x->x<=gmax[i],rtmp)
    mEs_tmp = zeros(size(rtmp,1))
    yield_tmp = zeros(size(rtmp,1))
    rate_tmp = zeros(size(rtmp,1))
    mEs_tmp = DEBmicroTrait.reserve_density_curve(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    yield_tmp, rate_tmp = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    for h in 1:size(rtmp,1)
        yield[i,h] = yield_tmp[h]
        rate[i,h] = rate_tmp[h]
        mEs[i,h] = mEs_tmp[h]
    end
end

FCR = zeros(size(Genome_size,1))
for i in 1:size(Genome_size,1)
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
end

FCR_filter = FCR[FCR .<= 1.0]
V_cs_filter = V_cs[FCR .<= 1.0]

p1 = plot(x=V_cs_filter, y=FCR_filter, Scale.x_log10, Scale.y_log10,
     Guide.xlabel("Cell volume [m^3]"), Guide.ylabel("FCR"))

p1 |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/FCR_terrestrial_genomes.png", dpi=300)
