using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using Gadfly
using JLD

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)
df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")

Genome_size = convert(Array{Float64,1}, df.Genome_size)
rrn_copies = convert(Array{Float64,1}, df.rRNA_genes)
Min_gen_time = df.Min_gen_time
gmax = log(2)./Min_gen_time
Gram_stain = convert(Array{String}, df.gram_stain)
#α = DEBmicroTrait.constrain_enzyme_production(ones(39), df.z_hydrolases)
α = zeros(39)
V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)

yield = zeros(size(V_cs,1), 100000)
rate = zeros(size(V_cs,1), 100000)
ρs = zeros(size(V_cs,1))
mEs = zeros(size(V_cs,1), 100000)

for i in 1:39
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
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

FCR = zeros(size(V_cs,1))
for i in 1:39
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
end

for i in 1:39
    rate = filter!(x->x>0.0,rate)
yield1 = filter!(x->x>0.0,yield[1,:])
plot(x=yield1, y=rate1, Scale.x_log10, Scale.y_log10, Guide.xlabel("yield"), Guide.ylabel("rate"))


JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/isolates_rate_yield_trade_off_2.jld",
         "yield", yield, "rate", rate, "mE", mEs, "Vc", V_cs, "fcr", FCR)
