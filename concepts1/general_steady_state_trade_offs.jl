using DEBmicroTrait, Roots
using Gadfly
using DataFrames, CSV, GLM

Genome_size  = convert(Array{Float64,1}, LinRange(2e6, 1e7, 100))
V_cs         = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies   = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax         = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain   = repeat(["+"], size(Genome_size,1))

#########################################################################################################
# Transporters and enzyme investment are not in trade-off conflict in growth rate optimization
α            = convert(Array{Float64,1}, LinRange(0.0, 1.0, 10))
ρs           = zeros(size(V_cs,1), size(α,1))
for j in 1:10
    for i in 1:100
        find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], α[j])
        ρ_p = Roots.find_zero(find_ρ, 1.0)
        ρs[i,j] = ρ_p
    end
end
plot(x=α, y=ρs[2,:], Scale.x_log10, Scale.y_log10)
#########################################################################################################
# Optimize transporters for α = 0.0
α            = fill(0.0, size(V_cs, 1))
ρs_zero      = zeros(size(V_cs,1))
for i in 1:100
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], α[i])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρs_zero[i] = ρ_p
end
plot(x=V_cs, y=ρs_zero, Scale.x_log10, Scale.y_log10)
#########################################################################################################
#
Vmax = zeros(size(V_cs,1))
K_Ds = zeros(size(V_cs,1))
k2p  = 180.0*60^2
for i in 1:100
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

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/general_steady_state_trade_off.jld",df)


ols = lm(@formula(log(Vmax) ~ log(Vc)), df)
layer1 = layer(df, x=:Vc, y=:Vmax)
y_Vmax = exp(coef(ols)[1]).*df.Vc.^(coef(ols)[2])
layer2 = layer(x=df.Vc, y=y_Vmax)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

ols = lm(@formula(log(K) ~ log(Vc)), df)
layer1 = layer(df, x=:Vc, y=:K)
y_K = exp(coef(ols)[1]).*df.Vc.^(coef(ols)[2])
layer2 = layer(x=df.Vc, y=y_K)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

ols = lm(@formula(log(Vmax) ~ log(K)), df)
layer1 = layer(df, x=:K, y=:Vmax)
y_VK = exp(coef(ols)[1]).*df.K.^(coef(ols)[2])
layer2 = layer(x=df.K, y=y_VK)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

ols = lm(@formula(log(VK) ~ log(Vc)), df)
layer1 = layer(df, x=:Vc, y=:VK)
y_V_K = exp(coef(ols)[1]).*df.Vc.^(coef(ols)[2])
layer2 = layer(x=df.Vc, y=y_V_K)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)


ols = lm(@formula(log2(VDNA/1e6) ~ log2(rrna)), df)
layer1 = layer(df, x=:VDNA, y=:rrna)
y_rnna = exp(coef(ols)[1]).*df.VDNA.^(coef(ols)[2])
layer2 = layer(x=df.VDNA, y=y_rnna)
plot(layer1, layer2, Scale.x_log2, Scale.y_log2)
