using DEBmicroTrait
using Gadfly, Cairo, Fontconfig
using CSV, DataFrames, Statistics
using GLM, HypothesisTests
using Roots

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame, missingstring="N/A")
df_isolates_out = DataFrame()
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")

V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)

y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C

y_EM            = ones(39)
gmax   = log(2)./Min_gen_time
V_p = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
V_r = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
y_EX = DEBmicroTrait.relative_translation_efficiency_regression(convert(Array{Float64},df_isolates.rRNA_genes))


z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

ρ_ps            = zeros(size(y_DE,1), size(V_cell,1))

for j in 1:size(y_DE,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[2]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[3]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[4]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[5]
        end
    else
        for i in 1:size(V_cell,1)
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[6]
        end
    end
end

ρ_ps[ρ_ps.==0.0] .= 1e-8


# ρ_ps            = zeros(size(V_cell,1))
# for i in 1:size(V_cell,1)
#     find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[1], N_C[1], [y_EM[i]])
#     ρ_p = Roots.find_zero(find_ρ, 1.0)
#     ρ_ps[i] = ρ_p[1]
# end

α  = DEBmicroTrait.constrain_enzyme_allocation(V_cell, Min_gen_time, Gram_stain, df_isolates.z_hydrolases./df_isolates.Genome_size*1e6)

rs  = zeros(size(y_DE,1), size(V_cell,1))
xs  = zeros(size(y_DE,1), size(V_cell,1))



r,x = DEBmicroTrait.steady_state_rates(ρ_ps[1,:], α, V_cell, Min_gen_time, Gram_stain, y_DE[1], N_C[1], y_EM, y_EX)

for j in 1:size(y_DE,1)
    r,x = DEBmicroTrait.steady_state_rates(ρ_ps[j,:], α, V_cell, Min_gen_time, Gram_stain, y_DE[j], N_C[j], y_EM, y_EX)
    rs[j,:] = r
    xs[j,:] = x
end

df_isolates_out.r = vec(rs)
df_isolates_out.x = vec(xs)
df_isolates_out.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))
response1 = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response1[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_isolates_out.response = vec(response1)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/transporter_enzyme_investment_1.csv", df_isolates_out)


yield, rate = DEBmicroTrait.steady_state_rate_yield(ρ_ps, α, V_cell, Min_gen_time, Gram_stain, y_DE[1], N_C[1], y_EM, y_EX)

plot(x=yield, y=rate, Geom.point, Scale.x_log10, Scale.y_log10)

FCR = zeros(size(V_cell,1))
for i in 1:size(V_cell,1)
    yield, rate = DEBmicroTrait.rate_yield_trade_off(ρ_ps[i]*ones(1,1), [α[i]], [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[1], N_C[1], [y_EM[i]], [y_EX[i]])
    FCR[i] = DEBmicroTrait.functional_control_region(rate, yield)
end

plot(x=V_cell, y=FCR, Geom.point, Scale.x_log10, Scale.y_log10)
plot(x=yield, y=rate, Geom.point, Scale.x_log10, Scale.y_log10)

df_isolates_out.Vcell = V_cell
df_isolates_out.FCR = FCR
linearRegressor = lm(@formula(log(FCR) ~ log(Vcell)), df_isolates_out)
r2(linearRegressor)
