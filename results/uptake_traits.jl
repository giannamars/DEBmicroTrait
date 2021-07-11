using DEBmicroTrait
using Gadfly, Cairo, Fontconfig
using CSV, DataFrames, Statistics
using GLM, HypothesisTests
using Roots

df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame, missingstring="N/A")
df_isolates_out = DataFrame()
#df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")
df_metabolites.Ontology[df_metabolites.Ontology.=="sugar"] .= "Sugars"
df_metabolites.Ontology[df_metabolites.Ontology.=="organic acid"] .= "Organic acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="amino acid"] .= "Amino acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="fatty acid"] .= "Fatty acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="nucleotide"] .= "Nucleotides"
df_metabolites.Ontology[df_metabolites.Ontology.=="auxin"] .= "Auxins"

V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

################################################################################
y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C
y_EM            = ones(size(V_cell,1))

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

N_SB = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax = @. 180.0*60^2*N_SB.*N_C
K_D = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, df_metabolites.diffusivity)

df_isolates_out.Vmax = vec(Vmax)
df_isolates_out.K = vec(K_D)
df_isolates_out.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))

response = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_isolates_out.response = vec(response)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/uptake_traits_avena_exudation.csv", df_isolates_out)
#################################################################################


df_p = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_isolates_out)
kw_vmax_sugars = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_sugars = KruskalWallisTest(df_n.K, df_p.K)
#
df_p = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_isolates_out)
kw_vmax_organics = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_organics = KruskalWallisTest(df_n.K, df_p.K)
#
df_p = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_isolates_out)
kw_vmax_aminos = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_aminos = KruskalWallisTest(df_n.K, df_p.K)
#
df_p = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_isolates_out)
kw_vmax_fattys = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_fattys = KruskalWallisTest(df_n.K, df_p.K)
#
df_p = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_isolates_out)
kw_vmax_nucleos = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_nucleos = KruskalWallisTest(df_n.K, df_p.K)
#
df_p = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_isolates_out)
df_n = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_isolates_out)
kw_vmax_auxins = KruskalWallisTest(df_n.Vmax, df_p.Vmax)
kw_K_auxins = KruskalWallisTest(df_n.K, df_p.K)
#################################################################################
plot(x=vec(df_isolates_out.K), y=vec(df_isolates_out.Vmax), Geom.point, Scale.x_log10, Scale.y_log10)
