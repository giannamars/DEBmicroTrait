using DEBmicroTrait
using CSV, DataFrames, Statistics
using Roots
using HypothesisTests
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")
########################################

########################################
# isolate traits
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
Gram_stain      = convert(Array{String,1}, df_isolates.gram_stain)
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)

z_sugars        = reshape(convert(Array{Float64,1}, df_isolates.z_sugars./df_isolates.Genome_size*1e6),1,39)
z_organics      = reshape(convert(Array{Float64,1}, df_isolates.z_organic_acids./df_isolates.Genome_size*1e6),1,39)
z_aminos        = reshape(convert(Array{Float64,1}, df_isolates.z_amino_acids./df_isolates.Genome_size*1e6), 1,39)
z_fattys        = reshape(convert(Array{Float64,1}, df_isolates.z_fatty_acids./df_isolates.Genome_size*1e6),1,39)
z_nucleos       = reshape(convert(Array{Float64,1}, df_isolates.z_nucleotides./df_isolates.Genome_size*1e6),1,39)
z_auxins        = reshape(convert(Array{Float64,1}, df_isolates.z_auxins./df_isolates.Genome_size*1e6),1,39)
genome_distr    = vcat(z_sugars, z_organics, z_aminos, z_fattys, z_nucleos, z_auxins)

y_EM            = ones(size(V_cell,1))
########################################

########################################
# monomer properties
y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C
########################################

########################################
# calc transporter density
ρ_ps            = zeros(size(y_DE,1), size(V_cell,1))

for j in 1:size(y_DE,1)
    if df_metabolites.Ontology[j] == "Sugars"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[1]
        end
    elseif df_metabolites.Ontology[j] == "Organic acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[2]
        end
    elseif df_metabolites.Ontology[j] == "Amino acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[3]
        end
    elseif df_metabolites.Ontology[j] == "Fatty acids"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[4]
        end
    elseif df_metabolites.Ontology[j] == "Nucleotides"
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], [rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[5]
        end
    else
        for i in 1:size(V_cell,1)
            #find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]],[rrn_copies[i]], y_DE[j], N_C[j], [y_EM[i]])
            find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[j], N_C[j], [y_EM[i]])
            ρ_p = Roots.find_zero(find_ρ, 1.0)
            closure = genome_distr[:,i]./sum(genome_distr[:,i])
            ρ_ps[j,i] = ρ_p[1].*closure[6]
        end
    end
end

ρ_ps[ρ_ps.==0.0] .= 1e-12
########################################

########################################
# uptake traits
N_SB              = DEBmicroTrait.transporter_density_to_monomer_uptake_sites(V_cell, ρ_ps, Min_gen_time, Gram_stain)
Vmax              = @. 180.0*60^2*N_SB.*N_C
K_D               = DEBmicroTrait.specific_reference_affinity(V_cell, ρ_ps, df_metabolites.diffusivity)
########################################

########################################
# I/O
df_assimilation = DataFrame()
df_assimilation.transporter_density = vec(ρ_ps)
df_assimilation.Vmax = vec(Vmax)
df_assimilation.KD = vec(K_D)
df_assimilation.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))
response = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_assimilation.response = vec(response)
genome_size = Array{Int}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     genome_size[:,i] .= df_isolates.Genome_size[i]
end
df_assimilation.genome_size = vec(genome_size)
min_gt = Array{Float64}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     min_gt[:,i] .= df_isolates.Min_gen_time[i]
end
df_assimilation.mingt = vec(min_gt)

#
save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_assimilation.jld", "rho", ρ_ps, "NSB", N_SB, "KD", K_D, "yEM", y_EM, "yDE", y_DE, "NC", N_C)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_assimilation.csv", df_assimilation)

########################################

########################################
# statistics
#sugars
df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_assimilation)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_assimilation)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_assimilation)
ρ_sugars    = KruskalWallisTest(df_p_sugars.transporter_density, df_n_sugars.transporter_density)
ρ_sugars_p  = median(df_p_sugars.transporter_density)
ρ_sugars_n  = median(df_n_sugars.transporter_density)
ρ_sugars_u  = median(df_u_sugars.transporter_density)
Vmax_sugars    = KruskalWallisTest(df_p_sugars.Vmax, df_n_sugars.Vmax)
Vmax_sugars_p  = median(df_p_sugars.Vmax)
Vmax_sugars_n  = median(df_n_sugars.Vmax)
Vmax_sugars_u  = median(df_u_sugars.Vmax)
K_sugars    = KruskalWallisTest(df_p_sugars.KD, df_n_sugars.KD)
K_sugars_p  = median(df_p_sugars.KD)
K_sugars_n  = median(df_n_sugars.KD)
K_sugars_u  = median(df_u_sugars.KD)
#organic acids
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_assimilation)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_assimilation)
df_u_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="undefined") , df_assimilation)
ρ_organics    = KruskalWallisTest(df_p_organics.transporter_density, df_n_organics.transporter_density)
ρ_organics_p  = median(df_p_organics.transporter_density)
ρ_organics_n  = median(df_n_organics.transporter_density)
ρ_organics_u  = median(df_u_organics.transporter_density)
Vmax_organics    = KruskalWallisTest(df_p_organics.Vmax, df_n_organics.Vmax)
Vmax_organics_p  = median(df_p_organics.Vmax)
Vmax_organics_n  = median(df_n_organics.Vmax)
Vmax_organics_u  = median(df_u_organics.Vmax)
K_organics    = KruskalWallisTest(df_p_organics.KD, df_n_organics.KD, df_u_organics.KD)
K_organics_p  = median(df_p_organics.KD)
K_organics_n  = median(df_n_organics.KD)
K_organics_u  = median(df_u_organics.KD)
#amino acids
df_p_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_assimilation)
df_n_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_assimilation)
df_u_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="undefined") , df_assimilation)
ρ_aminos    = KruskalWallisTest(df_p_aminos.transporter_density, df_n_aminos.transporter_density,  df_u_aminos.transporter_density)
ρ_aminos_p  = median(df_p_aminos.transporter_density)
ρ_aminos_n  = median(df_n_aminos.transporter_density)
ρ_aminos_u  = median(df_u_aminos.transporter_density)
Vmax_aminos    = KruskalWallisTest(df_p_aminos.Vmax, df_n_aminos.Vmax)
Vmax_aminos_p  = median(df_p_aminos.Vmax)
Vmax_aminos_n  = median(df_n_aminos.Vmax)
Vmax_aminos_u  = median(df_u_aminos.Vmax)
K_aminos    = KruskalWallisTest(df_p_aminos.KD, df_n_aminos.KD)
K_aminos_p  = median(df_p_aminos.KD)
K_aminos_n  = median(df_n_aminos.KD)
K_aminos_u  = median(df_u_aminos.KD)
#fatty acids
df_p_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_assimilation)
df_n_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_assimilation)
df_u_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="undefined") , df_assimilation)
ρ_fattys    = KruskalWallisTest(df_p_fattys.transporter_density, df_n_fattys.transporter_density, df_u_fattys.transporter_density)
ρ_fattys_p  = median(df_p_fattys.transporter_density)
ρ_fattys_n  = median(df_n_fattys.transporter_density)
ρ_fattys_u  = median(df_u_fattys.transporter_density)
Vmax_fattys    = KruskalWallisTest(df_p_fattys.Vmax, df_n_fattys.Vmax, df_u_fattys.Vmax)
Vmax_fattys_p  = median(df_p_fattys.Vmax)
Vmax_fattys_n  = median(df_n_fattys.Vmax)
Vmax_fattys_u  = median(df_u_fattys.Vmax)
K_fattys    = KruskalWallisTest(df_p_fattys.KD, df_n_fattys.KD, df_u_fattys.KD)
K_fattys_p  = median(df_p_fattys.KD)
K_fattys_n  = median(df_n_fattys.KD)
K_fattys_u  = median(df_u_fattys.KD)
#Nucleotides
df_p_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_assimilation)
df_n_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_assimilation)
df_u_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="undefined") , df_assimilation)
ρ_nucleos    = KruskalWallisTest(df_p_nucleos.transporter_density, df_n_nucleos.transporter_density, df_u_nucleos.transporter_density)
ρ_nucleos_p  = median(df_p_nucleos.transporter_density)
ρ_nucleos_n  = median(df_n_nucleos.transporter_density)
ρ_nucleos_u  = median(df_u_nucleos.transporter_density)
Vmax_nucleos    = KruskalWallisTest(df_p_nucleos.Vmax, df_n_nucleos.Vmax, df_u_nucleos.Vmax)
Vmax_nucleos_p  = median(df_p_nucleos.Vmax)
Vmax_nucleos_n  = median(df_n_nucleos.Vmax)
Vmax_nucleos_u  = median(df_u_nucleos.Vmax)
K_nucleos    = KruskalWallisTest(df_p_nucleos.KD, df_n_nucleos.KD, df_u_nucleos.KD)
K_nucleos_p  = median(df_p_nucleos.KD)
K_nucleos_n  = median(df_n_nucleos.KD)
K_nucleos_u  = median(df_u_nucleos.KD)
#Auxins
df_p_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_assimilation)
df_n_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_assimilation)
df_u_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="undefined") , df_assimilation)
ρ_auxins    = KruskalWallisTest(df_p_auxins.transporter_density, df_n_auxins.transporter_density, df_u_auxins.transporter_density)
ρ_auxins_p  = median(df_p_auxins.transporter_density)
ρ_auxins_n  = median(df_n_auxins.transporter_density)
ρ_auxins_u  = median(df_u_auxins.transporter_density)
Vmax_auxins    = KruskalWallisTest(df_p_auxins.Vmax, df_n_auxins.Vmax, df_u_auxins.Vmax)
Vmax_auxins_p  = median(df_p_auxins.Vmax)
Vmax_auxins_n  = median(df_n_auxins.Vmax)
Vmax_auxins_u  = median(df_u_auxins.Vmax)
K_auxins    = KruskalWallisTest(df_p_auxins.KD, df_n_auxins.KD, df_u_auxins.KD)
K_auxins_p  = median(df_p_auxins.KD)
K_auxins_n  = median(df_n_auxins.KD)
K_auxins_u  = median(df_u_auxins.KD)
########################################
