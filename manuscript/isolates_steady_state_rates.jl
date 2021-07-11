using DEBmicroTrait
using CSV, DataFrames
using JLD

########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props1.csv", DataFrame, missingstring="N/A")
d_assimilation  = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_assimilation.jld")
d_enzymes       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates_enzymes.jld")
########################################

########################################
# traits
ρ_ps            = d_assimilation["rho"]
α               = d_enzymes["alpha"]
V_cell          = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
Min_gen_time    = df_isolates.Min_gen_time
rrn_copies      = convert(Array{Float64,1}, df_isolates.rRNA_genes)
y_EM            = ones(size(V_cell,1))
y_EX            = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)
y_DE            = df_metabolites.y_DE
N_C             = df_metabolites.N_C
########################################

########################################
# loop over metabolites
m_E_maxs         = zeros(size(N_C,1), size(V_cell,1))
r_ms             = zeros(size(N_C,1), size(V_cell,1))
x_ms             = zeros(size(N_C,1), size(V_cell,1))

for i in 1:size(N_C,1)
    m_E_maxs[i,:]    = DEBmicroTrait.steady_state_reserve_density(ρ_ps[i,:], V_cell, Min_gen_time, Gram_stain, y_DE[i], N_C[i])
    r_ms[i,:], x_ms[i,:]  = DEBmicroTrait.steady_state_rates(ρ_ps[i,:], α, V_cell, Min_gen_time, rrn_copies, Gram_stain, y_DE[i], N_C[i], y_EM, y_EX)
end
########################################

########################################
# IO
df_out = DataFrame()
df_out.r_m = vec(r_ms)
df_out.x_m = vec(x_ms)
df_out.m_E_max = vec(m_E_maxs)
response = Array{String}(undef, (size(y_DE,1),size(V_cell,1)))
for i in 1:size(V_cell,1)
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_out.response = vec(response)
df_out.ontology = repeat(df_metabolites.Ontology, size(V_cell,1))
########################################

########################################
# statistics
#sugars
df_sugars   = filter(x->(x.ontology.=="Sugars") , df_out)
df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_out)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_out)
rm_sugars    = KruskalWallisTest(df_p_sugars.r_m, df_n_sugars.r_m, df_u_sugars.r_m)
rm_sugars_p  = median(df_p_sugars.r_m)
rm_sugars_n  = median(df_n_sugars.r_m)
rm_sugars_u  = median(df_u_sugars.r_m)
xm_sugars    = KruskalWallisTest(df_p_sugars.x_m, df_n_sugars.x_m, df_u_sugars.x_m)
xm_sugars_p  = median(df_p_sugars.x_m)
xm_sugars_n  = median(df_n_sugars.x_m)
xm_sugars_u  = median(df_u_sugars.x_m)

#organic acids
df_organics   = filter(x->(x.ontology.=="Organic acids") , df_out)
df_p_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="positive") , df_out)
df_n_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="negative") , df_out)
df_u_organics = filter(x->(x.ontology.=="Organic acids"&&x.response.=="undefined") , df_out)
rm_organics    = KruskalWallisTest(df_p_organics.r_m, df_n_organics.r_m, df_u_organics.r_m)
rm_organics_p  = median(df_p_organics.r_m)
rm_organics_n  = median(df_n_organics.r_m)
rm_organics_u  = median(df_u_organics.r_m)
xm_organics    = KruskalWallisTest(df_p_organics.x_m, df_n_organics.x_m, df_u_organics.x_m)
xm_organics_p  = median(df_p_organics.x_m)
xm_organics_n  = median(df_n_organics.x_m)
xm_organics_u  = median(df_u_organics.x_m)

#amino acids
df_aminos   = filter(x->(x.ontology.=="Amino acids") , df_out)
df_p_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="positive") , df_out)
df_n_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="negative") , df_out)
df_u_aminos = filter(x->(x.ontology.=="Amino acids"&&x.response.=="undefined") , df_out)
rm_aminos    = KruskalWallisTest(df_p_aminos.r_m, df_n_aminos.r_m, df_u_aminos.r_m)
rm_aminos_p  = median(df_p_aminos.r_m)
rm_aminos_n  = median(df_n_aminos.r_m)
rm_aminos_u  = median(df_u_aminos.r_m)
xm_aminos    = KruskalWallisTest(df_p_aminos.x_m, df_n_aminos.x_m, df_u_aminos.x_m)
xm_aminos_p  = median(df_p_aminos.x_m)
xm_aminos_n  = median(df_n_aminos.x_m)
xm_aminos_u  = median(df_u_aminos.x_m)

#fatty acids
df_fattys   = filter(x->(x.ontology.=="Fatty acids") , df_out)
df_p_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="positive") , df_out)
df_n_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="negative") , df_out)
df_u_fattys = filter(x->(x.ontology.=="Fatty acids"&&x.response.=="undefined") , df_out)
rm_fattys    = KruskalWallisTest(df_p_fattys.r_m, df_n_fattys.r_m, df_u_fattys.r_m)
rm_fattys_p  = median(df_p_fattys.r_m)
rm_fattys_n  = median(df_n_fattys.r_m)
rm_fattys_u  = median(df_u_fattys.r_m)
xm_fattys    = KruskalWallisTest(df_p_fattys.x_m, df_n_fattys.x_m, df_u_fattys.x_m)
xm_fattys_p  = median(df_p_fattys.x_m)
xm_fattys_n  = median(df_n_fattys.x_m)
xm_fattys_u  = median(df_u_fattys.x_m)

#nucleotides
df_out_nucleos = filter(:r_m => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df_out)
df_out_nucleos = filter(:x_m => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df_out)
df_nucleos   = filter(x->(x.ontology.=="Nucleotides") , df_out_nucleos)
df_p_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="positive") , df_out_nucleos)
df_n_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="negative") , df_out_nucleos)
df_u_nucleos = filter(x->(x.ontology.=="Nucleotides"&&x.response.=="undefined") , df_out_nucleos)
rm_nucleos    = KruskalWallisTest(df_p_nucleos.r_m, df_n_nucleos.r_m, df_u_nucleos.r_m)
rm_nucleos_p  = median(df_p_nucleos.r_m)
rm_nucleos_n  = median(df_n_nucleos.r_m)
rm_nucleos_u  = median(df_u_nucleos.r_m)
xm_nucleos    = KruskalWallisTest(df_p_nucleos.x_m, df_n_nucleos.x_m, df_u_nucleos.x_m)
xm_nucleos_p  = median(df_p_nucleos.x_m)
xm_nucleos_n  = median(df_n_nucleos.x_m)
xm_nucleos_u  = median(df_u_nucleos.x_m)

#auxins
df_out_auxins = filter(:r_m => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df_out)
df_out_auxins = filter(:x_m => x -> !any(f -> f(x), (ismissing, isnothing, isnan)), df_out)
df_auxins   = filter(x->(x.ontology.=="Auxins") , df_out_auxins)
df_p_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="positive") , df_out_auxins)
df_n_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="negative") , df_out_auxins)
df_u_auxins = filter(x->(x.ontology.=="Auxins"&&x.response.=="undefined") , df_out_auxins)
rm_auxins    = KruskalWallisTest(df_p_auxins.r_m, df_n_auxins.r_m, df_u_auxins.r_m)
rm_auxins_p  = median(df_p_auxins.r_m)
rm_auxins_n  = median(df_n_auxins.r_m)
rm_auxins_u  = median(df_u_auxins.r_m)
xm_auxins    = KruskalWallisTest(df_p_auxins.x_m, df_n_auxins.x_m, df_u_auxins.x_m)
xm_auxins_p  = median(df_p_auxins.x_m)
xm_auxins_n  = median(df_n_auxins.x_m)
xm_auxins_u  = median(df_u_auxins.x_m)


using Gadfly
using GLM

plot(x=df_nucleos.x_m, y=df_nucleos.r_m, Geom.point, Scale.x_log10, Scale.y_log10)
linearRegressor = lm(@formula(log(r_m) ~ log(x_m)), df_nucleos)
r2(linearRegressor)
