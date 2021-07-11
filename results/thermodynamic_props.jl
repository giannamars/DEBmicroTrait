using DEBmicroTrait
using CSV, DataFrames, Statistics

#df_metabolites     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/files/Saifuddin/Metabolite_Data_Saifuddin.csv", DataFrame, missingstring="N/A")
df_metabolites      = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/avena_exudation.csv", DataFrame, missingstring="N/A")
df_metabolites_out  = df_metabolites

chemFormBiom       = [1, 1.8, 0.2, 0.5, 0, 0, 0]

N_metabolites      = 84
dGcox              = zeros(N_metabolites)
dGcat              = zeros(N_metabolites)
dGAn               = zeros(N_metabolites)
λ_base             = zeros(N_metabolites)
y_DE               = zeros(N_metabolites)
N_C                = zeros(N_metabolites)


for i in 1:N_metabolites
    elementstring = df_metabolites.Formula[i]
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
    out = DEBmicroTrait.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
    dGcat[i] = out[2][5]
    dGAn[i]  = out[2][8]
    λ_base[i]     = out[1][1]
end

η                       = 0.43
y_DE                    = @. (λ_base*η*dGcat)/dGcox

df_metabolites.Ontology[df_metabolites.Ontology.=="sugar"] .= "Sugars"
df_metabolites.Ontology[df_metabolites.Ontology.=="organic acid"] .= "Organic acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="amino acid"] .= "Amino acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="fatty acid"] .= "Fatty acids"
df_metabolites.Ontology[df_metabolites.Ontology.=="nucleotide"] .= "Nucleotides"
df_metabolites.Ontology[df_metabolites.Ontology.=="auxin"] .= "Auxins"


df_metabolites_out.y_DE = y_DE
df_metabolites_out.N_C = N_C
df_metabolites_out.diffusivity = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_weight)
df_metabolites_out.lambda = λ_base

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/avena_exudation_props1.csv", df_metabolites_out)

df_sugars = filter(x->(x.Ontology.=="Nucleotides") , df_metabolites_out)
df_organics = filter(x->(x.Ontology.=="Amino acids") , df_metabolites_out)

kw_vmax_sugars = KruskalWallisTest(df_sugars.lambda, df_organics.lambda)
kw_K_sugars = KruskalWallisTest(df_n.K, df_p.K)
