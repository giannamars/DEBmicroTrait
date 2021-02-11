using DEBmicroTrait
using CSV, DataFrames, Statistics


df_metabolites     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")
chemFormBiom       = [1, 1.8, 0.2, 0.5, 0, 0, 0]

N_metabolites      = 25
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

df_metabolites.y_DE = y_DE
df_metabolites.N_C = N_C
df_metabolites.diffusivity = DEBmicroTrait.aqueous_diffusivity(df_metabolites.Molecular_Weight)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", df_metabolites)
