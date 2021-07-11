using DEBmicroTrait
using CSV, DataFrames
using Roots
using GLM
using Gadfly

################################################
# I/O
df_button = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/button_1998_2.csv", DataFrame, missingstring="N/A")
################################################

################################################
# substrate: thermodynamic props
N_metabolites      = 13
chemFormBiom       = [1, 1.8, 0.2, 0.5, 0, 0, 0]
dGcox              = zeros(N_metabolites)
dGcat              = zeros(N_metabolites)
dGAn               = zeros(N_metabolites)
λ_base             = zeros(N_metabolites)
y_DE               = zeros(N_metabolites)
N_C                = zeros(N_metabolites)

for i in 1:N_metabolites
    elementstring = df_button.Formula[i]
    N_C[i] = DEBmicroTrait.extract_composition(elementstring)[1]
    out = DEBmicroTrait.get_lambda(elementstring, chemFormBiom)
    dGcox[i] = out[2][3]
    dGcat[i] = out[2][5]
    dGAn[i]  = out[2][8]
    λ_base[i]     = out[1][1]
end

η                       = 0.43
df_button.y_DE          = @. (λ_base*η*dGcat)/dGcox
df_button.N_C           = N_C
df_button.D_S           = DEBmicroTrait.aqueous_diffusivity(df_button.Molecular_weight)
################################################

################################################
# Species: traits
V_cell = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_button.Genome_size))
gmax = DEBmicroTrait.gmax_regression(convert(Array{Float64,1}, df_button.rrn_copies))
Min_gen_time = log(2)./df_button.Vmax_h
#Min_gen_time = log(2)./gmax
Gram_stain = convert(Array{String,1}, df_button.Gram_stain)
y_EM = ones(N_metabolites)

ρ_ps            = zeros(size(V_cell,1))
K_D             = zeros(size(V_cell,1))
Vmax            = zeros(size(V_cell,1))

for i in 1:size(V_cell,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], convert(Array{Float64,1},[df_button.rrn_copies[i]]), y_DE[i], N_C[i], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps[i] = ρ_p[1].*100
    K_D_tmp  = DEBmicroTrait.specific_reference_affinity([V_cell[i]], reshape([ρ_ps[i]],1,1), [df_button.D_S[i]])
    K_D[i] = K_D_tmp[1]
    N_SB   = DEBmicroTrait.transporter_density_to_monomer_uptake_sites([V_cell[i]], reshape([ρ_ps[i]./100],1,1), [Min_gen_time[i]], [Gram_stain[i]])
    Vmax_tmp  = @. 180.0*60^2*N_SB.*N_C[i]
    Vmax[i] = Vmax_tmp[1]
end

df_button.Kpred = K_D
df_button.Km    = df_button.Km_mug_l.*1e-3./df_button.Molecular_weight
df_button.Vmax_pred  = Vmax

################################################
# statistics
linearRegressor = lm(@formula(log(Km) ~ log(Kpred)), df_button)
r2(linearRegressor)
# test for significant difference to 1-1 line
df_button.diff = df_button.Kpred-df_button.Km
linearRegressor = lm(@formula(diff ~ Kpred), df_button)
# RMSD
RMSD = sqrt(1/(N_metabolites-1)*sum((df_button.Kpred - df_button.Km).^2))
################################################
linearRegressor = lm(@formula(log(Vmax_h) ~ log(Vmax)), df_button)
r2(linearRegressor)
# test for significant difference to 1-1 line
df_button.diff = df_button.Vmax-df_button.Vmax_h
linearRegressor = lm(@formula(diff ~ Vmax), df_button)
# RMSD
RMSD = sqrt(1/(N_metabolites-1)*sum((df_button.Vmax- df_button.Vmax_h).^2))
################################################
# I/O
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/benchmark_assimilation.csv", df_button)
################################################

################################################
# plot
layer1 = layer(x=df_button.Kpred, y=df_button.Km, Geom.point)
layer2 = layer(x=df_button.Km, y=df_button.Km, Geom.line)
plt = plot(layer1, layer2, Guide.xlabel("Predicted K"), Guide.ylabel("Measured K"), Scale.x_log10, Scale.y_log10)
