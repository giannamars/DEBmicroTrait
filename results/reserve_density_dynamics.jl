using DEBmicroTrait
using Gadfly, Cairo, Fontconfig
using CSV, DataFrames, Statistics
using GLM, HypothesisTests

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/dethlefsen_schmidt_2008.csv", DataFrame, missingstring="N/A")
gdf = groupby(df, :Species_Name)
df = combine(gdf, :Protein_Quantity_Measurement => mean, :Size_Measurement => mean, :DNA_Quantity_Measurement => mean, :RNA_Quantity_Measurement => mean, :Specific_Growth_Rate_Measurement => mean, )
df.Species_Name
df_out = DataFrame()
df_isolates = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame, missingstring="N/A")
df_isolates_out = DataFrame()
df_terra = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/all_public_terrestrial.csv", DataFrame, missingstring="N/A")
df_terra_out = DataFrame()

##################################### Dethlefsen & Schmidt (2008): translational_power = μ*P/R
# Protein data
d_p = 1.37e6
P = df.Protein_Quantity_Measurement_mean
V_P = P/d_p
df_out.V_P = V_P
# RNA data
d_r = 1.79e6
R = df.RNA_Quantity_Measurement_mean
V_R = R/d_r
df_out.V_R = V_R
# Growth rate data
gmax = df.Specific_Growth_Rate_Measurement_mean
df_out.gmax = gmax
# Cell size data
V_cell = df.Size_Measurement_mean
df_out.V_cell = V_cell
# k_E calculation
k_E = gmax.*V_P./V_R
df_out.k_E = k_E
##################################### Kempes et al. (2016): allometric predictions
# Cell size prediction
d_DNA = 2e6
DNA = df.DNA_Quantity_Measurement_mean
V_DNA = DNA/d_DNA
v_N = 1.47e-27
L_DNA = V_DNA/v_N
V_cell_model = DEBmicroTrait.genome_size_to_cell_volume(L_DNA)
df_out.V_cell_model = V_cell_model
# Protein prediction
V_P_model = DEBmicroTrait.cell_volume_to_protein_volume(V_cell)
df_out.V_P = V_P
df_out.V_P_model = V_P_model
# RNA prediction
V_R_model = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell, gmax)
df_out.V_R = V_R
df_out.V_R_model = V_R_model
# tRNA predictions
V_tRNA_model = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell, gmax)
df_out.tRNA = V_tRNA_model
V_tR_model = V_R_model .+ V_tRNA_model
df_out.V_tR_model = V_tR_model
# mRNA predictions
V_mRNA_model = DEBmicroTrait.cell_volume_to_mRNA_volume(V_cell, gmax)
df_out.mRNA = V_mRNA_model
V_mtR_model = V_R_model
df_out.V_mtR_model = V_mtR_model
# k_E prediction
k_E_model= DEBmicroTrait.translation_power(V_P_model, V_mtR_model, log(2)./gmax)
df_out.k_E_model = k_E_model
##################################### Comparison to culture data
linearRegressor = lm(@formula(log(k_E_model) ~ log(k_E)), df_out)
r2(linearRegressor)
layer1 = layer(x=k_E_model, y=k_E, Geom.point)
layer2 = layer(x=k_E, y=k_E, Geom.line)
plt = plot(layer1, layer2, Guide.xlabel("Predicted translation rate [1/h]"), Guide.ylabel("Measured translation rate [1/h]"))
# test for significant difference to 1-1 line
df_out.diff = k_E_model-k_E
linearRegressor = lm(@formula(diff ~ k_E_model), df_out)
#plt |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/translation_rate_fit.png", dpi=300)
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/reserve_density_dynamics_benchmark.csv", df_out)

##################################### Prediction for Hopland isolates
V_cell_isolates = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_isolates.Genome_size))
gmax_isolates = log(2)./df_isolates.Min_gen_time
V_P_isolates = DEBmicroTrait.cell_volume_to_protein_volume(V_cell_isolates)
V_mtR_isolates = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell_isolates, gmax_isolates)
k_E_isolates= DEBmicroTrait.translation_power(V_P_isolates, V_mtR_isolates, log(2)./gmax_isolates)
#
y_EV_isolates = DEBmicroTrait.relative_translation_efficiency(V_P_isolates, V_mtR_isolates)
#plt = plot(x=y_EV_isolates, y=k_E_isolates, Geom.point, Scale.x_log10, Scale.y_log10)
y_EV_isolates_regr = DEBmicroTrait.relative_translation_efficiency_regression(convert(Array{Float64,1},df_isolates.rRNA_genes))
#
plt = plot(x=y_EV_isolates_regr, y=k_E_isolates, Geom.point, Scale.x_log10, Scale.y_log10)
#
df_isolates_out.Vcell = V_cell_isolates
df_isolates_out.kE = k_E_isolates
df_isolates_out.yEV = y_EV_isolates
df_isolates_out.yEVr = y_EV_isolates_regr
df_isolates_out.response = df_isolates.Rhizosphere_response
df_isolates_out.gmax = gmax_isolates
df_isolates_out.VP = V_P_isolates
df_isolates_out.VR = V_mtR_isolates
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/reserve_density_dynamics_isolates.csv", df_isolates_out)

linearRegressor = lm(@formula(log(kE) ~ log(yEVr)), df_isolates_out)
r2(linearRegressor)


df_positive = df_isolates_out[df_isolates_out.response .== "positive", :]
df_negative = df_isolates_out[df_isolates_out.response .== "negative", :]
df_undefined = df_isolates_out[df_isolates_out.response .== "undefined", :]

# only driven by differences in gmax
kE_kw = KruskalWallisTest(df_positive.kE, df_negative.kE)
yEV_kw = KruskalWallisTest(df_positive.yEV, df_negative.yEV)

# how does Kempes equation scale with μ?!


##################################### Comparison across all soil genomes
df_terra = df_terra[df_terra[!,:mingentime] .>= 0.3, :]
V_cell_model_terra = DEBmicroTrait.genome_size_to_cell_volume(df_terra.genomesize)
gmax_terra = log(2)./df_terra.mingentime
rrn_copies_terra = DEBmicroTrait.genome_size_to_rRNA_copy_number(df_terra.genomesize)
V_R_model_terra = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell_model_terra, gmax_terra)
V_tR_model_terra = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell_model_terra, gmax_terra)
V_mR_model_terra = DEBmicroTrait.cell_volume_to_mRNA_volume(V_cell_model_terra, gmax_terra)
V_P_model_terra = DEBmicroTrait.cell_volume_to_protein_volume(V_cell_model_terra)
V_mtR_model_terra = @. V_R_model_terra + V_tR_model_terra + V_mR_model_terra
k_E_terra = DEBmicroTrait.translation_power(V_P_model_terra, V_mtR_model_terra, df_terra.mingentime)

df_terra_out.k_E = k_E_terra
df_terra_out.rrn = rrn_copies_terra
df_terra_out.mingt = df_terra.mingentime
df_terra_out.gmax = gmax_terra
df_terra_out.L_DNA = df_terra.genomesize
df_terra_out.V_cell = V_cell_model_terra

linearRegressor = lm(@formula(log(k_E) ~ log(mingt)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))
plt = plot(x=df_terra.mingentime, y=k_E_terra, Geom.point, Scale.x_log10, Scale.y_log10, Guide.xlabel("Max. specific growth rate [1/h]"), Guide.ylabel("Translation rate [1/h]"))
plt |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/translation_rate_global.png", dpi=300)

##################################### Isolates
V_cell_model_isolates = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1},df_isolates.Genome_size))
gmax_isolates = log(2)./df_isolates.Min_gen_time
rrn_copies_isolates = df_isolates.rRNA_genes
V_R_model_isolates = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell_model_isolates, gmax_isolates)
V_tR_model_isolates = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell_model_isolates, gmax_isolates)
V_mR_model_isolates = DEBmicroTrait.cell_volume_to_mRNA_volume(V_cell_model_isolates, gmax_isolates)
V_P_model_isolates = DEBmicroTrait.cell_volume_to_protein_volume(V_cell_model_isolates)
V_mtR_model_isolates = @. V_R_model_isolates + V_tR_model_isolates + V_mR_model_isolates
k_E_isolates = DEBmicroTrait.translation_power(V_P_model_isolates, V_mtR_model_isolates, df_isolates.Min_gen_time)

df_isolates_out.k_E = k_E_isolates
df_isolates_out.rrn = rrn_copies_isolates
df_isolates_out.mingt = df_isolates.Min_gen_time
df_isolates_out.gmax = gmax_isolates
df_isolates_out.L_DNA = df_isolates.Genome_size
df_isolates_out.V_cell = V_cell_model_isolates
df_isolates_out.response = df_isolates.Rhizosphere_response

linearRegressor = lm(@formula(log(k_E) ~ log(mingt)), df_isolates_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))
plt = plot(x=df_isolates.Min_gen_time, y=k_E_isolates, Geom.point, Scale.x_log10, Scale.y_log10, Guide.xlabel("Max. specific growth rate [1/h]"), Guide.ylabel("Translation rate [1/h]"))
plt |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/translation_rate_global.png", dpi=300)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/results/files/reserve_density_dynamics_isolates.csv", df_isolates_out)
