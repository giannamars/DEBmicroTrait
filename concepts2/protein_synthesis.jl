using DEBmicroTrait
using CSV, DataFrames
using Gadfly, Cairo, Fontconfig
using GLM

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/dethlefsen_schmidt_2008.csv", DataFrame, missingstring="N/A")
df_out = DataFrame()
df_terra = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/all_public_terrestrial.csv", DataFrame, missingstring="N/A")
df_terra_out = DataFrame()

##################################### Dethlefsen & Schmidt (2008): translational_power = μ*P/R
# Protein data
d_p = 1.37e6
P = df.Protein_Quantity_Measurement
V_P = P/d_p
df_out.V_P = V_P
# RNA data
d_r = 1.79e6
R = df.RNA_Quantity_Measurement
V_R = R/d_r
df_out.V_R = V_R
# Growth rate data
gmax = df.Specific_Growth_Rate_Measurement
df_out.gmax = gmax
# Cell size data
V_cell = df.Size_Measurement
df_out.V_cell = V_cell
# Cell size prediction
d_DNA = 2e6
DNA = df.DNA_Quantity_Measurement
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
# Protein/RNA prediction
V_P_R = V_P./V_R
df_out.V_P_R = V_P_R
V_P_R_model = V_P_model./V_R_model
df_out.V_P_R_model = V_P_R_model
# tRNA predictions
V_tRNA_model = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell, gmax)
df_out.tRNA = V_tRNA_model
V_tR_model = V_R_model .+ V_tRNA_model
df_out.V_tR_model = V_tR_model
V_P_tR_model = V_P_model./V_tR_model
df_out.V_P_tR_model = V_P_tR_model
# k_E predictions
k_E = gmax.*V_P_R
df_out.k_E = k_E
#k_E_model = gmax.*V_P_R_model
k_E_model = DEBmicroTrait.translation_power(V_P_model, V_R_model, log(2)./gmax)
df_out.k_E_model = k_E_model
k_E_model_t = DEBmicroTrait.translation_power(V_P_model, V_tR_model, log(2)./gmax)
df_out.k_E_model_t = k_E_model_t

# Regressions
linearRegressor = lm(@formula(log(k_E_model) ~ log(k_E)), df_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))
layer1 = layer(x=k_E, y=k_E_model, Geom.point)
layer2 = layer(x=k_E, y=linearFit, Geom.line)
plt = plot(layer1, layer2, Scale.x_log10, Scale.y_log10, Guide.xlabel("Measured translation rate [1/h]"), Guide.ylabel("Predicted translation rate [1/h]"))
plt |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/translation_rate_fit.png", dpi=300)

maximum(k_E_model)
log(2)/2.76
# terrestrial predictions
V_cell_model_terra = DEBmicroTrait.genome_size_to_cell_volume(df_terra.genomesize)
gmax_terra = log(2)./df_terra.mingentime
rrn_copies_terra = DEBmicroTrait.genome_size_to_rRNA_copy_number(df_terra.genomesize)
V_R_model_terra = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cell_model_terra, gmax_terra)
V_tR_model_terra = DEBmicroTrait.cell_volume_to_tRNA_volume(V_cell_model_terra, gmax_terra)
V_P_model_terra = DEBmicroTrait.cell_volume_to_protein_volume(V_cell_model_terra)
#k_E_terra = gmax_terra.*V_P_model_terra./V_R_model_terra
#k_E_t_terra = gmax_terra.*V_P_model_terra./(V_R_model_terra.+V_tR_model_terra)
k_E_terra = DEBmicroTrait.translation_power(V_P_model_terra, V_R_model_terra, df_terra.mingentime)
k_E_t_terra = DEBmicroTrait.translation_power(V_P_model_terra, V_R_model_terra.+V_tR_model_terra, df_terra.mingentime)

df_terra_out.k_E = k_E_terra
df_terra_out.k_E_t = k_E_t_terra
df_terra_out.rrn = rrn_copies_terra
df_terra_out.mingt = df_terra.mingentime
df_terra_out.gmax = gmax_terra
df_terra_out.V_R = V_R_model_terra
df_terra_out.V_P = V_P_model_terra
df_terra_out.L_DNA = df_terra.genomesize
df_terra_out.V_cell = V_cell_model_terra
df_terra_out.V_t_R = V_tR_model_terra./V_R_model_terra

linearRegressor = lm(@formula(log(k_E) ~ log(rrn)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=rrn_copies_terra, y=k_E_terra, Geom.point)
layer2 = layer(x=rrn_copies_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

linearRegressor = lm(@formula(log(k_E) ~ log(gmax)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

plt = plot(x=gmax_terra, y=V_P_model_terra./V_R_model_terra, Geom.point, Scale.x_log10, Scale.y_log10, Guide.xlabel("Max. specific growth rate [1/h]"), Guide.ylabel("Translation rate [1/h]"))
plt |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/translation_rate_global.png", dpi=300)


linearRegressor = lm(@formula(log(k_E_t) ~ log(gmax)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=gmax_terra, y=k_E_t_terra, Geom.point)
layer2 = layer(x=gmax_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)


linearRegressor = lm(@formula(log(V_t_R) ~ log(gmax)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=gmax_terra, y=V_tR_model_terra./V_R_model_terra, Geom.point)
layer2 = layer(x=gmax_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)


linearRegressor = lm(@formula(log(k_E) ~ log(V_R)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=V_R_model_terra, y=k_E_terra, Geom.point)
layer2 = layer(x=V_R_model_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

linearRegressor = lm(@formula(log(k_E) ~ log(V_P)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=V_P_model_terra, y=k_E_terra, Geom.point)
layer2 = layer(x=V_P_model_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

linearRegressor = lm(@formula(log(k_E) ~ log(V_cell)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=V_cell_model_terra, y=k_E_terra, Geom.point)
layer2 = layer(x=V_cell_model_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)

linearRegressor = lm(@formula(log(gmax) ~ log(V_cell)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=V_cell_model_terra, y=gmax_terra, Geom.point)
layer2 = layer(x=V_cell_model_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)


linearRegressor = lm(@formula(log(V_R) ~ log(V_cell)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=V_cell_model_terra, y=V_R_model_terra, Geom.point)
layer2 = layer(x=V_cell_model_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)


linearRegressor = lm(@formula(log(V_R) ~ log(gmax)), df_terra_out)
r2(linearRegressor)
linearFit = exp.(predict(linearRegressor))

layer1 = layer(x=gmax_terra, y=V_R_model_terra, Geom.point)
layer2 = layer(x=gmax_terra, y=linearFit, Geom.line)
plot(layer1, layer2, Scale.x_log10, Scale.y_log10)
