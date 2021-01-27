using DEBmicroTrait
using CSV, DataFrames
using Gadfly, Cairo, Fontconfig

df = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/isolates2traits.csv", DataFrame)

#Genome_size = convert(Array{Float64,1}, LinRange(2e6, 1e7, 100))
Genome_size = convert(Array{Float64,1}, df.Genome_size)
#rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
#gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = df.Min_gen_time
#Gram_stain = repeat(["+"], 100)
Gram_stain = convert(Array{String,1}, df.gram_stain)

p = IsolateComposition(Genome_size, Min_gen_time, Gram_stain)

Gadfly.push_theme(style(major_label_font_size=6pt))
pltCN = plot(x=Genome_size, y=p.CNP[1,:]./p.CNP[2,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("C:N"))
pltCP = plot(x=Genome_size, y=p.CNP[1,:]./p.CNP[3,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("C:P"))
pltNP = plot(x=Genome_size, y=p.CNP[2,:]./p.CNP[3,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("N:P"))
fig1  = hstack(pltCN, pltCP, pltNP)
fig1 |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/biomass_stoichiometry.png", dpi=300)

df = DataFrame()
df.Genome_size = Genome_size
df.Min_gen_time = Min_gen_time
df.percent_C = p.CNP[1,:]
df.percent_N = p.CNP[2,:]
df.percent_P = p.CNP[3,:]
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/CNP_isolates.csv", df)

pltyECV = plot(x=Genome_size, y=p.YEV[1,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("y_ECV"))
pltyENV = plot(x=Genome_size, y=p.YEV[2,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("y_ENV"))
pltyEPV = plot(x=Genome_size, y=p.YEV[3,:], Geom.point, Guide.xlabel("Genome size [bp]"), Guide.ylabel("y_EPV"))
fig2  = hstack(pltyECV, pltyENV, pltyEPV)
fig2 |> PNG("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/stoichiometric_yields.png", dpi=300)
