using DEBmicroTrait
using CSV, DataFrames
using Roots

df = CSV.read("/Users/glmarschmann/Data/Zhen/IsogenieGenomes.ecosysguilds.xls - IsogenieGenomes.ecosysguilds.xls.csv", DataFrame)

Genome_size = convert(Array{Float64}, df.bp)
V_cell = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
Min_gen_time = df.minimum_generation_time_hours
Gram_stain = repeat(["[+]"],1529)
y_EM = ones(1529)
y_DE = [0.24]
N_C  = [6]


ρ_ps = zeros(size(V_cell,1))

for i in 1:size(V_cell,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, [V_cell[i]], [Min_gen_time[i]], [Gram_stain[i]], y_DE[1], N_C[1], [y_EM[i]])
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    ρ_ps[i] = ρ_p[1]
end

D_S  = [1e-10]
K_D  = DEBmicroTrait.specific_reference_affinity(V_cell, reshape(ρ_ps, 1529,1), D_S)


using Gadfly

plot(x=K_D, Geom.histogram)
