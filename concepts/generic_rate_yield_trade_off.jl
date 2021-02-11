using DEBmicroTrait
using Roots
using Gadfly, Cairo, Fontconfig
using JLD
using CSV, DataFrames, Statistics


df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")


Genome_size = convert(Array{Float64,1}, LinRange(2e6, 9e6, 100))
V_cs = DEBmicroTrait.genome_size_to_cell_volume(Genome_size)
rrn_copies = DEBmicroTrait.genome_size_to_rRNA_copy_number(Genome_size)
gmax = DEBmicroTrait.gmax_regression(rrn_copies)
Min_gen_time = log(2)./gmax
Gram_stain = repeat(["+"], size(Genome_size,1))
α = fill(0.0, size(Genome_size,1))

V_p = DEBmicroTrait.cell_volume_to_protein_volume(V_cs)
V_r = DEBmicroTrait.cell_volume_to_ribosome_volume(V_cs, gmax)
k_E = DEBmicroTrait.translation_power(V_p, V_r, Min_gen_time)
k_M = DEBmicroTrait.cell_volume_to_specific_maintenance_rate(V_cs, Min_gen_time, Gram_stain)
y_EV = DEBmicroTrait.relative_translation_efficiency_regression(rrn_copies)


yield = zeros(size(V_cs,1), 5109)
rate = zeros(size(V_cs,1), 5109)
ρs = zeros(size(V_cs,1))
mEs = zeros(size(V_cs,1), 8424)

i=56
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
mEs[i,:] = DEBmicroTrait.reserve_density_curve(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_56 = yield[i,:]
rate_56 = rate[i,:]
mE_56 = mEs[i,:]

plot(x=yield_56, y=mE_56, Geom.point, Guide.xlabel("yield"), Guide.ylabel("reserve density"), Scale.y_log10, Scale.x_log10)
plot(x=mE_56, y=rate_56, Geom.point, Guide.xlabel("reserve density"), Guide.ylabel("Rate"), Scale.y_log10, Scale.x_log10)
plot(p1,p,Scale.y_log10)


j_EC  = mE_56.*(k_E[56] .- rate_56)
j_E = mE_56.*(k_E[56])
plot(x=rate_56, y=j_EC)
j_EM  = k_M[56].*ones(1)
jEM   = min.(j_EC, j_EM)
j_EG  = (j_EC - jEM)
plot(x=rate_56, y=j_EG)
rG_CO2  = (1 .- 1 ./y_EV[56]).*j_EG
p1=layer(x=rate_56, y=j_EG, color=["red"])
p2=layer(x=rate_56, y=0.4*j_EC, color=["blue"])
plot(p1,p2, Scale.y_log10)
y_VM  = ones(1)./y_EV[56]
jVM   = (j_EM .- jEM).*y_VM./ones(1)
rM_CO2  = @. (jEM + jVM)
p1=layer(x=mE_56, y=rate_56, color=["red"])
p2=layer(x=mE_56, y=yield_56, color=["blue"])
p3=layer(x=rate_56, y=rM_CO2, color=["green"])
plot(p1, Scale.y_log10, Scale.x_log10)

i=50
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_50 = yield[i,:]
rate_50 = rate[i,:]

yield = zeros(size(V_cs,1), 5109)
rate = zeros(size(V_cs,1), 5109)
i=45
y_DE = df_metabolites.y_DE[1]
N_C = df_metabolites.N_C[1]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
yield_45_glucose = yield[i,:]
rate_45_glucose = rate[i,:]

yield = zeros(size(V_cs,1), 5109)
rate = zeros(size(V_cs,1), 5109)
y_DE = df_metabolites.y_DE[3]
N_C = df_metabolites.N_C[3]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
yield_45_acetate = yield[i,:]
rate_45_acetate = rate[i,:]

yield = zeros(size(V_cs,1), 5109)
rate = zeros(size(V_cs,1), 5109)
y_DE = df_metabolites.y_DE[14]
N_C = df_metabolites.N_C[14]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
yield_45_lysine = yield[i,:]
rate_45_lysine = rate[i,:]

yield = zeros(size(V_cs,1), 5109)
rate = zeros(size(V_cs,1), 5109)
y_DE = df_metabolites.y_DE[22]
N_C = df_metabolites.N_C[22]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]], y_DE, N_C)
yield_45_xylose = yield[i,:]
rate_45_xylose = rate[i,:]

JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/generic_rate_yield_trade_off_mtm_3.jld",
         "yield_45_glucose", yield_45_glucose, "rate_45_glucose", rate_45_glucose,
         "yield_45_acetate", yield_45_acetate, "rate_45_acetate", rate_45_acetate,
         "yield_45_lysine", yield_45_lysine, "rate_45_lysine", rate_45_lysine,
         "yield_45_xylose", yield_45_xylose, "rate_45_xylose", rate_45_xylose)


i=35
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_35 = yield[i,:]
rate_35 = rate[i,:]

i=25
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_25 = yield[i,:]
rate_25 = rate[i,:]

i=15
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_15 = yield[i,:]
rate_15 = rate[i,:]

i=5
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
ρs[i] = ρ_p
mEs[i] = DEBmicroTrait.reserve_density(ρs[i], α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α[i], V_cs[i]*ones(1), [Min_gen_time[i]], [Gram_stain[i]])
yield_5 = yield[i,:]
rate_5 = rate[i,:]


rmax, rid = findmax(rate_5)
vcat(mEs[5], mEs[15], mEs[25], mEs[35], mEs[45])

Levin = zeros(7,2)
Ymax, Yid = findmax(yield_5)
Levin[1,1] = Ymax
Levin[1,2] = rate_5[Yid]
Ymax, Yid = findmax(yield_15)
Levin[2,1] = Ymax
Levin[2,2] = rate_15[Yid]
Ymax, Yid = findmax(yield_25)
Levin[3,1] = Ymax
Levin[3,2] = rate_25[Yid]
Ymax, Yid = findmax(yield_35)
Levin[4,1] = Ymax
Levin[4,2] = rate_35[Yid]
Ymax, Yid = findmax(yield_45)
Levin[5,1] = Ymax
Levin[5,2] = rate_45[Yid]
Ymax, Yid = findmax(yield_50)
Levin[6,1] = Ymax
Levin[6,2] = rate_50[Yid]
Ymax, Yid = findmax(yield_56)
Levin[7,1] = Ymax
Levin[7,2] = rate_56[Yid]


FCR = zeros(7)
rmax, rid = findmax(rate_5)
Ymax, Yid = findmax(yield_5)
FCR[1] = Ymax - yield_5[rid]
rmax, rid = findmax(rate_15)
Ymax, Yid = findmax(yield_15)
FCR[2] = Ymax - yield_15[rid]
rmax, rid = findmax(rate_25)
Ymax, Yid = findmax(yield_25)
FCR[3] = Ymax - yield_25[rid]
rmax, rid = findmax(rate_35)
Ymax, Yid = findmax(yield_35)
FCR[4] = Ymax - yield_35[rid]
rmax, rid = findmax(rate_45)
Ymax, Yid = findmax(yield_45)
FCR[5] = Ymax - yield_45[rid]
rmax, rid = findmax(rate_50)
Ymax, Yid = findmax(yield_50)
FCR[6] = Ymax - yield_50[rid]
rmax, rid = findmax(rate_56)
Ymax, Yid = findmax(yield_56)
FCR[7] = Ymax - yield_56[rid]

Vcs = zeros(7)
Vcs[1] = V_cs[5]
Vcs[2] = V_cs[15]
Vcs[3] = V_cs[25]
Vcs[4] = V_cs[35]
Vcs[5] = V_cs[45]
Vcs[6] = V_cs[50]
Vcs[7] = V_cs[56]

yield_5a = zeros(1, 6000)
rate_5a = zeros(1, 6000)
α = 1.0
V_c5 = [V_cs[5]]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[1]*ones(1), [Min_gen_time[5]], [Gram_stain[5]])
ρ_p = Roots.find_zero(find_ρ, 1.0)
yield_5a, rate_5a = DEBmicroTrait.rate_yield_trade_off(ρ_p, [α], V_cs[1]*ones(1), [Min_gen_time[5]], [Gram_stain[5]])





JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/generic_rate_yield_trade_off_7.jld",
         "yield_5", yield_5, "rate_5", rate_5, "yield_15", yield_15, "rate_15", rate_15,
         "yield_25", yield_25, "rate_25", rate_25, "yield_35", yield_35, "rate_35", rate_35,
         "yield_45", yield_45, "rate_45", rate_45, "yield_50", yield_50, "rate_50", rate_50,
         "yield_56", yield_56, "rate_56", rate_56, "levin", Levin, "fcr", FCR, "Vc", Vcs,
         "yield_5a", yield_5a, "rate_5a", rate_5a)



#α = 1.0
#V_cs = [1.05e-20]
#find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[1]*ones(1), Min_gen_time, Gram_stain)
#ρ_p = Roots.find_zero(find_ρ, 1.0)
#yield[6,:], rate[6,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α, V_cs[1]*ones(1), Min_gen_time, Gram_stain, r_rates)

# FCR = zeros(100)
# for i in 1:size(FCR,1)
#     rmax, rid = findmax(rate[i,:])
#     Ymax, Yid = findmax(yield[i,:])
#     FCR[i] = Ymax - yield[i,rid]
# end





#################
DEBmicroTrait.g_Φ_berg(0.0)
log10(DEBmicroTrait.activity_coefficient_berg(0.4))


using SymPy

P = symbols("P")

R_w = 0.138
R_e = 3.06
R_i = 5.57
P_w_0 = 0.363

a = R_e/R_w
b = R_i/R_e
f_a = (3*(1-P_w_0)/(a*(1+2*P_w_0))^2)*(a*(1+2*P_w_0) + (1-P_w_0)*(1+(1-P_w_0)/(3*a*P_w_0)))
P_w = P_w_0*(1-P-P*f_a)
S_3 = P_w + P
S_2 = 0.5*(a*P_w+P)/R_e
S_1 = 0.25*(a^2*P_w+P)/R_e^2

g_P = (-log(1-S_3) + 6*S_2*R_i/(1-S_3) + (12*S_1/(1-S_3) + 18*S_2^2/(1-S_3)^2)*R_i^2)
g_P_0 = g_P.subs(P,0.0)
gamma_E = exp(g_P - g_P_0)

a,b = symbols("a b")

K_m_i = 1/a*1/gamma_E + 1/b*exp(P)

dKdP = diff(K_m_i, P)

simplify(solveset(Eq(dKdP.subs(P,0.34), 0.0), a))
