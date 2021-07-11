using DEBmicroTrait, Test
using DifferentialEquations
using CSV, DataFrames

########################################
# setup
n_polymers = 1
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 1
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
########################################

########################################
# metabolism
k_E     = rand(n_microbes)
y_EV    = rand(n_microbes)
k_M     = rand(n_microbes)
y_EM    = rand(n_microbes)
α_X     = rand(n_microbes)
y_EX    = rand(n_microbes)
f_αX    = rand(n_enzymes)
p_met   = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX)
########################################

########################################
# assimilation
K_D_0        = rand(n_monomers,n_microbes)
K_M_0        = rand(n_monomers,n_minerals)
K_D          = hcat(K_D_0, K_M_0)
N_SB_D       = rand(n_monomers,n_microbes)
N_SB_M       = ones(n_monomers,n_minerals)
N_SB         = hcat(N_SB_D, N_SB_M)
N_C          = rand(n_monomers)
y_DE         = rand(n_monomers)
M            = rand(n_minerals)
p_ass        = AssimilationCM(N_SB,K_D,y_DE,N_C,M)
########################################

########################################
# turnover
γ_V_0        = rand(n_microbes)
γ_V_1        = 1000.0.*ones(n_microbes)
γ_X          = rand(n_enzymes)
γ_D_ads      = rand(n_monomers)
γ_X_ads      = rand(n_enzymes)
f_ED         = ones(n_monomers)
f_VD         = ones(n_monomers)
f_VP         = ones(n_polymers)
f_V          = rand(n_microbes)
f_XD         = ones(n_monomers)
f_XP         = ones(n_polymers)
f_X          = rand(n_enzymes)
p_turn       = Turnover(γ_V_0,γ_V_1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
########################################

########################################
# depolymerization
Χ_0          = rand(n_polymers)
V_E          = rand(n_polymers)
α_kin_P      = rand(n_polymers, n_enzymes)
K_EP_P       = rand(n_polymers, n_enzymes)
M            = rand(n_minerals)
α_kin_M      = ones(n_minerals, n_enzymes)
K_MX         = rand(n_minerals, n_enzymes)
α_kin        = vcat(α_kin_P, α_kin_M)
K_EP         = vcat(K_EP_P,K_MX)
p_depoly     = DepolymerizationM(Χ_0, V_E, α_kin, K_EP, M)
########################################

########################################
# run model
p        = Params(p_set, p_met, p_ass, p_depoly, p_turn)
u0       = 1e-3*ones(p_set.dim)
tspan    = (0.0, 100000.0)
prob     = ODEProblem(DEBmicroTrait.rhizo_model_linear!, u0, tspan, p)
sol      = solve(prob, alg_hints=[:stiff])
