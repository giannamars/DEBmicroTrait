using DEBmicroTrait, Test

########################################
n_polymers = 1
n_monomers = 1
n_microbes = 1
n_enzymes  = 1
n_minerals = 1
p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
@test p_set.dim == n_polymers + 2*n_monomers + 3*n_microbes + 2*n_enzymes
########################################

########################################
# split_state_rhizo
u     = ones(p_set.dim)
P     = u[1:n_polymers]
D     = u[1+n_polymers:n_polymers+n_monomers]
E     = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
V     = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
X     = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
D_ads = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers]
X_ads = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers]
CO2   = u[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes]
@test size(P,1) == n_polymers
@test size(D,1) == n_monomers
@test size(E,1) == n_microbes
@test size(V,1) == n_microbes
@test size(X,1) == n_enzymes
@test size(D_ads,1) == n_monomers
@test size(X_ads,1) == n_enzymes
@test size(CO2,1) == n_microbes
########################################

########################################
# metabolism
k_E     = 0.2*ones(n_microbes)
y_EV    = 1.0*ones(n_microbes)
k_M     = 0.1*ones(n_microbes)
y_EM    = 1.0*ones(n_microbes)
α_X     = 0.0*ones(n_microbes)
y_EX    = 1.0*ones(n_microbes)
f_αX    = ones(n_enzymes)

p_met                     = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX, f_αX)
r                         = DEBmicroTrait.growth!(0.0*ones(1), p_met, E, V)
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p_met, E, V)
J_EX                      = DEBmicroTrait.enzyme_production!(x, p_met, V)
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
J_DE         = DEBmicroTrait.assimilation!(zeros(p_set.n_monomers,p_set.n_microbes+p_set.n_minerals), p_ass, D, V)
J_DE_CO2     = DEBmicroTrait.assimilation_production!(zeros(p_set.n_monomers,p_set.n_microbes+p_set.n_minerals), p_ass, D, V)
J_D, J_DM    = DEBmicroTrait.uptake!(zeros(p_set.n_monomers,p_set.n_microbes+p_set.n_minerals), p_ass, D, V)
########################################

########################################
# turnover
γ_V_0        = ones(n_microbes)
γ_V_1        = ones(n_microbes)
γ_X          = ones(n_enzymes)
γ_D_ads      = ones(n_monomers)
γ_X_ads      = ones(n_enzymes)
f_ED         = ones(n_monomers)
f_VD         = ones(n_monomers)
f_VP         = ones(n_polymers)
f_V          = ones(n_microbes)
f_XD         = ones(n_monomers)
f_XP         = ones(n_polymers)
f_X          = ones(n_enzymes)

p_turn       = Turnover(γ_V_0,γ_V_1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)
J_ED         = DEBmicroTrait.reserve_recycling!(zeros(p_set.n_monomers), p_turn, E)
J_X          = DEBmicroTrait.enzyme_decay!(zeros(p_set.n_enzymes), p_turn, X)
J_XD, J_XP   = DEBmicroTrait.enzyme_recycling!(zeros(p_set.n_monomers), p_turn, X)
J_V          = DEBmicroTrait.biomass_turnover!(zeros(p_set.n_microbes), p_turn, V)
J_VD, J_VP   = DEBmicroTrait.biomass_recycling!(zeros(p_set.n_monomers), p_turn, V)
J_E          = DEBmicroTrait.biomass_turnover!(zeros(p_set.n_microbes), p_turn, E)
J_D_ads      = DEBmicroTrait.adsorbed_monomer_decay!(zeros(p_set.n_monomers), p_turn, D_ads)
J_X_ads      = DEBmicroTrait.adsorbed_enzyme_decay!(zeros(p_set.n_enzymes), p_turn, X_ads)
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
J_P, J_XM    = DEBmicroTrait.depolymerization!(zeros(n_polymers), p_depoly, P, X)
########################################

########################################
# rhizo_model
du = zeros(p_set.dim)


dP = @. du[1:n_polymers] = - J_P + J_VP + J_XP
dD = @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_P - J_DM + J_ED + J_VD + J_XD + J_D_ads
dE = @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p_met.k_E - r)*E - J_E
dV = @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
dX = @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X - J_XM + J_X_ads
dD_ads = @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_DM - J_D_ads
dX_ads = @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers] = J_XM - J_X_ads
dCO2 = @. du[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

mass_balance = @. dP + dD + dE + dV + dX + dD_ads + dX_ads + dCO2
@test mass_balance ≈ [0.0] atol = 1e-10
########################################

########################################
# rhs function
p = Params(p_set, p_met, p_ass, p_depoly, p_turn)
du = zeros(p_set.dim)
du = DEBmicroTrait.rhizo_model!(du, u, p, 0.0)
@test sum(du) ≈ 0.0 atol = 1e-10
########################################