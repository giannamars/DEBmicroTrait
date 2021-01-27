using DEBmicroTrait, Test

n_microbes = rand(1:10)
n_monomers = rand(1:10)
n_polymers = rand(1:10)
n_enzymes = rand(1:10)

## MetabolismC
k_E = 0.2*ones(n_microbes)
y_EV = 1.0*ones(n_microbes)
k_M = 0.1*ones(n_microbes)
y_EM = 1.0*ones(n_microbes)
α_X = 0.0*ones(n_microbes)
y_EX = 1.0*ones(n_microbes)

p_met = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)

# AssimilationC
D = rand(n_monomers)
V = rand(n_microbes)
J_DE = zeros(n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_monomers)
K            = rand(n_monomers,n_microbes)
N_SB         = rand(n_monomers,n_microbes)
y_DE         = rand(n_monomers)
N_C          = rand(n_monomers)

p_ass = AssimilationC(N_SB,K,y_DE,N_C)

# Turnover
γ_V_0 = ones(n_microbes)
γ_V_1 = ones(n_microbes)
γ_X   = ones(n_enzymes)
f_ED  = ones(n_monomers)

p_tu  = Turnover(γ_V_0,γ_V_1,γ_X,f_ED)

# Depolymerization
Χ_0 = rand(n_polymers)
V_E = rand(n_polymers)
α_kin = rand(n_polymers, n_enzymes)
K_EP = rand(n_polymers, n_enzymes)

p_dep = Depolymerization(Χ_0, V_E, α_kin, K_EP)


# Microbe
p = Params(p_met, p_ass, p_dep, p_tu)

DEBmicroTrait.test_params(p)


@test 1==1


###### MetabolismCN
k_E = 0.2*ones(n_microbes)
y_EV = fill(1.0, 2, n_microbes)
k_M = 0.1*ones(n_microbes)
y_EM = fill(1.0, 2, n_microbes)
κ = fill(1.0, 2, n_microbes)
y_EX = fill(1.0, 2, n_microbes)

p_met = MetabolismCN(k_E, y_EV, k_M, y_EM, κ, y_EX)

# AssimilationCN

D = rand(n_monomers)
V = rand(n_microbes)
J_DE = zeros(n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_monomers)
K            = rand(n_monomers,n_microbes)
N_SB         = rand(n_monomers,n_microbes)
y_DE         = rand(n_monomers)
N_C          = rand(n_monomers)
N_N          = rand(n_monomers)

p_ass = AssimilationCN(N_SB,K,y_DE,N_C,N_N)


# Turnover
γ_V_0 = ones(n_microbes)
γ_V_1 = ones(n_microbes)
γ_X   = ones(n_enzymes)
f_ED  = ones(n_monomers)

p_tu  = Turnover(γ_V_0,γ_V_1,γ_X,f_ED)

# Microbe
p = Params(p_met, p_ass, p_dep, p_tu)
@test 1==1

#### MetabolismCNP
k_E = 0.2*ones(n_microbes)
y_EV = fill(1.0, 3, n_microbes)
k_M = 0.1*ones(n_microbes)
y_EM = fill(1.0, 3, n_microbes)
κ = fill(1.0, 3, n_microbes)
y_EX = fill(1.0, 3, n_microbes)

p_met = MetabolismCNP(k_E, y_EV, k_M, y_EM, κ, y_EX)

# AssimilationCNP
N_P          = rand(n_monomers)
J_DE         = zeros(3, n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_monomers)
p_ass = AssimilationCNP(N_SB,K,y_DE,N_C,N_N,N_P)

# Microbe
p = Params(p_met, p_ass, p_tu)
@test 1==1
