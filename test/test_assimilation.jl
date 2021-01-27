using DEBmicroTrait, Test

n_substrates = rand(1:10)
n_microbes   = rand(1:10)

D = rand(n_substrates)
V = rand(n_microbes)
J_DE = zeros(n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_substrates)
K            = rand(n_substrates,n_microbes)
N_SB         = rand(n_substrates,n_microbes)
y_DE         = rand(n_substrates)
N_C          = rand(n_substrates)

p = AssimilationC(N_SB,K,y_DE,N_C)
J_DE = DEBmicroTrait.assimilation!(J_DE, p, D, V)
J_DE_CO2 = DEBmicroTrait.assimilation_production!(J_DE_CO2, p, D, V)
J_D = DEBmicroTrait.uptake!(J_D, p, D, V)
@test size(J_DE,1) == n_microbes
@test size(J_DE_CO2,1) == n_microbes
@test size(J_D,1) == n_substrates

n_reserves   = 2
N_N          = rand(n_substrates)
J_DE         = zeros(n_reserves, n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_substrates)
p = AssimilationCN(N_SB,K,y_DE,N_C,N_N)
J_DE = DEBmicroTrait.assimilation!(J_DE, p, D, V)
J_DE_CO2 = DEBmicroTrait.assimilation_production!(J_DE_CO2, p, D, V)
J_D = DEBmicroTrait.uptake!(J_D, p, D, V)
@test size(J_DE) == (n_reserves, n_microbes)
@test size(J_DE_CO2,1) == n_microbes
@test size(J_D,1) == n_substrates

n_reserves   = 3
N_P          = rand(n_substrates)
J_DE         = zeros(n_reserves, n_microbes)
J_DE_CO2 = zeros(n_microbes)
J_D = zeros(n_substrates)
p = AssimilationCNP(N_SB,K,y_DE,N_C,N_N,N_P)
J_DE = DEBmicroTrait.assimilation!(J_DE, p, D, V)
J_DE_CO2 = DEBmicroTrait.assimilation_production!(J_DE_CO2, p, D, V)
J_D = DEBmicroTrait.uptake!(J_D, p, D, V)
@test size(J_DE) == (n_reserves, n_microbes)
@test size(J_DE_CO2,1) == n_microbes
@test size(J_D,1) == n_substrates
