using DEBmicroTrait, Test

n_microbes = rand(1:10)
n_enzymes  = rand(1:10)
n_substrates = rand(1:10)

γ_V_0 = ones(n_microbes)
γ_V_1 = ones(n_microbes)
γ_X   = ones(n_enzymes)
f_ED  = ones(n_substrates)

p     = Turnover(γ_V_0,γ_V_1,γ_X,f_ED)

B     = rand(n_microbes)
J_B   = zeros(n_microbes)
J_B   = DEBmicroTrait.biomass_turnover!(J_B, p, B)
@test size(J_B,1) == n_microbes

X     = rand(n_enzymes)
J_X   = zeros(n_enzymes)
J_X   = DEBmicroTrait.enzyme_decay!(J_X, p, X)
@test size(J_X,1) == n_enzymes

E     = rand(n_microbes)
J_ED  = rand(n_microbes)
J_ED  = DEBmicroTrait.reserve_recycling!(J_ED, p, E)
@test size(J_ED,1) == n_substrates
