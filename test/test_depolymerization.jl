using DEBmicroTrait, Test

n_polymers  = rand(1:10)
n_enzymes   = rand(1:10)

P = rand(n_polymers)
E = rand(n_enzymes)

Χ_0 = rand(n_polymers)
V_E = rand(n_polymers)
α_kin = rand(n_polymers, n_enzymes)
K_EP = rand(n_polymers, n_enzymes)

p = Depolymerization(Χ_0, V_E, α_kin, K_EP)
J_P = zeros(n_polymers)
J_P = DEBmicroTrait.depolymerization!(J_P, p, P, E)
@test size(J_P,1) == n_polymers
