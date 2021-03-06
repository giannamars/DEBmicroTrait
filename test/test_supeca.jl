using DEBmicroTrait, Test

# Assimilation
n_substrates = rand(1:10)
n_consumers  = rand(1:10)

F_c          = zeros(n_consumers)
F_r          = zeros(n_substrates)
ECA          = zeros(n_substrates,n_consumers)
S            = rand(n_substrates)
E            = rand(n_consumers)
K            = rand(n_substrates,n_consumers)
N_SB         = rand(n_substrates,n_consumers)
k2p          = rand(n_substrates)

F_c = DEBmicroTrait.normalized_substrate_flux!(F_c, S, K)
@test size(F_c) == (n_consumers,)
F_r = DEBmicroTrait.conjugate_substrate_flux!(F_r, E, K, N_SB)
@test size(F_r) == (n_substrates,)
ECA = DEBmicroTrait.ECA_kinetics!(ECA, S, E, K, k2p, N_SB)
@test size(ECA) == (n_substrates,n_consumers)

# Depolymerization
n_substrates  = rand(1:10)
n_consumers   = rand(1:10)

P = rand(n_substrates)
E = rand(n_consumers)

Χ_0 = rand(n_substrates)
V_E = rand(n_substrates)
α_kin = rand(n_substrates, n_consumers)
K_EP = rand(n_substrates, n_consumers)

F_c = zeros(n_consumers)
F_c = DEBmicroTrait.normalized_substrate_flux!(F_c, P, K_EP, α_kin)
@test size(F_c) == (n_consumers,)
F_r = zeros(n_substrates)
F_r = DEBmicroTrait.conjugate_substrate_flux!(F_r, E, K_EP)
@test size(F_r) == (n_substrates,)
ECA = zeros(n_substrates, n_consumers)
ECA = DEBmicroTrait.ECA_kinetics!(ECA, P, E, K_EP, V_E, α_kin, "depoly")
@test size(ECA) == (n_substrates, n_consumers)
