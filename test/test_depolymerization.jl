using DEBmicroTrait, Test

n_polymers = 5
n_enzymes  = 6
P = rand(n_polymers)
E = rand(n_enzymes)

Χ_0 = [10,10,10,10,10]
ρ_P = [1.35, 1.5, 1.4, 0.4, 1.0]
R_E = [3.29e-9, 4e-9, 4e-9, 4.4e-9, 5e-9,3.9e-9]
V_E = [5, 0.24, 4,0.19, 3,4]
MW_E = [48e3, 65e3, 59.7e3, 70e3, 50e3, 52e3]
D_E = DEBmicroTrait.aqueous_diffusivity(MW_E)
α_kin_P = DEBmicroTrait.binding_sites(ρ_P, R_E)
K_EP_P =  zeros(n_polymers,n_enzymes)
tmp = DEBmicroTrait.cellulase_affinity(R_E, V_E, D_E)
for i in 1:5
    K_EP_P[i,:] = tmp
end
K_EP_P

p = Depolymerization(Χ_0, V_E, α_kin_P, K_EP_P)
J_P = zeros(n_polymers)
J_P = DEBmicroTrait.depolymerization!(J_P, p, P, E)
@test size(J_P,1) == n_polymers

n_minerals   = 1

M            = rand(n_minerals)
α_kin_M      = ones(n_minerals, n_enzymes)
K_MX         = rand(n_minerals, n_enzymes)
α_kin        = vcat(α_kin_P, α_kin_M)
K_EP         = vcat(K_EP_P,K_MX)

p_depoly     = DepolymerizationM(Χ_0, V_E, α_kin, K_EP, M)
J_P, J_XM    = DEBmicroTrait.depolymerization!(zeros(n_polymers), p_depoly, P, E)
