abstract type AbstractDepolymerization end

@columns struct Depolymerization{X0,VE,AKIN,KEP} <: AbstractDepolymerization
    Χ_0::X0     | _         | "Number of C atoms per product molecule"
    V_E::VE     | 1/h       | "Specific hydrolysis rate"
    α_kin::AKIN | mol/mol   | "Enzyme binding sites per polymer molecule"
    K_EP::KEP   | mol/m^3   | "Depolymerization half-saturation constant"
end

for fn in fieldnames(Depolymerization)
    @eval $fn(p::Depolymerization) = p.$fn
end

depolymerization!(J_P::Vector{Float64}, p::Depolymerization, P::Vector{Float64}, E::Vector{Float64}) = begin
    n_substrates = size(P,1)
    n_consumers  = size(E,1)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers),P, E, K_EP(p), V_E(p), α_kin(p), "depoly")
    J_P = vcat(sum(ECA[1:n_substrates, :], dims=2)...)
end
