abstract type AbstractAssimilation end

abstract type AbstractAssimilationC <: AbstractAssimilation end
abstract type AbstractAssimilationCN <: AbstractAssimilation end
abstract type AbstractAssimilationCNP <: AbstractAssimilation end


@columns struct AssimilationC{NSB,KD,YDE,NC} <: AbstractAssimilationC
    N_SB::NSB   | mol/mol   | "Mol transporters per mol of biomass carbon"
    K_D::KD     | mol/m^3   | "Reference half-saturation constant"
    y_DE::YDE   | mol/mol   | "Yield on assimilation"
    N_C::NC     | _         | "Number of C atoms per monomer"
end

for fn in fieldnames(AssimilationC)
    @eval $fn(p::AssimilationC) = p.$fn
end

@columns struct AssimilationCM{NSB,KD,YDE,NC,M} <: AbstractAssimilationC
    N_SB::NSB   | mol/mol   | "Mol transporters per mol of biomass carbon"
    K_D::KD     | mol/m^3   | "Reference half-saturation constant"
    y_DE::YDE   | mol/mol   | "Yield on assimilation"
    N_C::NC     | _         | "Number of C atoms per monomer"
    M::M        | _         | "Mineral surfaces"
end

for fn in fieldnames(AssimilationCM)
    @eval $fn(p::AssimilationCM) = p.$fn
end

assimilation!(J_DE::Vector{Float64}, p::AssimilationC, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    #J_DE = vcat(sum((1.0./y_DE(p) .- 1.0).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE = vcat(sum((1.0 .- y_DE(p)).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end

assimilation!(J_DE::Matrix{Float64}, p::AssimilationCM, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    n_minerals   = size(M(p),1)
    k2p_D        = 180.0*60*60*ones(n_substrates, n_consumers)
    k2p_M        = ones(n_substrates, n_minerals)
    k2p          = hcat(k2p_D, k2p_M)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers+n_minerals), D, vcat(V, M(p)), K_D(p), k2p, N_SB(p))
    #J_DE = vcat(sum((1.0./y_DE(p) .- 1.0).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE = vcat(sum((1.0 .- y_DE(p)).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end


@columns struct AssimilationCN{NSB,KD,YDE,NC,NN} <: AbstractAssimilationCN
    N_SB::NSB   | mol/mol   | "Mol transporters per mol of biomass carbon"
    K_D::KD     | mol/m^3   | "Reference half-saturation constant"
    y_DE::YDE   | mol/mol   | "Yield on assimilation"
    N_C::NC     | _         | "Number of C atoms per monomer"
    N_N::NN     | _         | "Number of N atoms per monomer"
end

for fn in fieldnames(AssimilationCN)
    @eval $fn(p::AssimilationCN) = p.$fn
end

assimilation!(J_DE::Matrix{Float64}, p::AssimilationCN, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_DE[1,:] = vcat(sum((1.0./y_DE(p) .- 1.0).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE[2,:] = vcat(sum((1.0./y_DE(p) .- 1.0).*N_N(p).*ECA[:, 1:n_consumers], dims=1)...)
    return J_DE
end

@columns struct AssimilationCNP{NSB,KD,YDE,NC,NN,NP} <: AbstractAssimilationCNP
    N_SB::NSB   | mol/mol   | "Mol transporters per mol of biomass carbon"
    K_D::KD     | mol/m^3   | "Reference half-saturation constant"
    y_DE::YDE   | mol/mol   | "Yield on assimilation"
    N_C::NC     | _         | "Number of C atoms per monomer"
    N_N::NN     | _         | "Number of N atoms per monomer"
    N_P::NP     | _         | "Number of P atoms per monomer"
end

for fn in fieldnames(AssimilationCNP)
    @eval $fn(p::AssimilationCNP) = p.$fn
end

assimilation!(J_DE::Matrix{Float64}, p::AssimilationCNP, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_DE[1,:] = vcat(sum((1.0./y_DE(p) .- 1.0).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE[2,:] = vcat(sum((1.0./y_DE(p) .- 1.0).*N_N(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE[3,:] = vcat(sum((1.0./y_DE(p) .- 1.0).*N_P(p).*ECA[:, 1:n_consumers], dims=1)...)
    return J_DE
end

uptake!(J_D::Vector{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    J_D = vcat(sum(ECA[:, 1:n_consumers], dims=2)...).*N_C(p)
end

uptake!(J_D::Matrix{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    n_minerals   = size(M(p),1)
    k2p_D        = 180.0*60*60*ones(n_substrates, n_consumers)
    k2p_M        = ones(n_substrates, n_minerals)
    k2p          = hcat(k2p_D, k2p_M)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers+n_minerals), D, vcat(V, M(p)), K_D(p), k2p, N_SB(p))
    J_D  = vcat(sum(ECA[:, 1:n_consumers], dims=2)...).*N_C(p)
    J_DM  = vcat(sum(ECA[:, 1+n_consumers:n_consumers+n_minerals], dims=2)...)
    return J_D, J_DM
end

assimilation_production!(J_DE_CO2::Vector{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    k2p = 180.0*60*60*ones(n_substrates)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers), D, V, K_D(p), k2p, N_SB(p))
    #J_DE_CO2 = vcat(sum((1.0./y_DE(p)).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE_CO2 = vcat(sum(y_DE(p).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end

assimilation_production!(J_DE_CO2::Matrix{Float64}, p::AbstractAssimilation, D::Vector{Float64}, V::Vector{Float64}) = begin
    n_substrates = size(D,1)
    n_consumers  = size(V,1)
    n_minerals   = size(M(p),1)
    k2p_D        = 180.0*60*60*ones(n_substrates, n_consumers)
    k2p_M        = ones(n_substrates, n_minerals)
    k2p          = hcat(k2p_D, k2p_M)
    ECA  = ECA_kinetics!(zeros(n_substrates,n_consumers+n_minerals), D, vcat(V, M(p)), K_D(p), k2p, N_SB(p))
    #J_DE_CO2 = vcat(sum((1.0./y_DE(p)).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
    J_DE_CO2 = vcat(sum(y_DE(p).*N_C(p).*ECA[:, 1:n_consumers], dims=1)...)
end
