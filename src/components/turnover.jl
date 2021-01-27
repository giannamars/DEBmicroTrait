abstract type AbstractTurnover end

@columns struct Turnover{gV0,gV1,gX,FED} <: AbstractTurnover
    γ_V_0::gV0     | 1/h       | "Max. reserve and structure turnover rate"
    γ_V_1::gV1     | mol/m^3   | "Half-saturation constant for reserve and structure turnover rate"
    γ_X::gX       | 1/h       | "Enzyme turnover rate"
    f_ED::FED      | _         | "Fraction of decayed reserve recycling to monomers"
end

for fn in fieldnames(Turnover)
    @eval $fn(p::Turnover) = p.$fn
end

biomass_turnover!(J_B::Vector{Float64}, p::AbstractTurnover, B::Vector{Float64}) = begin
    J_B = γ_V_0(p).*B./(γ_V_1(p) .+ B)
end

enzyme_decay!(J_X::Vector{Float64}, p::AbstractTurnover, X::Vector{Float64}) = begin
    J_X = γ_X(p).*X
end

reserve_recycling!(J_ED::Vector{Float64}, p::AbstractTurnover, E::Vector{Float64}) = begin
    J_ED = f_ED(p).*sum(γ_V_0(p).*E./(γ_V_1(p) .+ E), dims=1)
end

# reserve_recycling!(J_ED::Vector{Float64}, p::AbstractTurnover, E::Matrix{Float64}) = begin
#     J_ED = f_ED(p).*sum(γ_V_0(p).*E./(γ_V_1(p) .+ E), dims=1)
# end
