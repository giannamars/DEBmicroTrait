abstract type AbstractParams end
abstract type AbstractModel end
abstract type AbstractSetup end


struct Setup <: AbstractSetup
    n_polymers::Int64
    n_monomers::Int64
    n_microbes::Int64
    n_enzymes::Int64
    n_minerals::Int64
    dim::Int64
end

function Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)
    if n_minerals == 0
        dim = n_polymers + n_monomers + 3*n_microbes + n_enzymes
    else
        dim = n_polymers + 2*n_monomers + 3*n_microbes + 2*n_enzymes
    end
    return Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals, dim)
end



@selectable struct Params{SE,ME,AS,DEP,TU} <: AbstractParams
    # Field                     | Selectable Types
    setup_pars::SE              | AbstractSetup
    metabolism_pars::ME         | AbstractMetabolism
    assimilation_pars::AS       | AbstractAssimilation
    depolymerization_pars::DEP  | Union{nothing, AbstractDepolymerization}
    turnover_pars::TU           | AbstractTurnover
end

for fn in fieldnames(Params)
    @eval $fn(p::Params) = p.$fn
end


# @description mutable struct BatchModel{P} <: AbstractModel
#     params::P       | "Model parameters"
# end


function split_state_batch(u::AbstractVector{<:Real}, p::AbstractParams)
    n_polymers = p.setup_pars.n_polymers
    n_monomers = p.setup_pars.n_monomers
    n_microbes = p.setup_pars.n_microbes
    n_enzymes  = p.setup_pars.n_enzymes

    D    = u[1+n_polymers:n_polymers+n_monomers]
    E    = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
    V    = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
    X    = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
    CO2  = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes]
    return D, E, V, X, CO2
end


function split_state_minerals(u::AbstractVector{<:Real}, p::AbstractParams)
    n_polymers = p.setup_pars.n_polymers
    n_monomers = p.setup_pars.n_monomers
    n_microbes = p.setup_pars.n_microbes
    n_enzymes  = p.setup_pars.n_enzymes

    D    = u[1+n_polymers:n_polymers+n_monomers]
    E    = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
    V    = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
    X    = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
    D_ads = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers]
    CO2  = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers+n_microbes]

    return D, E, V, X, D_ads, CO2
end

function split_state_rhizo(u::AbstractVector{<:Real}, p::AbstractParams)
    n_polymers = p.setup_pars.n_polymers
    n_monomers = p.setup_pars.n_monomers
    n_microbes = p.setup_pars.n_microbes
    n_enzymes  = p.setup_pars.n_enzymes

    P     = u[1:n_polymers]
    D     = u[1+n_polymers:n_polymers+n_monomers]
    E     = u[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]
    V     = u[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes]
    X     = u[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes]
    D_ads = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers]
    X_ads = u[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers]
    CO2   = u[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes]
    return P, D, E, V, X, D_ads, X_ads, CO2
end
