abstract type AbstractParams end
abstract type AbstractModel end


@selectable struct Params{ME,AS,DEP,TU} <: AbstractParams
    # Field                     | Selectable Types
    metabolism_pars::ME         | AbstractMetabolism
    assimilation_pars::AS       | AbstractAssimilation
    depolymerization_pars::DEP  | AbstractDepolymerization
    turnover_pars::TU           | AbstractTurnover
end

for fn in fieldnames(Params)
    @eval $fn(p::Params) = p.$fn
end


@description mutable struct BatchModel{P} <: AbstractModel
    params::P       | "Model parameters"
end
