module DEBmicroTrait

using FieldMetadata,
      Unitful,
      Roots,
      LinearAlgebra

import FieldMetadata: @default, @description, @units, @bounds, @logscaled, @flattenable, @plottable, @selectable,
                      default, description, units, bounds, logscaled, flattenable, plottable, selectable

# Field metadata columns
 @chain columns @units @description

using Unitful: J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol
Unitful.register(@__MODULE__);
@unit bp "bp" Basepair 1u"m" true;

export AbstractMetabolism, AbstractMetabolismC, MetabolismC, AbstractMetabolismCN, MetabolismCN, AbstractMetabolismCNP, MetabolismCNP
export AbstractAssimilation, AbstractAssimilationC, AssimilationC, AbstractAssimilationCN, AssimilationCN, AbstractAssimilationCNP, AssimilationCNP
export AbstractDepolymerization, Depolymerization
export AbstractTurnover, Turnover
export AbstractParams, Params
export AbstractIsolate, IsolateComposition
export AbstractSetup, Setup

include("components/metabolism.jl")
include("components/supeca.jl")
include("components/assimilation.jl")
include("components/depolymerization.jl")
include("components/synthesizing_units.jl")
include("components/turnover.jl")
include("components/thermostoichwizard.jl")
include("allometry/cell_composition.jl")
include("allometry/cell_assimilation.jl")
include("allometry/cell_metabolism.jl")
include("setup.jl")
include("model.jl")


end
