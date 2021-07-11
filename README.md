# DEBmicroTrait

[![Build Status](https://travis-ci.com/giannamars/DEBmicroTrait.jl.svg?branch=master)](https://travis-ci.com/giannamars/DEBmicroTrait.jl)
[![Coverage](https://codecov.io/gh/giannamars/DEBmicroTrait.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/giannamars/DEBmicroTrait.jl)


DEBmicroTrait is a Julia package with documentation on GitHub pages. In the current version, DEBmicroTrait can be built in batch mode, chemostat mode, or with time-dependent substrate forcings for any number of substrates (polymers, monomers), microbes, and enzymes. Representations of the soil mineral matrix and microscale-environment feedback on process rates (soil temperature, soil saturation and texture) can be added based on metadata availability. 

To set up, run julia in this folder, then run:

]
instantiate
activate .

DEBmicroTrait can be coupled with DEBplant (https://github.com/rafaqz/DEBplant) to provide rhizosphere boundary conditions.
