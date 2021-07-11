using DEBmicroTrait
using DifferentialEquations, JLD
using Gadfly

d = load("/Users/glmarschmann/.julia/dev/DEBplant/plant_out.jld")
ps = d["exudation_rate"]
ps_pos = abs.(ps[1:7780])
ps_neg = ps[7781:end]
ex_rate = vcat(ps_pos, ps_neg)
plot(y=ps_pos)

rate(t) = ex_rate[t]
