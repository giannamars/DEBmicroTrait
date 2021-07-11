using DEBmicroTrait
using DifferentialEquations
using CSV, DataFrames, DataFramesMeta, Statistics
using Gadfly, Cairo, Fontconfig
using Dierckx
using JLD

root_density  = 0.6e-3                                                      # g/g, Rillig (2002)
bulk_density  = 1.21                                                        # g/cm^3, Thea Whitman

d = load("/Users/glmarschmann/.julia/dev/DEBplant/plant_out.jld")
ps = d["exudation_rate"]
ps_pos = abs.(ps[1:8050])
ps_neg = ps[8051:end]
ps_aux = abs.(ps[7761:end])./4
ex_rate = vcat(ps_pos, ps_neg, ps_aux).*root_density.*bulk_density
plot(y=ex_rate)

class_norm = DEBmicroTrait.exudation_class_io("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/avena_exudation.csv")
class_norm1 = reshape(repeat(class_norm[:,1], 3800), 6, 3800)
class_norm2 = reshape(repeat(class_norm[:,2], 3300), 6, 3300)
class_norm2[1,:].=0.064462
class_norm2[2,:].=0.0191312
class_norm2[3,:].=0.0199917
class_norm2[4,:].=0.609242
class_norm3 = reshape(repeat(class_norm[:,3], 964), 6, 964)
class_norm3[1,:].+=0.03333333333333333
class_norm3[2,:].+=0.03333333333333333
class_norm3[3,:].+=0.03333333333333333
class_norm3[4,:].=0.035254
class_norm4 = reshape(repeat(class_norm[:,4], 696), 6, 696)
class_norm4[1,:].+=0.03333333333333333
class_norm4[2,:].+=0.03333333333333333
class_norm4[3,:].+=0.03333333333333333
class_norm4[4,:].=0.035254
class_norm5 = reshape(repeat(class_norm[:,5], 1000), 6, 1000)
class_norm5[1,:].*=-1.0
class_norm5[2,:].*=-1.0
class_norm5[3,:].*=-1.0
class_norm5[4,:].*=-1.0
class_norm5[6,:].*=-1.0
class_norm_ex  = hcat(class_norm1, class_norm2, class_norm3, class_norm4, class_norm5)




max_exudation_rate_amino(t) = ex_rate[round(Int,t)+1].*class_norm_ex[1,round(Int,t)+1]
max_exudation_rate_organic(t) =  ex_rate[round(Int,t)+1].*class_norm_ex[2,round(Int,t)+1]
max_exudation_rate_nucleotide(t) =  ex_rate[round(Int,t)+1].*class_norm_ex[3,round(Int,t)+1]
max_exudation_rate_sugar(t) =  ex_rate[round(Int,t)+1].*class_norm_ex[4,round(Int,t)+1]
max_exudation_rate_auxin(t) =  ex_rate[round(Int,t)+1].*class_norm_ex[5,round(Int,t)+1]
max_exudation_rate_fatty(t) =  ex_rate[round(Int,t)+1].*class_norm_ex[6,round(Int,t)+1]


function max_exudation_class!(du,u,p,t)
 du[1] = 0.0
 du[2] = 0.0
 du[3] = 0.0
 du[4] = 0.0
 du[5] = 0.0
 du[6] = 0.0
end


ts = collect(range(0,12*7*24,step=1))
tmp_ts = zeros(size(ts,1))
for i in 1:size(ts,1)
    if mod(ts[i], 24) == 0
        tmp_ts[i] = ts[i]
    end
end
diurnal_ts = filter(!iszero, tmp_ts)
diurnal_ts = collect(1:8759)
diurnal_auxin = collect(8759:9759)
diurnal_t = vcat(diurnal_ts, diurnal_auxin)
cb_amino = FunctionCallingCallback((u,t,integrator)->integrator.u[1] += max_exudation_rate_amino(t), funcat= diurnal_t, func_everystep=false)
cb_organic = FunctionCallingCallback((u,t,integrator)->integrator.u[2] += max_exudation_rate_organic(t), funcat= diurnal_t, func_everystep=false)
cb_nucleo = FunctionCallingCallback((u,t,integrator)->integrator.u[3] += max_exudation_rate_nucleotide(t), funcat= diurnal_t, func_everystep=false)
cb_sugar = FunctionCallingCallback((u,t,integrator)->integrator.u[4] += max_exudation_rate_sugar(t), funcat= diurnal_t, func_everystep=false)
cb_auxin = FunctionCallingCallback((u,t,integrator)->integrator.u[5] += max_exudation_rate_auxin(t), funcat= diurnal_t, func_everystep=false)
cb_fatty = FunctionCallingCallback((u,t,integrator)->integrator.u[6] += max_exudation_rate_fatty(t), funcat= diurnal_t, func_everystep=false)

u0 = zeros(6)
tspan = (0.0,9760)
prob = ODEProblem(max_exudation_class!,u0,tspan)
sol = solve(prob, callback=CallbackSet(cb_amino,cb_organic,cb_nucleo,cb_sugar, cb_auxin,cb_fatty), tstops=diurnal_t)
amino   = [sol[i][1] for i in 1:length(sol.t)]
organic = [sol[i][2] for i in 1:length(sol.t)]
nucleo  = [sol[i][3] for i in 1:length(sol.t)]
sugar   = [sol[i][4] for i in 1:length(sol.t)]
auxin   = [sol[i][5] for i in 1:length(sol.t)]
fatty   = [sol[i][6] for i in 1:length(sol.t)]

layer1 = layer(x=sol.t, y=amino, Geom.line)
layer2 = layer(x=sol.t, y=organic, Geom.line)
layer3 = layer(x=sol.t, y=nucleo, Geom.line)
layer4 = layer(x=sol.t, y=sugar, Geom.line)
layer5 = layer(x=sol.t, y=auxin, Geom.line)
layer6 = layer(x=sol.t, y=fatty, Geom.line)

plot(layer1, layer2, layer3, layer4, layer5, layer6)
