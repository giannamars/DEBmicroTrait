using DEBmicroTrait
using CSV, DataFrames, Statistics
using DifferentialEquations
using Gadfly
using JLD

df_metabolites = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/saifuddin_segre_2019.csv", DataFrame, missingstring="N/A")
df_pold        = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/data/IsolateMetadata190309.csv", DataFrame, missingstring="N/A")
df_out         = DataFrame()
df_traits      = DataFrame()


df_volume               = filter(:cylindrical_volume_pm3 => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df_pold)
df_out.cell_volume_pred = DEBmicroTrait.genome_size_to_cell_volume(convert(Array{Float64,1}, df_volume.actual_genome_size_bp))
df_out.cell_volume_meas = df_volume.cylindrical_volume_pm3.*1e-18

df_out.rrn_copies       = df_volume.rrN
df_out.gmax             = df_volume.MinGT
df_out.gmax_roller      = DEBmicroTrait.gmax_regression(convert(Array{Float64,1}, df_volume.rrN))

# Batch model
df_batch                = filter(:MBCug => x -> !(ismissing(x) || isnothing(x) || isnan(x)), df_pold)
genome_size             = convert(Array{Float64,1}, df_batch.actual_genome_size_bp)
rrn_copies              = convert(Array{Float64,1}, df_batch.rrN)
min_gen_time            = convert(Array{Float64,1}, df_batch.MinGT)
gram_stain              = convert(Array{String}, df_batch.Gram_stain)
α                       = df_batch.totalEEA./sum(df_batch.totalEEA).*3

n_polymers              = 0
n_monomers              = 1
n_microbes              = 1
n_enzymes               = 1
n_minerals              = 0
p_set                   = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)

BGEs                    = zeros(22,12)
rates                   = zeros(22,12)
BRs                     = zeros(22,12)

ts                      = zeros(22,12,500)
xs                      = zeros(22,12,500)
rG_CO2s                 = zeros(22,12,500)
rM_CO2s                 = zeros(22,12,500)
rX_CO2s                 = zeros(22,12,500)
J_DE_CO2s               = zeros(22,12,500)
Ds                      = zeros(22,12,500)
Es                      = zeros(22,12,500)
Vs                      = zeros(22,12,500)
Xs                      = zeros(22,12,500)
CO2s                    = zeros(22,12,500)

id_monomers             = [1,11,12]

for l = 1:22
    for m in id_monomers
        id_microbe = l
        id_monomer = m
        p                         = DEBmicroTrait.init_batch_model(id_microbe, id_monomer, df_metabolites, genome_size, rrn_copies, min_gen_time, gram_stain, α, p_set)

        u0                        = zeros(p_set.dim)
        u0[1+n_polymers:n_polymers+n_monomers] .= 1/1e-3/12.011
        u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] .= 0.5*df_batch.MBCug[id_microbe].*1e-6/12.011
        u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.5*df_batch.MBCug[id_microbe].*1e-6/12.011

        tspan                     = (0.0,df_batch.tend[id_microbe])

        prob                      = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol                       = solve(prob, alg_hints=[:stiff])

        D                         = [sol[i][1] for i in 1:length(sol.t)]
        E                         = [sol[i][2] for i in 1:length(sol.t)]
        V                         = [sol[i][3] for i in 1:length(sol.t)]
        X                         = [sol[i][4] for i in 1:length(sol.t)]
        CO2                       = [sol[i][5] for i in 1:length(sol.t)]


        du                        = zeros(p_set.dim)
        BR                        = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BRs[l,m]                  = BR[end]
        BP                        = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE                       = @. BP/(BP + BR)
        BGEs[l,m]                 = BGE[end]

        r                         = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
        rates[l,m]                = r[end]

        for i in 1:size(sol.t,1)
            x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r[i], p.metabolism_pars, [E[i]], [V[i]])
            xs[l,m,i] = x[1]
            rG_CO2s[l,m,i] = rG_CO2[1]
            rM_CO2s[l,m,i] = rM_CO2[1]
            rX_CO2s[l,m,i] = rX_CO2[1]
            J_DE_CO2s[l,m,i] = DEBmicroTrait.assimilation_production!(zeros(p.setup_pars.n_microbes), p.assimilation_pars, [D[i]], [V[i]])[1]
            ts[l,m,i] = sol.t[i]
            Ds[l,m,i] = D[i]
            Es[l,m,i] = E[i]
            Vs[l,m,i] = V[i]
            Xs[l,m,i] = X[i]
            CO2s[l,m,i] = CO2[i]
        end
    end
end

plot(x=BGEs, y=rates, Geom.point, Scale.y_log10, Scale.x_log10)

plot(x=BRs, y=BGEs, Geom.point, Scale.y_log10, Scale.x_log10)

plot(x=ts[2,1,:], y=xs[2,1,:])


JLD.save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/pold_batch_timeseries.jld", "ts", ts,
     "xs", xs, "rGCO2s", rG_CO2s, "rMCO2s", rM_CO2s, "rXCO2s", rX_CO2s, "JDECO2", J_DE_CO2s,
     "Ds", Ds, "Es", Es, "Vs", Vs, "Xs", Xs, "CO2s", CO2s)


df_traits.BGEs_out = vcat(BGEs[:,1], BGEs[:,11], BGEs[:,12])
df_traits.rates_out = vcat(rates[:,1], rates[:,11], rates[:,12])
df_traits.media_out = vcat(repeat(["Glucose"], 22), repeat(["Pyruvate"], 22), repeat(["Succinate"], 22))
CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/pold_batch_tradeoff.csv", df_traits)


CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/pold_batch_traits.csv", df_out)
