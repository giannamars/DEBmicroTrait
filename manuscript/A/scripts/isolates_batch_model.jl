using DEBmicroTrait
using CSV, DataFrames, Statistics
using JLD
using DifferentialEquations
using Gadfly
using HypothesisTests

assimilation      = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_assimilation.jld")
enzymes           = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_enzymes.jld")
maintenance       = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_metabolism.jld")
protein_synthesis = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_protein_synthesis.jld")
turnover          = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_turnover.jld")
init              = load("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_init.jld")

condition(u,t,integrator) = u[1] - 1e-6
affect!(integrator)       = terminate!(integrator)
cb                        = ContinuousCallback(condition,affect!)


BGE_tseries       = zeros(39, 84, 500)
r_tseries         = zeros(39, 84, 500)
x_tseries         = zeros(39, 84, 500)

for i in 1:39
    for j in 1:84
        id_isolate = i
        id_monomer = j
        p                 = DEBmicroTrait.init_batch_model(id_isolate, id_monomer, assimilation, enzymes, maintenance, protein_synthesis, turnover)
        n_polymers        = p.setup_pars.n_polymers
        n_monomers        = p.setup_pars.n_monomers
        n_microbes        = p.setup_pars.n_microbes

        u0                                                                         = zeros(p.setup_pars.dim)
        # u0[1+n_polymers:n_polymers+n_monomers]                                    .= 1.0
        # u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 0.9*init["Bio0"][id_isolate]
        # u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 0.1*init["Bio0"][id_isolate]
        u0[1+n_polymers:n_polymers+n_monomers]                                    .= 25.0
        u0[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes]              .= 1e-3
        u0[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] .= 1e-3

        tspan             = (0.0,500.0)
        prob              = ODEProblem(DEBmicroTrait.batch_model!,u0,tspan,p)
        sol               = solve(prob, alg_hints=[:stiff], callback=cb)

        D    = [sol[i][1] for i in 1:length(sol.t)]
        E    = [sol[i][2] for i in 1:length(sol.t)]
        V    = [sol[i][3] for i in 1:length(sol.t)]
        X    = [sol[i][4] for i in 1:length(sol.t)]
        CO2  = [sol[i][5] for i in 1:length(sol.t)]

        du = zeros(p.setup_pars.dim)
        BR   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
        BP   = [DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[2] + DEBmicroTrait.batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
        BGE  = @. BP/(BP + BR)
        BGE_med = median(BGE[BGE.>0.0])
        BGE_all[i,j] = BGE_med

        r    = [DEBmicroTrait.growth!(0.0*ones(1), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
        r_med = median(r[r.>0.0])
        r_all[i,j] = r_med

        # BGE = DEBmicroTrait.calc_BGE(p,sol)
        # for k in 1:size(BGE,1)
        #     BGE_tseries[i,j,k] = BGE[k]
        # end
        # r_tmp = DEBmicroTrait.calc_growth_rate(p,sol)
        # for k in 1:size(r_tmp,1)
        #     r_tseries[i,j,k] = r_tmp[k]
        # end
        #
        # x_tmp = DEBmicroTrait.calc_enzyme_production_rate(p, sol)
        # for k in 1:size(x_tmp,1)
        #     x_tseries[i,j,k] = x_tmp[k]
        # end

    end
end


BGE_all = zeros(39, 84)
for i in 1:39
    for j in 1:84
        try
            #BGE_all[i,j] = median(BGE_tseries[i,j,1:499][(diff(BGE_tseries[i,j,:]) .> 0.0) .& (diff(BGE_tseries[i,j,:]) .< 1e-5)])
            #BGE_all[i,j] = BGE_tseries[i,j,1]
            BGE_med = median(BGE_tseries[i,j,:][BGE_tseries[i,j,:].>0.0])
            BGE_all[i,j] = BGE_med
        catch
            BGE_all[i,j] = NaN
        end
    end
end


BGE_all[BGE_all.>1.0].=NaN
BGE_all[BGE_all.<0.05].=NaN


r_all = zeros(39, 84)
for i in 1:39
    for j in 1:84
        try
            #r_all[i,j] = median(r_tseries[i,j,1:499][(diff(r_tseries[i,j,:]) .< 0.0) .& (diff(r_tseries[i,j,:]) .< -1e-6)])
            r_all[i,j] =  median(r_tseries[i,j,:][r_tseries[i,j,:].>0.0])
            #r_all[i,j] = r_tseries[i,j,1]

        catch
            r_all[i,j] = NaN
        end
    end
end


x_all = zeros(39, 84)
for i in 1:39
    for j in 1:84
        try
            x_all[i,j] = median(x_tseries[i,j,1:499][(diff(x_tseries[i,j,:]) .<0.0) .& (diff(x_tseries[i,j,:]) .< -1e-6)])
        catch
            x_all[i,j] = NaN
        end
    end
end



########################################
# I/O
df_isolates     = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/isolates2traits.csv", DataFrame, missingstring="N/A")
df_metabolites  = CSV.read("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/files/avena_exudation_props.csv", DataFrame, missingstring="N/A")
df_out = DataFrame()
df_out.BGE = vec(BGE_all)
df_out.rate = vec(r_all)
df_out.enz = vec(x_all)
df_out.ontology = repeat(df_metabolites.Ontology, 39)

response = Array{String}(undef, (84,39))
for i in 1:39
     response[:,i] .= df_isolates.Rhizosphere_response[i]
 end
df_out.response = vec(response)

species = Array{String}(undef, 84,39)
for i in 1:39
     species[:,i] .= df_isolates.Isolate[i]
 end
df_out.species = vec(species)

class = Array{String}(undef, 84,39)
for i in 1:39
     class[:,i] .= df_isolates.Class[i]
 end
df_out.class = vec(class)
df_out.monomer = repeat(df_metabolites.Name, 39)

CSV.write("/Users/glmarschmann/.julia/dev/DEBmicroTrait/manuscript/A/files/isolates_batch_model.csv", df_out)
########################################

########################################
# statistics
#sugars
df_p_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="negative") , df_out)
df_u_sugars = filter(x->(x.ontology.=="Sugars"&&x.response.=="undefined") , df_out)
BGE_sugars_pn  = KruskalWallisTest(df_p_sugars.BGE, df_n_sugars.BGE)
BGE_sugars_pu  = KruskalWallisTest(df_p_sugars.BGE, df_u_sugars.BGE)
rate_sugars_pn = KruskalWallisTest(df_p_sugars.rate, df_n_sugars.rate)
rate_sugars_pu = KruskalWallisTest(df_p_sugars.rate, df_u_sugars.rate)


df_p_sugars = filter(x->(x.monomer.=="salicylic acid"&&x.response.=="positive") , df_out)
df_n_sugars = filter(x->(x.monomer.=="salicylic acid"&&x.response.=="negative") , df_out)
df_u_sugars = filter(x->(x.monomer.=="salicylic acid"&&x.response.=="undefined") , df_out)
BGE_sugars_pn  = KruskalWallisTest(df_p_sugars.BGE, df_n_sugars.BGE)
BGE_sugars_pu  = KruskalWallisTest(df_p_sugars.BGE, df_u_sugars.BGE)
rate_sugars_pn = KruskalWallisTest(df_p_sugars.rate, df_n_sugars.rate)
rate_sugars_pu = KruskalWallisTest(df_p_sugars.rate, df_u_sugars.rate)
#D                 = [sol[i][1] for i in 1:length(sol.t)]
#E                 = [sol[i][2] for i in 1:length(sol.t)]
#V                 = [sol[i][3] for i in 1:length(sol.t)]
#X                 = [sol[i][4] for i in 1:length(sol.t)]
#CO2               = [sol[i][5] for i in 1:length(sol.t)]
#Bio               = E.+V
#N_cells           = [@. Bio[i]*1e-6*12.011/(init["rhoB"]*init["Md"])[1] for i in 1:length(sol.t)]
