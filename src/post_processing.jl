function calc_BGE(p, sol)
    du   = zeros(p.setup_pars.dim)
    BR   = abs.([batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)])
    BP   = [batch_model!(du, sol.u[i], p, 0)[2] + batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
    BGE  = zeros(size(sol.t,1))
    # for i in 1:size(sol.t,1)
    #     if BP[i] <= 0.0
    #         BGE[i] = NaN
    #     else
    #         BGE[i]  =  BP[i]/(BP[i] + BR[i])
    #     end
    # end
    # return BGE
    BGE = @. BP/(BP + BR)
end

function calc_growth_rate(p, sol)
    E                 = [sol[i][2] for i in 1:length(sol.t)]
    V                 = [sol[i][3] for i in 1:length(sol.t)]
    r                 = [growth!(0.0*ones(p.setup_pars.n_microbes), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
end

function calc_enzyme_production_rate(p, sol)
    E                 = [sol[i][2] for i in 1:length(sol.t)]
    V                 = [sol[i][3] for i in 1:length(sol.t)]
    r                 = [growth!(0.0*ones(p.setup_pars.n_microbes), p.metabolism_pars, [E[i]], [V[i]])[1] for i in 1:size(sol.t,1)]
    x                 = [growth_production!(r[i], p.metabolism_pars, [E[i]], [V[i]])[1][1] for i in 1:size(sol.t,1)]
end



function calc_mass_balance(p, sol)
    du     = zeros(p.setup_pars.dim)
    dD     = [batch_model!(du, sol.u[i], p, 0)[1] for i in 1:size(sol.t,1)]
    dE     = [batch_model!(du, sol.u[i], p, 0)[2] for i in 1:size(sol.t,1)]
    dV     = [batch_model!(du, sol.u[i], p, 0)[3] for i in 1:size(sol.t,1)]
    dX     = [batch_model!(du, sol.u[i], p, 0)[4] for i in 1:size(sol.t,1)]
    dCO2   = [batch_model!(du, sol.u[i], p, 0)[end] for i in 1:size(sol.t,1)]
    mass_balance = dD .+ dE .+ dV .+ dX .+ dCO2
end
