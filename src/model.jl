# (model::BatchModel)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Real) = begin
#     D, E, V, X, CO2 = split_state(du, p)
# end



function init_batch_model(id_microbe, id_monomer, df_metabolites, Genome_size, p_set)
    V_cs = genome_size_to_cell_volume([Genome_size[id_microbe]])
    rrn_copies = genome_size_to_rRNA_copy_number([Genome_size[id_microbe]])
    gmax = gmax_regression(rrn_copies)
    Min_gen_time = log(2)./gmax
    Gram_stain = ["+"]
    α = zeros(1)

    V_p = cell_volume_to_protein_volume(V_cs)
    V_r = cell_volume_to_ribosome_volume(V_cs, gmax)
    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)
    k_M = cell_volume_to_specific_maintenance_rate(V_cs, Min_gen_time, Gram_stain)
    y_EM = ones(1)
    α = zeros(1)
    y_EX = ones(1)

    p_met = MetabolismC(k_E,y_EV,k_M,y_EM,α,y_EX)

    y_DE = df_metabolites.y_DE[id_monomer]
    N_C = df_metabolites.N_C[id_monomer]
    D_S = [df_metabolites.diffusivity[id_monomer]]

    find_ρ(x) = constrain_transporter_density(x, V_cs, Min_gen_time, Gram_stain, y_DE, N_C)
    ρ_p = Roots.find_zero(find_ρ, 1.0)

    N_SB = transporter_density_to_monomer_uptake_sites(V_cs, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    K_D = specific_reference_affinity(V_cs, ρ_p*ones(1,1), D_S)

    p_ass = AssimilationC(N_SB,K_D,[y_DE],[N_C])

    γ_V0 = 4e-5*ones(1)
    γ_V1 = 1e-2*ones(1)
    γ_X = 2.5e-4*ones(1)
    f_ED = ones(1)

    p_turn = Turnover(γ_V0,γ_V1,γ_X,f_ED)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end

function init_batch_model(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set)
    V_cs = genome_size_to_cell_volume([Genome_size[id_microbe]])
    rrn_copies = [rrn_copies[id_microbe]]
    Min_gen_time = [Min_gen_time[id_microbe]]
    gmax = log(2)./Min_gen_time
    Gram_stain = [Gram_stain[id_microbe]]
    α = [α[id_microbe]]

    V_p = cell_volume_to_protein_volume(V_cs)
    V_r = cell_volume_to_ribosome_volume(V_cs, gmax)
    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)
    k_M = cell_volume_to_specific_maintenance_rate(V_cs, Min_gen_time, Gram_stain)
    y_EM = ones(1)
    y_EX = ones(1)

    p_met = MetabolismC(k_E,y_EV,k_M,y_EM,α,y_EX)

    y_DE = df_metabolites.y_DE[id_monomer]
    N_C = df_metabolites.N_C[id_monomer]
    D_S = [df_metabolites.diffusivity[id_monomer]]

    find_ρ(x) = constrain_transporter_density(x, V_cs, Min_gen_time, Gram_stain, y_DE, N_C)
    ρ_p = Roots.find_zero(find_ρ, 1.0)

    N_SB = transporter_density_to_monomer_uptake_sites(V_cs, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    K_D = specific_reference_affinity(V_cs, ρ_p*ones(1,1), D_S)

    p_ass = AssimilationC(N_SB,K_D,[y_DE],[N_C])

    γ_V0 = 4e-5*ones(1)
    γ_V1 = 1e-2*ones(1)
    γ_X = 2.5e-4*ones(1)
    f_ED = ones(1)

    p_turn = Turnover(γ_V0,γ_V1,γ_X,f_ED)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end




function batch_model!(du, u, p, t)
    D, E, V, X, CO2 = DEBmicroTrait.split_state(u, p)

    # metabolism
    r = growth!(0.0*ones(1), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)

    # assimilation
    J_DE = assimilation!(zeros(p.setup_pars.n_microbes), p.assimilation_pars, D, V)
    J_DE_CO2 = assimilation_production!(zeros(p.setup_pars.n_microbes), p.assimilation_pars, D, V)
    J_D = uptake!(zeros(p.setup_pars.n_monomers), p.assimilation_pars, D, V)

    # turnover
    J_ED = reserve_recycling!(zeros(p.setup_pars.n_microbes), p.turnover_pars, E)
    J_X  = enzyme_decay!(zeros(p.setup_pars.n_enzymes), p.turnover_pars, X)
    J_V  = biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, V)
    J_E  = biomass_turnover!(zeros(p.setup_pars.n_microbes), p.turnover_pars, E)

    @. du[1+p.setup_pars.n_polymers:p.setup_pars.n_polymers+p.setup_pars.n_monomers] = -J_D + J_ED +  J_X
    @. du[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers:p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes] = J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes] = r*V - J_V
    @. du[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes] = x*V - J_X
    @. du[1+p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes:p.setup_pars.n_polymers+p.setup_pars.n_monomers+2*p.setup_pars.n_microbes+p.setup_pars.n_enzymes+p.setup_pars.n_microbes]  = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end
