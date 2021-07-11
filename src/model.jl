########################################
function init_batch_model(id_isolate, id_monomer, assimilation, enzymes, maintenance, protein_synthesis, turnover)

    n_polymers = 0
    n_monomers = 1
    n_microbes = 1
    n_enzymes  = 1
    n_minerals = 0
    p_set      = Setup(n_polymers, n_monomers, n_microbes, n_enzymes, n_minerals)

    k_E        = protein_synthesis["kE"][id_isolate]
    y_EV       = protein_synthesis["yEV"][id_isolate]
    k_M        = maintenance["kM"][id_isolate]
    y_EM       = assimilation["yEM"][id_isolate]
    α_X        = enzymes["alpha"][id_isolate]
    y_EX       = y_EV
    f_αX       = ones(n_enzymes)
    min_gt     = protein_synthesis["mingt"][id_isolate]

    p_met      = MetabolismC([k_E], [y_EV], [k_M], [y_EM], [α_X], [y_EX], f_αX, min_gt)

    N_SB       = assimilation["NSB"][id_monomer,id_isolate]
    K_D        = assimilation["KD"][id_monomer,id_isolate]
    y_DE       = assimilation["yDE"][id_monomer]
    N_C        = assimilation["NC"][id_monomer]

    p_ass      = AssimilationC(N_SB*ones(1,1),K_D*ones(1,1),[y_DE],[N_C])

    #γ_V0       = turnover["gV0"][id_isolate]
    γ_V0        = 1e-5
    γ_V1        = 1e-2
    #γ_V1       = turnover["gV1"]
    γ_X        = 1/(7*24)*ones(n_enzymes)

    γ_D_ads    = zeros(n_monomers)
    γ_X_ads    = zeros(n_enzymes)
    f_ED       = ones(n_monomers)
    f_V        = ones(n_microbes)    # all structural biomass recycling to monomers
    f_VD       = ones(n_microbes)
    f_VP       = ones(n_microbes)
    f_X        = zeros(n_enzymes)    # all enzymes recycling to monomers
    f_XD       = ones(n_enzymes)
    f_XP       = ones(n_enzymes)

    p_turn     = Turnover([γ_V0],[γ_V1],γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end
########################################

########################################
function init_batch_model_old(id_microbe, id_monomer, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)
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
    y_EX = relative_translation_efficiency_regression(rrn_copies)


    f_αX       = ones(1)
    p_met = MetabolismC(k_E,y_EV,k_M,y_EM,α,y_EX,f_αX,Min_gen_time)

    y_DE = df_metabolites.y_DE[id_monomer]
    N_C = df_metabolites.N_C[id_monomer]
    D_S = [df_metabolites.diffusivity[id_monomer]]

    K_D = K_D[id_monomer,id_microbe]*ones(1,1)
    N_SB = N_SB[id_monomer,id_microbe]*ones(1,1)
    p_ass = AssimilationC(N_SB,K_D,y_DE,N_C)

    γ_V0 = 2e-5*ones(1)
    γ_V1 = 1e-2*ones(1)
    γ_X = 2.5e-4*ones(1)
    γ_D_ads = 1e-6*ones(1)

    f_ED       = ones(1)
    f_V        = ones(1)    # all structural biomass recycling to monomers
    f_VD       = ones(1)
    f_VP       = ones(1)
    f_X        = zeros(1)    # all enzymes recycling to monomers
    f_XD       = ones(1)
    f_XP       = ones(1)

    p_turn = Turnover(γ_V0,γ_V1,γ_X,γ_D_ads,γ_X_ads,f_ED,f_VD,f_VP,f_V,f_XD,f_XP,f_X)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end

########################################


########################################
# batch model
function batch_model!(du, u, p, t)
    D, E, V, X, CO2 = DEBmicroTrait.split_state_batch(u, p)

    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes

    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)

    # assimilation
    J_DE         = assimilation!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_D          = uptake!(zeros(n_monomers), p.assimilation_pars, D, V)

    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E)
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V)
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)

    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_ED + J_VD + J_XD
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes] = rG_CO2 + rX_CO2 + rM_CO2 + J_DE_CO2
    return du
end
########################################

########################################
function init_batch_model_minerals(id_microbe, id_monomer, id_mineral, df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)
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
    y_EX = relative_translation_efficiency_regression(rrn_copies)

    p_met = MetabolismC(k_E,y_EV,k_M,y_EM,α,y_EX)

    y_DE = df_metabolites.y_DE[id_monomer]
    N_C = df_metabolites.N_C[id_monomer]
    D_S = [df_metabolites.diffusivity[id_monomer]]

    K_Dass = hcat(K_D[id_monomer,id_microbe]*ones(1,1), K_D[id_monomer,end]*ones(1,1))
    N_SBass = hcat(N_SB[id_monomer,id_microbe]*ones(1,1), N_SB[id_monomer,end]*ones(1,1))
    M_glucose = [35,173,126].*1e-6
    M_alanine = [98,527,444].*1e-6
    M_salicyclic = [240,502,481].*1e-6
    M_sinpyl = [2031, 214, 521].*1e-6
    M_oxalic = [886, 471, 955].*1e-6
    M = hcat(M_alanine, M_glucose, M_salicyclic, M_sinpyl, M_oxalic)
    p_ass = AssimilationCM(N_SBass,K_Dass,y_DE,N_C,[M[id_mineral,id_monomer]])

    γ_V0 = 2e-5*ones(1)
    γ_V1 = 1e-2*ones(1)
    γ_X = 2.5e-4*ones(1)
    γ_D_ads = 1e-6*ones(1)
    f_ED = ones(1)
    f_XD = ones(size(y_DE,1))./size(y_DE,1)
    p_turn = Turnover(γ_V0,γ_V1,γ_X,γ_D_ads,f_ED,f_XD)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end



function batch_model_minerals!(du, u, p, t)
    D, E, V, X, D_ads, CO2 = split_state_minerals(u, p)

    # metabolism
    r = growth!(0.0*ones(1), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)

    # assimilation
    J_DE = assimilation!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_DE_CO2 = assimilation_production!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_D, J_M = uptake!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)

    # turnover
    J_ED = reserve_recycling!(zeros(n_microbes), p.turnover_pars, E)
    J_X  = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD = enzyme_recycling!(zeros(n_enzymes), p.turnover_pars, X)
    J_V  = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_E  = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)
    J_D_ads = adsorbed_monomer_decay!(zeros(n_monomers), p.turnover_pars, D_ads)

    # model
    @. du[1+n_polymers:n_polymers+n_monomers] = -J_D - J_M + J_ED +  J_XD + J_D_ads
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] = J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = x*V - J_X
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_M - J_D_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers+n_microbes]  = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end



################################################################################
function init_hydroponics(df_metabolites, Genome_size, rrn_copies, Min_gen_time, Gram_stain, α, p_set, N_SB, K_D)
    V_cs = genome_size_to_cell_volume(Genome_size)
    gmax = log(2)./Min_gen_time

    V_p = cell_volume_to_protein_volume(V_cs)
    V_r = cell_volume_to_ribosome_volume(V_cs, gmax)
    k_E = translation_power(V_p, V_r, Min_gen_time)
    y_EV = relative_translation_efficiency_regression(rrn_copies)
    k_M = cell_volume_to_specific_maintenance_rate(V_cs, Min_gen_time, Gram_stain)
    y_EM = ones(39)
    y_EX = relative_translation_efficiency_regression(rrn_copies)
    p_met = MetabolismC(k_E,y_EV,k_M,y_EM,α,y_EX)

    y_DE = df_metabolites.y_DE
    N_C = df_metabolites.N_C
    D_S = df_metabolites.diffusivity
    p_ass = AssimilationC(N_SB,K_D,y_DE,N_C)

    γ_V0 = 1e-6*ones(39)
    γ_V1 = 1e-1*ones(39)
    γ_X = 1/(24*7*3600)*ones(1)
    f_ED = ones(size(y_DE,1))./size(y_DE,1)
    f_XD = ones(size(y_DE,1))./size(y_DE,1)
    p_turn = Turnover(γ_V0,γ_V1,γ_X,f_ED,f_XD)

    p = Params(p_set,p_met,p_ass,nothing,p_turn)
end


function hydroponics!(du, u, p, t)
    D, E, V, X, CO2 = DEBmicroTrait.split_state(u, p)

    # metabolism
    r = growth!(0.0*ones(39), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)

    # assimilation
    J_DE = assimilation!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_DE_CO2 = assimilation_production!(zeros(n_microbes), p.assimilation_pars, D, V)
    J_D = uptake!(zeros(n_monomers), p.assimilation_pars, D, V)
    # turnover
    J_ED = reserve_recycling!(zeros(n_microbes), p.turnover_pars, E)
    J_X  = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD = enzyme_recycling!(zeros(n_enzymes), p.turnover_pars, X)
    J_V  = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_E  = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)

    # enzyme production
    J_EX = sum(x.*V, dims=1)

    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_ED +  J_XD
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] = J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_microbes]  = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end


function rhizo_model!(du, u, p, t)

    # setup
    P, D, E, V, X, D_ads, X_ads, CO2 = split_state_rhizo(u, p)
    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes
    n_minerals                = p.setup_pars.n_minerals

    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)

    # assimilation
    J_DE         = assimilation!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_D, J_DM    = uptake!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)

    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E)
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V)
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V)
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E)
    J_D_ads      = adsorbed_monomer_decay!(zeros(n_monomers), p.turnover_pars, D_ads)
    J_X_ads      = adsorbed_enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X_ads)

    # depolymerization
    J_P, J_XM    = depolymerization!(zeros(n_polymers), p.depolymerization_pars, P, X)

    # rhs
    @. du[1:n_polymers] = - J_P + J_VP + J_XP
    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_P - J_DM + J_ED + J_VD + J_XD + J_D_ads
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X - J_XM + J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_DM - J_D_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers] = J_XM - J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end


function rhizo_model_linear!(du, u, p, t)

    # setup
    P, D, E, V, X, D_ads, X_ads, CO2 = split_state_rhizo(u, p)
    n_polymers                = p.setup_pars.n_polymers
    n_monomers                = p.setup_pars.n_monomers
    n_microbes                = p.setup_pars.n_microbes
    n_enzymes                 = p.setup_pars.n_enzymes
    n_minerals                = p.setup_pars.n_minerals

    # metabolism
    r                         = growth!(0.0*ones(n_microbes), p.metabolism_pars, E, V)
    x, rG_CO2, rM_CO2, rX_CO2 = growth_production!(r, p.metabolism_pars, E, V)
    J_EX                      = enzyme_production!(x, p.metabolism_pars, V)

    # assimilation
    J_DE         = assimilation!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_DE_CO2     = assimilation_production!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)
    J_D, J_DM    = uptake!(zeros(n_monomers,n_microbes+n_minerals), p.assimilation_pars, D, V)

    # turnover
    J_ED         = reserve_recycling!(zeros(n_monomers), p.turnover_pars, E, "lin")
    J_X          = enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X)
    J_XD, J_XP   = enzyme_recycling!(zeros(n_monomers), p.turnover_pars, X)
    J_V          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, V, "lin")
    J_VD, J_VP   = biomass_recycling!(zeros(n_monomers), p.turnover_pars, V, "lin")
    J_E          = biomass_turnover!(zeros(n_microbes), p.turnover_pars, E, "lin")
    J_D_ads      = adsorbed_monomer_decay!(zeros(n_monomers), p.turnover_pars, D_ads)
    J_X_ads      = adsorbed_enzyme_decay!(zeros(n_enzymes), p.turnover_pars, X_ads)

    # depolymerization
    J_P, J_XM    = depolymerization!(zeros(n_polymers), p.depolymerization_pars, P, X)

    # rhs
    @. du[1:n_polymers] = - J_P + J_VP + J_XP
    @. du[1+n_polymers:n_polymers+n_monomers] = - J_D + J_P - J_DM + J_ED + J_VD + J_XD + J_D_ads
    @. du[1+n_polymers+n_monomers:n_polymers+n_monomers+n_microbes] =  J_DE - (p.metabolism_pars.k_E - r)*E - J_E
    @. du[1+n_polymers+n_monomers+n_microbes:n_polymers+n_monomers+2*n_microbes] = r*V - J_V
    @. du[1+n_polymers+n_monomers+2*n_microbes:n_polymers+n_monomers+2*n_microbes+n_enzymes] = J_EX - J_X - J_XM + J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes:n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers] = J_DM - J_D_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers] = J_XM - J_X_ads
    @. du[1+n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers:n_polymers+n_monomers+2*n_microbes+2*n_enzymes+n_monomers+n_microbes] = rG_CO2 + rM_CO2 + rX_CO2 + J_DE_CO2

    return du
end
