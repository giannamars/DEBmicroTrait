function cell_volume_to_specific_maintenance_rate(V_c::Vector{Float64}, Min_gen_time::Vector{Float64}, Gram_stain::Vector{String})
    # Lynch & Marinov (2015), Eq. 1a
    gmax = log(2)./Min_gen_time
    n_0 = cell_volume_to_cell_number_density(V_c, gmax, Gram_stain)
    V_cell = V_c./1e-18         #  V_cell in μm^3
    E_M = 0.39*V_cell.^0.88     # [1e9 ATP/(cell h)]
    E_G = 1e9*E_M/26            # 26 ATP per glucose molecule: [glucose/(cell h)]
    E_C = E_G*6*24              # [C atoms/(cell/d)]
    mol_C = E_C/(6.022e23)      # [mol C/(cell d)]
    k_M = mol_C./n_0./24        # 1/h
end

function translation_power(Protein_volume::Vector{Float64}, Ribosome_volume::Vector{Float64}, Min_gen_time::Vector{Float64})
    V_p = Protein_volume
    V_r = Ribosome_volume
    gmax = log(2)./Min_gen_time
    transl_power = @. gmax*V_p/V_r
    return max.(transl_power, 1.01*gmax)
end

function relative_translation_efficiency_regression(rrn_copies::Vector{Float64})
    log2_rna = log2.(rrn_copies)
    y = @. 9.5 - 1.22*log2_rna
    return 9.6./y
end

function gmax_regression(rrn_copies::Vector{Float64})
    log2_rna = log2.(rrn_copies)
    y = @. log2(0.06) + 1.03*log2_rna
    return 2.0.^y
end

function constrain_transporter_density(ρ_p, V_c, Min_gen_time, Gram_stain)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    y_DE = 0.22
    N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB

    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    y_EM = 1.0
    #k_M = zeros(1)

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    #
    Φ = cell_volume_to_volume_fraction(V_c, gmax, Gram_stain)
    #
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+y_EV)
    return r_m[1] - gmax[1]
end

function constrain_transporter_density(ρ_p, V_c, Min_gen_time, Gram_stain, y_DE, N_C)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    # y_DE = 0.22
    # N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB

    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    #k_M = zeros(1)
    y_EM = 1.0

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    #
    Φ = cell_volume_to_volume_fraction(V_c, gmax, Gram_stain)
    #
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+y_EV)
    return r_m[1] - gmax[1]
end

function reserve_density(ρ_p, α, V_c, Min_gen_time, Gram_stain)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    y_DE = 0.22
    N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB

    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    y_EM = 1.0

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    #
    Φ = cell_volume_to_volume_fraction(V_c, gmax, Gram_stain)
    #
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m[1]-k_M*y_EM)/(j_EA_m[1]/k_E+y_EV)
    m_E = @. (k_M*y_EM + (1+α)*r_m*y_EV)/(k_E-r_m)
    return m_E[1]
end

function reserve_density_curve(ρ_p, α, V_c, Min_gen_time, Gram_stain)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    y_DE = 0.22
    N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB

    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    y_EM = 1.0

    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    #
    Φ = cell_volume_to_volume_fraction(V_c, gmax, Gram_stain)
    #
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)

    r_m = @. (j_EA_m[1]-k_M*y_EM)/(j_EA_m[1]/k_E+y_EV)
    rtmp = convert(Array{Float64,1}, LinRange(1e-5, 1e1, 100000))
    filter!(x->x<=gmax[1],rtmp)
    r = rtmp
    m_E = @. (k_M*y_EM + (1+α)*r*y_EV)/(k_E-r)
    return m_E
end

function rate_yield_trade_off(ρ_p, α, V_c, Min_gen_time, Gram_stain)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    y_DE = 0.22
    N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB
    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    #k_M = zeros(1)
    y_EM = 1.0
    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)
    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+(1+α)*y_EV)

    #rtmp = [1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
    rtmp = convert(Array{Float64,1}, LinRange(1e-5, 1e1, 100000))
    filter!(x->x<=gmax[1],rtmp)
    #filter!(x->x>=k_M[1],rtmp)
    #r = vcat(k_M[1], rtmp, r_m[1])
    r = rtmp
    m_E = @. (k_M*y_EM + (1+α)*r*y_EV)/(k_E-r)
    j_EA = @. m_E*k_E
    yield_coeff = r./j_EA
    return yield_coeff, convert(Array{Float64,1}, r)
end

function rate_yield_trade_off(ρ_p, α, V_c, Min_gen_time, Gram_stain, y_DE, N_C)
    N_SB = transporter_density_to_monomer_uptake_sites(V_c, ρ_p*ones(1,1), Min_gen_time, Gram_stain)
    k2p = 180.0*60^2
    #y_DE = 0.22
    #N_C = 6
    j_EA_m = (1.0/y_DE - 1.0)*N_C*k2p*N_SB
    k_M = cell_volume_to_specific_maintenance_rate(V_c, Min_gen_time, Gram_stain)
    #k_M = zeros(1)
    y_EM = 1.0
    gmax   = log(2)./Min_gen_time
    V_p = cell_volume_to_protein_volume(V_c)
    V_r = cell_volume_to_ribosome_volume(V_c, gmax)
    k_E = translation_power(V_p, V_r, Min_gen_time)
    L_DNA = cell_volume_to_genome_size(V_c)
    rrn_copies = genome_size_to_rRNA_copy_number(L_DNA)
    y_EV = relative_translation_efficiency_regression(rrn_copies)
    r_m = @. (j_EA_m-k_M*y_EM)/(j_EA_m/k_E+(1+α)*y_EV)

    #rtmp = [1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2]
    rtmp = convert(Array{Float64,1}, LinRange(1e-5, 1e1, 100000))
    filter!(x->x<=gmax[1],rtmp)
    #filter!(x->x>=k_M[1],rtmp)
    #r = vcat(k_M[1], rtmp, r_m[1])
    r = rtmp
    m_E = @. (k_M*y_EM + (1+α)*r*y_EV)/(k_E-r)
    j_EA = @. y_DE.*m_E*k_E
    yield_coeff = r./j_EA
    return yield_coeff, convert(Array{Float64,1}, r)
end


function g_Φ_berg(Φ)
    R_w = 0.138
    R_e = 3.06
    R_i = 5.57
    Φ_w_0 = 0.363
    a = R_e/R_w
    b = R_i/R_e
    f_a = (3*(1-Φ_w_0)/(a*(1+2*Φ_w_0))^2)*(a*(1+2*Φ_w_0) + (1-Φ_w_0)*(1+(1-Φ_w_0)/(3*a*Φ_w_0)))
    Φ_w = Φ_w_0*(1-Φ-Φ*f_a)
    S_3 = Φ_w + Φ
    S_2 = 0.5*(a*Φ_w+Φ)/R_e
    S_1 = 0.25*(a^2*Φ_w+Φ)/R_e^2
    g_Φ = (-log(1-S_3) + 6*S_2*R_i/(1-S_3) + (12*S_1/(1-S_3) + 18*S_2^2/(1-S_3)^2)*R_i^2)
    return g_Φ
end

function activity_coefficient_berg(Φ)
    ln_gamma = g_Φ_berg(Φ) - g_Φ_berg(0.0)
    return exp(ln_gamma)
end


function constrain_enzyme_production(α, hydrolase_distr)
  genome_distr_norm = hydrolase_distr./sum(hydrolase_distr)
  closure = α.*genome_distr_norm
end
