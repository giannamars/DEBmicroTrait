function transporter_density_to_monomer_uptake_sites(V_c::Vector{Float64}, ρ_closure::Matrix{Float64}, Min_gen_time::Vector{Float64}, Gram_stain::Vector{String})
    gmax = log(2)./Min_gen_time
    n_0 = cell_volume_to_cell_number_density(V_c, gmax, Gram_stain)
    r_p = 1e-9                    # 1 μm
    r_c = (3*V_c/(4*pi)).^(1/3)
    N_SB = zeros(size(ρ_closure,1), size(ρ_closure,2))
    for i in 1:size(V_c,1)
        N_SB[:,i] = n_0[i]*4*r_c[i]^2/(r_p^2)*ρ_closure[:,i]/12.011
    end
    N_SB[N_SB.<1e-12].=0.0
    return N_SB
end
