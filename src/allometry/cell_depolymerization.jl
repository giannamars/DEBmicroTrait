function binding_sites(ρ_P::Vector{Float64}, R_E::Vector{Float64})
    σ_EP = 0.1
    #ρ_P  = 1.35     # g/cm^3
    #R_E  = 3.29e-7  # cm
    R_P  = 20e-3
    N_A  = 6.022e23
    α_kin = zeros(size(ρ_P,1), size(R_E,1))
    for i in 1:size(ρ_P,1)
        for j in 1:size(R_E,1)
            α_kin[i,j] = @. 3*σ_EP/(pi*ρ_P[i]*R_P*R_E[j]^2*N_A)
        end
    end
    return α_kin
end

function cellulase_affinity(R_E::Vector{Float64}, V_E::Vector{Float64}, D_E::Vector{Float64})
    σ_EP = 0.1
    R_P  = 20e-3
    N_A  = 6.022e23
    p_inv = @. (4*σ_EP*R_P + pi*R_E)/(4*σ_EP*R_P)
    K_EP_0 = @. σ_EP*R_P*V_E/(pi*D_E*R_E^2*N_A)*p_inv
end
