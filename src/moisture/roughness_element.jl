function effective_film_thickness(P::Vector{Float64})
    A_s = -1.9e-19
    ρ = 998
    σ = 0.07275
    I_P = @. (A_s/(6*pi*ρ*P))^(1/3)
    r_P = @. -σ/(ρ*P)
    γ = 4
    H = 0.5*1e-5
    θ = 130
    d_var = @. (I_P*(γ*H + 2*(H/cos(θ/2)) - r_P/tan(θ/2)))/(H*(γ + 2*tan(θ/2)))
end

function max_immersed_diameter(P::Vector{Float64})
    σ = 0.07275 #[kg/s], Or&Tuller 2000
    θ = 130
    angle = (1-sin(θ/2))/(1+sin(θ/2))
    D_star =  @. -4*σ/(P*1000)*angle
end

function lamba_drag(R::Vector{Float64}, d::Vector{Float64})
    λ_p = @. (1-9/16*(R/d) + 1/8*(R/d)^3)^(-1)
    λ_n = @. 1+9/8*(R/d) + (9*R/(8*d))^2
    λ = @. sqrt(λ_n^2 + λ_p^2)
end

function bulk_velocity(R::Vector{Float64})
    v_0 = @. 38.542*(2*R)^0.5424  # Flynn 2018 [μm/s]
end

function propulsion_force(R::Vector{Float64})
    η = 0.001002 #kg/(ms)
    v_0 = bulk_velocity(R)
    F_M = @. 6*π*η*R*v_0*1e-6
end

function hydrodynamic_force(R::Vector{Float64}, d::Vector{Float64})
    F_M = propulsion_force(R)
    λ = lamba_drag(R, d)
    F_λ = @. (1-1/λ)*F_M
end



function capillary_force(R::Vector{Float64}, d::Vector{Float64}, P::Vector{Float64})
        σ = 0.07275 #[kg/s], Or&Tuller 2000
        δ = 0.01
        r_c = max_immersed_diameter(P)./2
        F_c = zeros(size(r_c,1))
        for i in 1:size(F_c,1)
            if r_c[i] > R[i]
                F_c[i] = 0.0
            else
                tmp =  @. (2*π*σ/R-π*P)*r_c^4*δ
                F_c[i] = tmp[1]
            end
        end
        return F_c
end

function cell_velocity(R::Vector{Float64},d::Vector{Float64}, P::Vector{Float64})
    v_0 = @. 38.542*(2*R)^0.5424  # Flynn 2018 [μm/s]
    F_M = propulsion_force(R)
    F_λ = hydrodynamic_force(R,d)
    F_c = capillary_force(R,d,P)
    tmp = @. F_M - F_λ - F_c
    v = zeros(size(tmp,1))
    for i in 1:size(tmp,1)
        if tmp[i] < 0
            v[i] = 0
        else
            tmp1 = @. v_0*tmp/F_M
            v[i] = tmp1[1]
        end
    end
    return v
end


#function colony_expansion_rate(ξ::Vector{Float64},gmax::Vector{Float64},v::Vector{Float64},N_p::Vector{Float64},K::Vector{Float64}, )
