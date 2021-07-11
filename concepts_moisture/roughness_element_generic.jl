using DEBmicroTrait
using Gadfly

matric_potential = convert(Array{Float64,1}, range(0.01, 100, length=100))
d = DEBmicroTrait.effective_film_thickness(-matric_potential)
D_star = DEBmicroTrait.max_immersed_diameter(-matric_potential)

#
v_0 = DEBmicroTrait.bulk_velocity(cell_radius)
matric_potentials = convert(Array{Float64,1}, range(0.01, 100, length=100))

plot(x=matric_potentials, y=D_star, Scale.y_log10)
d = DEBmicroTrait.effective_film_thickness(-matric_potential)
λ = DEBmicroTrait.lamba_drag(cell_radius, d)

F_M = DEBmicroTrait.propulsion_force(cell_radius)
F_λ = DEBmicroTrait.hydrodynamic_force(cell_radius, d)
F_c = DEBmicroTrait.capillary_force(cell_radius, d, -matric_potential)
vs = zeros(100)
for i in 1:100
    v = DEBmicroTrait.cell_velocity(cell_radius, d, -[matric_potentials[i]])
    vs[i] = v[1]
end

vs

A_s = 1.9e-19
ρ = 998
σ = 0.07275
I_P = @. (A_s/(6*pi*ρ*P))^(1/3)
r_P = @. -σ/(ρ*P)
γ = 4
H = 0.5*1e-3
θ = 130
d_var = @. I_P*(γ*H + 2*(H/cos(θ/2)) - r_P/tan(θ/2))/(H*(γ + 2*tan(θ/2)))
