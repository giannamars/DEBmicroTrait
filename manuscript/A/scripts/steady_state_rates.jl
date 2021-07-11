using DEBmicroTrait
using Gadfly

j_ED_max = 0.01
k_M      = 1e-4
y_EM     = 1.0
k_E      = 1.1
y_EV     = 1.4
y_EX     = y_EV
α_X      = range(0.01, 1.5, length=100)

r_ss     = @. (j_ED_max-k_M*y_EM)/(j_ED_max/k_E + (1+α_X)*y_EV)
x_ss     = @. α_X*y_EV*r_ss/y_EX

plot(x=x_ss, y=r_ss, Geom.line)
