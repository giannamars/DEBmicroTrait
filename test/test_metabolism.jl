using DEBmicroTrait, Test

## MetabolismC
### 1 microbe
k_E = 0.2*ones(1)
y_EV = 1.0*ones(1)
k_M = 0.1*ones(1)
y_EM = 1.0*ones(1)
α_X = 0.0*ones(1)
y_EX = 1.0*ones(1)

p = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)
E = 1.0*ones(1)
V = 1.0*ones(1)

r = DEBmicroTrait.growth!(0.0*ones(1), p, E, V)
r_analytic = @. (k_E*E/V - k_M*y_EM)/(E/V + y_EV)
@test r ≈ r_analytic atol=1e-3
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p, E, V)
@test x == [0.0]

α_X = 0.1*ones(1)
y_EX = 1.0*ones(1)
p = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)
E = 1.0*ones(1)
V = 1.0*ones(1)
r = DEBmicroTrait.growth!(0.0*ones(1), p, E, V)
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p, E, V)

### 2 microbes
k_E = 0.2*ones(2)
y_EV = 1.0*ones(2)
k_M = 0.1*ones(2)
y_EM = 1.0*ones(2)
α_X = 0.0*ones(2)
y_EX = 1.0*ones(2)

p = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)
E = 1.0*ones(2)
V = 1.0*ones(2)

r = DEBmicroTrait.growth!(0.0*ones(2), p, E, V)
r_analytic = @. (k_E*E/V - k_M*y_EM)/(E/V + y_EV)
@test r ≈ r_analytic atol=1e-3

α_X = 0.1*ones(2)
y_EX = 1.0*ones(2)
p = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)
E = 1.0*ones(2)
V = 1.0*ones(2)
r = DEBmicroTrait.growth!(0.0*ones(2), p, E, V)
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p, E, V)

### case: k_M*y_EM > j_EC
k_E = 0.2*ones(1)
y_EV = 1.0*ones(1)
k_M = 1.0*ones(1)
y_EM = 1.0*ones(1)
α_X = 0.1*ones(1)
y_EX = 1.0*ones(1)

p = MetabolismC(k_E, y_EV, k_M, y_EM, α_X, y_EX)
E = 1.0*ones(1)
V = 1.0*ones(1)

r = DEBmicroTrait.growth!(0.0*ones(1), p, E, V)
r_analytic = @. (k_E*E/V - k_M*y_EM)/(E/V + y_EV)
@test r[1] ≈ r_analytic[1] atol=1e-3
x, rG_CO2, rM_CO2, rX_CO2 = DEBmicroTrait.growth_production!(r, p, E, V)
@test x == [0.0]

### Metabolism CN
### 1 microbe
k_E = 0.2*ones(1)
y_EV = hcat([1.0, 1.0])
k_M = 0.1*ones(1)
y_EM = hcat([1.0, 1.0])
κ = hcat([1.0, 1.0])
y_EX = hcat([1.0, 1.0])

p = MetabolismCN(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = hcat([1.0, 1.0])
V = 1.0*ones(1)
r0 = 0.0*ones(1)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)
@test size(r,1) == 1
@test size(rM_CO2) == (2,1)
@test size(rEG_rej) == size(J_EE) == (2,1)
@test size(J_X,1) == 1
### 2 microbes
k_E = 0.2*ones(2)
y_EV = fill(1.0, 2, 2)
k_M = 0.1*ones(2)
y_EM = fill(1.0, 2, 2)
κ = fill(1.0, 2, 2)
y_EX = fill(1.0, 2, 2)

p = MetabolismCN(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = fill(1.0, 2, 2)
V = 1.0*ones(2)
r0 = 0.0*ones(2)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)
@test size(r,1) == 2
@test size(rM_CO2) == (2,2)
@test size(rEG_rej) == size(J_EE) == (2,2)
@test size(J_X,1) == 2

### 3 microbes
k_E = 0.2*ones(3)
y_EV = fill(1.0, 2, 3)
k_M = 0.1*ones(3)
y_EM = fill(1.0, 2, 3)
κ = fill(1.0, 2, 3)
y_EX = fill(1.0, 2, 3)

p = MetabolismCN(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = fill(1.0, 2, 3)
V = 1.0*ones(3)
r0 = 0.0*ones(3)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)

@test size(r,1) == 3
@test size(rM_CO2) == (2,3)
@test size(rEG_rej) == size(J_EE) == (2,3)
@test size(J_X,1) == 3

### Metabolism CNP
### 1 microbe

k_E = 0.2*ones(1)
y_EV = hcat([0.1, 0.2, 1.0])
k_M = 0.1*ones(1)
y_EM = hcat([1.0, 1.0, 1.0])
κ = hcat([1.0, 1.0, 1.0])
y_EX = hcat([1.0, 1.0, 1.0])

p = MetabolismCNP(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = hcat([1.0, 1.0, 1.0])
V = 1.0*ones(1)
r0 = 0.0*ones(1)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)

@test size(r,1) == 1
@test size(rM_CO2) == (3,1)
@test size(rEG_rej) == size(J_EE) == (3,1)
@test size(J_X,1) == 1


## 2 microbes
k_E = 0.2*ones(2)
y_EV = fill(1.0, 3, 2)
k_M = 0.1*ones(2)
y_EM = fill(1.0, 3, 2)
κ = fill(1.0, 3, 2)
y_EX = fill(1.0, 3, 2)

p = MetabolismCNP(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = fill(1.0, 3, 2)
V = 1.0*ones(2)
r0 = 0.0*ones(2)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)

@test size(r,1) == 2
@test size(rM_CO2) == (3,2)
@test size(rEG_rej) == size(J_EE) == (3,2)
@test size(J_X,1) == 2

## 3 microbes
k_E = 0.2*ones(3)
y_EV = fill(1.0, 3, 3)
k_M = 0.1*ones(3)
y_EM = fill(1.0, 3, 3)
κ = fill(1.0, 3, 3)
y_EX = fill(1.0, 3, 3)

p = MetabolismCNP(k_E, y_EV, k_M, y_EM, κ, y_EX)
E = fill(1.0, 3, 3)
V = 1.0*ones(3)
r0 = 0.0*ones(3)
r = DEBmicroTrait.growth!(r0, p, E, V)
rM_CO2, rEG_rej = DEBmicroTrait.growth_production!(r, p, E, V)
J_X = DEBmicroTrait.enzyme_production!(rEG_rej, p, V)
J_EE = DEBmicroTrait.reserve_feedback!(rEG_rej, p, V)

@test size(r,1) == 3
@test size(rM_CO2) == (3,3)
@test size(rEG_rej) == size(J_EE) == (3,3)
@test size(J_X,1) == 3
