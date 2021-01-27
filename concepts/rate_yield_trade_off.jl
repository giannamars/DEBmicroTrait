using DEBmicroTrait
using Roots
using Gadfly
using JLD


rrn_copies = DEBmicroTrait.
gmax = DEBmicroTrait.gmax_regression(rrn_copies)

Min_gen_time = 2.0*ones(1)
Gram_stain = ["+"]
α = 0.0

V_cs = [1.05e-20, 1e-19, 1e-18, 1e-17, 1e-16]
r_rates = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]
yield = zeros(size(V_cs,1)+1, size(r_rates,1)+1)
rate = zeros(size(V_cs,1)+1, size(r_rates,1)+1)

for i in 1:size(V_cs,1)
    find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[i]*ones(1), Min_gen_time, Gram_stain)
    ρ_p = Roots.find_zero(find_ρ, 1.0)
    yield[i,:], rate[i,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α, V_cs[i]*ones(1), Min_gen_time, Gram_stain, r_rates)
end

α = 1.0
V_cs = [1.05e-20]
find_ρ(x) = DEBmicroTrait.constrain_transporter_density(x, V_cs[1]*ones(1), Min_gen_time, Gram_stain)
ρ_p = Roots.find_zero(find_ρ, 1.0)
yield[6,:], rate[6,:] = DEBmicroTrait.rate_yield_trade_off(ρ_p, α, V_cs[1]*ones(1), Min_gen_time, Gram_stain, r_rates)


FCR = zeros(5)
for i in 1:5
    rmax, rid = findmax(rate[i,:])
    Ymax, Yid = findmax(yield[i,:])
    FCR[i] = Ymax - yield[i,rid]
end

save("/Users/glmarschmann/.julia/dev/DEBmicroTrait/plots/files/rate_yield_trade_off_1.jld", "yield", yield, "rate", rate, "FCR", FCR, "Vc", V_cs)
