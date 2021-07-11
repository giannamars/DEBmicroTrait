abstract type AbstractMetabolism end
abstract type AbstractMetabolismC <: AbstractMetabolism end
abstract type AbstractMetabolismCN <: AbstractMetabolism end
abstract type AbstractMetabolismCNP <: AbstractMetabolism end

@columns struct MetabolismC{KE,YEV,KM,YEM,AX,YEX,FAX,MGT} <: AbstractMetabolismC
    k_E::KE     | 1/h       | "Reserve export rate"
    y_EV::YEV   | mol/mol   | "Yield of structure on reserve"
    k_M::KM     | 1/h       | "Specific maintenance rate"
    y_EM::YEM   | mol/mol   | "Maintenance yield on reserve vs structure"
    α_X::AX     | _         | "Fraction of growth flux allocated to enzyme production"
    y_EX::YEX   | mol/mol   | "Yield on enzyme production"
    f_αX::FAX   | _         | "Fractional allocation to different enzymes"
    min_gt::MGT | 1/h       | "Minimum generation time"
end

for fn in fieldnames(MetabolismC)
    @eval $fn(p::MetabolismC) = p.$fn
end

growth!(r0, p::MetabolismC, E::Vector{Float64}, V::Vector{Float64}) = begin
    y_VM = y_EM(p)./y_EV(p)

    function f(r, j)
      m_E   = max(1e-8, E[j]/V[j])
      j_EC  = m_E*(k_E(p)[j] - r)
      j_EM  = k_M(p)[j]*y_EM(p)[j]
      jEM   = max(1e-8, min(j_EC, j_EM))
      jVM   = (j_EM - jEM)*y_VM[j]/y_EM(p)[j]
      j_EG  = (j_EC - jEM)/(1.0 + α_X(p)[j])
      j_G   = j_EG/y_EV(p)[j]
      res   = r - (j_G - jVM)
    end

    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
    gmax = log(2)./min_gt(p)
    r_b = min.(r, gmax)
end

growth_production!(r, p::MetabolismC, E::Vector{Float64}, V::Vector{Float64}) = begin
      y_VM  = y_EM(p)./y_EV(p)

      m_E   = E./V
      j_EC  = m_E.*(k_E(p) .- r)
      j_EM  = k_M(p).*y_EM(p)
      jEM   = min(j_EC, j_EM)
      jVM   = (j_EM - jEM).*y_VM./y_EM(p)
      j_EG  = (j_EC - jEM)./(1.0 .+ α_X(p))
      x     = α_X(p).*j_EG./y_EX(p)

      rG_CO2  = (1 .- 1 ./y_EV(p)).*j_EG.*V
      rM_CO2  = @. (jEM + jVM)*V
      rX_CO2  = (1 .- 1 ./y_EX(p)).*α_X(p).*j_EG.*V
      return x, rG_CO2, rM_CO2, rX_CO2
end

enzyme_production!(x, p::MetabolismC, V::Vector{Float64}) = begin
    J_EX = f_αX(p).*sum(x.*V, dims=1)
end

@columns struct MetabolismCN{KE,YEV,KM,YEM,KAP,YEX} <: AbstractMetabolismCN
    k_E::KE         | 1/h       | "Reserve export rate"
    y_EV::YEV       | mol/mol   | "Yield of structure on C- and N-reserve"
    k_M::KM         | 1/h       | "Specific maintenance rate"
    y_EM::YEM       | mol/mol   | "Maintenance yield on C and N-reserve vs structure"
    κ::KAP          | _         | "Fraction of rejection flux allocated to enzyme production"
    y_EX::YEX       | mol/mol   | "Yield of enzyme production on reserves"
end

for fn in fieldnames(MetabolismCN)
    @eval $fn(p::MetabolismCN) = p.$fn
end

growth!(r0, p::MetabolismCN, E::Matrix{Float64}, V::Vector{Float64}) = begin
    y_VM = y_EM(p)./y_EV(p)

    function f(r, j)
      m_E     = @. max(1e-8, E[:,j]/V[j])
      j_EC    = m_E.*(k_E(p)[j] .- r)
      j_EM    = k_M(p)[j].*y_EM(p)[:,j]
      jEM     = @. max(1e-8, min(j_EC, j_EM))
      jVM     = (j_EM - jEM).*y_VM[:,j]./y_EM(p)[:,j]
      j_EG    = @. (j_EC - jEM)
      j_EG_V  = max.(1e-8, j_EG./y_EV(p)[:,j])
      j_G     = parallel_complementary_su(p, j_EG_V)
      #j_G     = 1.0./(sum(1.0./j_EG_V) - 1.0./sum(j_EG_V))
      res     = r - (j_G - sum(jVM))
    end

    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
end

growth_production!(r, p::MetabolismCN, E::Matrix{Float64}, V::Vector{Float64}) = begin
  y_VM    = y_EM(p)./y_EV(p)

  rM_CO2 = zeros(size(y_VM))
  rEG_rej = zeros(size(y_VM))

  for j in 1:size(r,1)
      m_E     = @. E[:,j]/V[j]
      j_EC    = m_E.*(k_E(p)[j] .- r[j])
      j_EM    = k_M(p)[j].*y_EM(p)[:,j]
      jEM     = @. min(j_EC, j_EM)
      jVM     = (j_EM .- jEM).*y_VM[:,j]./y_EM(p)[:,j]

      rM_CO2[:,j]  = @. (jEM + jVM)*V[j]
      rEG_rej[:,j] = (k_E(p)[j] .- r[j]).*m_E .- jEM .- y_EV(p)[:,j].*r[j]
  end
  return rM_CO2, rEG_rej
end

enzyme_production!(rEG_rej, p::MetabolismCN, V::Vector{Float64}) = begin
    rEG_X = κ(p).*rEG_rej
    rEG_X_V = rEG_X./y_EX(p)
    x = parallel_complementary_su(p, rEG_X_V)
    return x.*V
end

reserve_feedback!(rEG_rej, p::MetabolismCN, V::Vector{Float64}) = begin
    rEG_E = (1.0 .- κ(p)).*rEG_rej
end


@columns struct MetabolismCNP{KE,YEV,KM,YEM,KAP,YEX} <: AbstractMetabolismCNP
    k_E::KE         | 1/h       | "Reserve export rate"
    y_EV::YEV       | mol/mol   | "Yield of structure on C- and N-reserve"
    k_M::KM         | 1/h       | "Specific maintenance rate"
    y_EM::YEM       | mol/mol   | "Maintenance yield on C and N-reserve vs structure"
    κ::KAP          | _         | "Fraction of rejection flux allocated to enzyme production"
    y_EX::YEX       | mol/mol   | "Yield of enzyme production on reserves"
end

for fn in fieldnames(MetabolismCNP)
    @eval $fn(p::MetabolismCNP) = p.$fn
end

growth!(r0, p::MetabolismCNP, E::Matrix{Float64}, V::Vector{Float64}) = begin
    y_VM = y_EM(p)./y_EV(p)

    function f(r, j)
      m_E     = @. max(1e-8, E[:,j]/V[j])
      j_EC    = m_E.*(k_E(p)[j] .- r)
      j_EM    = k_M(p)[j].*y_EM(p)[:,j]
      jEM     = @. max(1e-8, min(j_EC, j_EM))
      jVM     = (j_EM - jEM).*y_VM[:,j]./y_EM(p)[:,j]
      j_EG    = @. (j_EC - jEM)
      j_EG_V  = max.(1e-8, j_EG./y_EV(p)[:,j])
      j_G     = parallel_complementary_su(p, j_EG_V)
      #j_G     = 1.0./(sum(1.0./j_EG_V)-1.0./(j_EG_V[1]+j_EG_V[2])-1.0./(j_EG_V[1]+j_EG_V[3])-1.0./(j_EG_V[2]+j_EG_V[3])+1.0./sum(j_EG_V))
      res     = r - (j_G - sum(jVM))
    end

    f1 = [r->f(r,j) for j in 1:size(r0,1)]
    r = Roots.find_zero.(f1, r0)
end


growth_production!(r, p::MetabolismCNP, E::Matrix{Float64}, V::Vector{Float64}) = begin
  y_VM    = y_EM(p)./y_EV(p)

  rM_CO2 = zeros(size(y_VM))
  rEG_rej = zeros(size(y_VM))

  for j in 1:size(r,1)
      m_E     = @. E[:,j]/V[j]
      j_EC    = m_E.*(k_E(p)[j] .- r[j])
      j_EM    = k_M(p)[j].*y_EM(p)[:,j]
      jEM     = @. min(j_EC, j_EM)
      jVM     = (j_EM .- jEM).*y_VM[:,j]./y_EM(p)[:,j]

      rM_CO2[:,j]  = @. (jEM + jVM)*V[j]
      rEG_rej[:,j] = (k_E(p)[j] .- r[j]).*m_E .- jEM .- y_EV(p)[:,j].*r[j]
  end
  return rM_CO2, rEG_rej
end


enzyme_production!(rEG_rej, p::MetabolismCNP, V::Vector{Float64}) = begin
    rEG_X = κ(p).*rEG_rej
    rEG_X_V = rEG_X./y_EX(p)
    x = parallel_complementary_su(p, rEG_X_V)
    return x.*V
end

reserve_feedback!(rEG_rej, p::MetabolismCNP, V::Vector{Float64}) = begin
    rEG_E = (1.0 .- κ(p)).*rEG_rej
end
