

(model::BatchModel)(du::AbstractVector{<:Real}, u::AbstractVector{<:Real}, p::AbstractVector{<:Real}, t::Real) = begin
    D, E, V, X, CO2 = split_state(du, p)
end
