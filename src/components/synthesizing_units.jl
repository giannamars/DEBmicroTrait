function parallel_complementary_su(p::AbstractMetabolismCN, v)
    return 1.0./(sum(1.0./v) - 1.0./sum(v))
end

function parallel_complementary_su(p::AbstractMetabolismCNP, v)
    return 1.0./(sum(1.0./v)-1.0./(v[1]+v[2])-1.0./(v[1]+v[3])-1.0./(v[2]+v[3])+1.0./sum(v))
end
