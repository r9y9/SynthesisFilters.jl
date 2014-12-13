function mc2b!(coef::Vector{Float64}, mc::Vector{Float64}, α::Float64)
    @assert length(coef) == length(mc)
    coef[end] = mc[end]
    for i=length(mc)-1:-1:1
        coef[i] = mc[i] - α*coef[i+1]
    end    
    coef
end

function mc2b(mc::Vector{Float64}, α::Float64)
    coef = Array(Float64, length(mc))
    mc2b!(coef, mc, α)
end
