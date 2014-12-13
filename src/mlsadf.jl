# Mel-Log Spectrum Approximation filter (MLSADF)

# MLSABaseFilter represents a base filter of the MLSADF.
type MLSABaseFilter <: Filter
    order::Int
    α::Float64
    delay::Vector{Float64}
    
    function MLSABaseFilter(order::Int, α::Float64)
        new(order, α, zeros(order+1))
    end
end

delay(f::Filter) = f.delay
alpha(f::MLSABaseFilter) = f.α

function filter!(bf::MLSABaseFilter, x::Float64, coef::Vector{Float64})
    d = delay(bf)
    d[1] = x
    d[2] = (1.0-bf.α*bf.α)*d[1] + bf.α*d[2]

    result = 0.0
    for i=3:length(coef)
        @inbounds d[i] = d[i] + bf.α*(d[i+1] - d[i-1])
        @inbounds result += d[i] * coef[i]
    end

    # special case
    if length(coef) == 2
        result += d[2] * coef[2]
    end

    # t <- t+1 in time
    for i=length(d):-1:3
        @inbounds d[i] = d[i-1] 
    end

    result
end

# MLSACascadeFilter represents a cascade filter which contains MLSA base 
# filters.
type MLSACascadeFilter <: Filter
    filters::Vector{MLSABaseFilter}
    padecoef::Vector{Float64}
    delay::Vector{Float64}

    function MLSACascadeFilter(order::Int, α::Float64, pade::Int)
        padecoef = Array(Float64, pade+1)

        if pade == 4
            padecoef = [1.0, 4.999273e-1, 1.067005e-1, 1.170221e-2, 5.656279e-4]
        elseif pade == 5
            padecoef = [1.0, 4.999391e-1, 1.107098e-1, 1.369984e-2, 9.564853e-4,
                        3.041721e-5]
        else
            error("MLSADF: Order of pade approximation 4 or 5 is only supported.")
        end

        filters = Array(MLSABaseFilter, pade+1)
        for i=1:length(filters)
            filters[i] = MLSABaseFilter(order, α)
        end

        new(filters, padecoef, zeros(pade+1))
    end
end

function filter!(cf::MLSACascadeFilter, x::Float64, coef::Vector{Float64})
    d = delay(cf)
    result, feedback = 0.0, 0.0

    for i=length(cf.padecoef):-1:2
        @inbounds d[i] = filter!(cf.filters[i], d[i-1], coef)
        @inbounds val = d[i] * cf.padecoef[i]
        if iseven(i)
            feedback += val
        else
            feedback -= val
        end
        result += val
    end

    d[1] = feedback + x
    result += d[1]
    result
end

type MLSADF <: MelGeneralizedSynthesisFilter
    filters::Vector{MLSACascadeFilter}

    function MLSADF(order::Int, α::Float64; pade::Int=5)
        filters = Array(MLSACascadeFilter, 2)
        filters[1] = MLSACascadeFilter(2, α, pade)
        filters[2] = MLSACascadeFilter(order+1, α, pade)
        new(filters)
    end
end

first(f::MLSADF) = f.filters[1]
last(f::MLSADF) = f.filters[2]

function filter!(f::MLSADF, x::Float64, coef::Vector{Float64})
    filter!(last(f), filter!(first(f), x, [0.0, coef[2]]), coef)
end
