# Mel Generalized-Log Spectrum Approximation Digital Filter

# MGLSABaseFilter represents a base filter of MGLSAFilter.
type MGLSABaseFilter <: Filter
    order::Int   # order of mel generalized cepstrum
    α::Float64  # all-pass constant
    delay::Vector{Float64}

    function MGLSABaseFilter(order::Int, α::Float64)
        new(order, α, zeros(order+1))
    end
end

alpha(f::MGLSABaseFilter) = f.α

function filter!(bf::MGLSABaseFilter, x::Float64, coef::Vector{Float64})
    d = delay(bf)
    y = d[1] * coef[2]
    α = alpha(bf)

    for i=2:length(coef)-1
        @inbounds d[i] += α*(d[i+1] - d[i-1])
        @inbounds y += d[i] * coef[i+1]
    end

    result = x - y

    # t <- t+1 in time
    for i=length(d):-1:2
        @inbounds d[i] = d[i-1]
    end
    d[1] = α*d[1] + (1.0 - α*α)*result

    result
end

# MGLSAFilter represents a Mel Generalized Log Spectrum Approximation Digital
# Filter.
type MGLSADF <: MelGeneralizedSynthesisFilter
    filters::Vector{MGLSABaseFilter}

    function MGLSADF(order::Int, α::Float64, nstage::Int)
        filters = Array(MGLSABaseFilter, nstage)
        for i=1:length(filters)
            filters[i] = MGLSABaseFilter(order, α)
        end
        new(filters)
    end
end

alpha(f::MGLSADF) = alpha(f.filters[1])
gamma(f::MGLSADF) = -1.0/length(filters)

function filter!(f::MGLSADF, x::Float64, coef::Vector{Float64})
    y = x
    for i=1:length(f.filters)
        y = filter!(f.filters[i], y, coef)
    end
    y
end

function filtercoef_from_mgc(f::MGLSADF, mgc::Vector{Float64})
    α= alpha(f)
    γ = gamma(f)
    b = mc2b(mc, α)

    if γ == 0.0
        return b
    end

    b = gnorm(b, γ)

    # scale by gamma
    b[1] = log(b[1])
    b[2:end] *= γ

    b
end
