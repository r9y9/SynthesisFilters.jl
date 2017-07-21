# Log Magnitude Approximation (LMA) digital filter

mutable struct LMABaseFilter <: Filter
    order::Int
    delay::Vector{Float64}

    function LMABaseFilter(order::Int)
        new(order, zeros(order+1))
    end
end

function filt!(f::LMABaseFilter, x::Float64, coef::Vector{Float64},
               m1::Int, m2::Int)
    @assert length(coef) == f.order+1
    @assert m1 >= 1
    @assert m2 < length(coef)

    d = delay(f)

    for i=length(coef)-1:-1:2
        @inbounds d[i] = d[i-1]
    end
    d[1] = x

    result = 0.0
    for i=m1+1:m2+1
        @inbounds result += coef[i] * d[i-1]
    end

    result
end

# LMACascadeFilter represents a cascade filter which contains LMA base
# filters.
mutable struct LMACascadeFilter <: Filter
    filters::Vector{LMABaseFilter}
    padecoef::Vector{Float64}
    delay::Vector{Float64}

    function LMACascadeFilter(order::Int, pade::Int)
        padecoef = Vector{Float64}(pade+1)

        if pade == 4
            padecoef = padecoef4th
        elseif pade == 5
            padecoef = padecoef5th
        else
            error("LMADF: Order of pade approximation 4 or 5 is only supported.")
        end

        filters = Vector{LMABaseFilter}(pade+1)
        for i=1:length(filters)
            filters[i] = LMABaseFilter(order)
        end

        new(filters, padecoef, zeros(pade+1))
    end
end

padecoef(f::LMACascadeFilter) = f.padecoef

function filt!(f::LMACascadeFilter, x::Float64, coef::Vector{Float64},
                 m1::Int, m2::Int)
    d = delay(f)
    pade = padecoef(f)
    result, feedback = 0.0, 0.0

    for i=length(pade):-1:2
        @inbounds d[i] = filt!(f.filters[i], d[i-1], coef, m1, m2)
        @inbounds val = d[i] * pade[i]
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

mutable struct LMADF <: MelGeneralizedCepstrumSynthesisFilter
    filters::Vector{LMACascadeFilter}

    function LMADF(order::Int; pade::Int=5)
        filters = Vector{LMACascadeFilter}(2)
        filters[1] = LMACascadeFilter(order, pade)
        filters[2] = LMACascadeFilter(order, pade)
        new(filters)
    end
end

first(f::LMADF) = f.filters[1]
last(f::LMADF) = f.filters[2]
allpass_alpha(f::LMADF) = 0.0
glog_gamma(f::LMADF) = 0.0

function filt!(f::LMADF, x::Float64, coef::Vector{Float64})
    m = length(coef)-1
    y = filt!(first(f), x, coef, 1, 1)
    filt!(last(f), y, coef, 2, m)
end

function to_filtcoef(f::LMADF, c::SpectralParamState{LinearCepstrum})
    copy(rawdata(c))
end
