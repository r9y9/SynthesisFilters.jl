# LSP speech synthesis digital filter
mutable struct LSPDF <: LinearPredictionVariantsSynthesisFilter
    delay::Vector{Float64}

    function LSPDF(order)
        new(SPTK.lspdf_delay(order))
    end
end

function filt!(f::LSPDF, x::Float64, coef::Vector{Float64})
    SPTK.lspdf(x, coef, delay(f))
end

function to_filtcoef(f::LSPDF, c::SpectralParamState{LineSpectralPair})
    b = copy(c)
    rawdata(loggain!(b))
end
