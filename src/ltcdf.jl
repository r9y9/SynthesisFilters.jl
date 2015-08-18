# All Pole lattice digital filter
type AllPoleLatticeDF <: LinearPredictionVariantsSynthesisFilter
    delay::Vector{Float64}

    function AllPoleLatticeDF(order)
        new(SPTK.ltcdf_delay(order))
    end
end

function filt!(f::AllPoleLatticeDF, x::Float64, coef::Vector{Float64})
    SPTK.ltcdf(x, coef, delay(f))
end

function to_filtcoef(f::AllPoleLatticeDF,
                     c::SpectralParamState{PartialAutoCorrelation})
    b = copy(c)
    rawdata(loggain!(b))
end
