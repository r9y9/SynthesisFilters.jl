# high-level interface

### Synthesis for LPC variants ###

function synthesis(excitation::AbstractVector,
                   param::SpectralParamState{LinearPredictionCoef},
                   hopsize::Integer)
    f = AllPoleDF(param_order(paramdef(param)))
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end

function synthesis(excitation::AbstractVector,
                   param::SpectralParamState{PartialAutoCorrelation},
                   hopsize::Integer)
    f = AllPoleLatticeDF(param_order(paramdef(param)))
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end

function synthesis(excitation::AbstractVector,
                   param::SpectralParamState{LineSpectralPair},
                   hopsize::Integer)
    f = LSPDF(param_order(paramdef(param)))
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end

### Synthesis for mel-generalized cepstrums

function synthesis{T<:MelGeneralizedCepstrum}(excitation::AbstractVector,
                                              param::SpectralParamState{T},
                                              hopsize::Integer)
    def = paramdef(param)
    ns = try
        # InexactError could happen
        convert(Int, -1/glog_gamma(def))
    catch e
        ns = -1/glog_gamma(def)
        error("-1/γ must be able to convert to Int, but got $(ns)")
    end
    f = MGLSADF(param_order(def), allpass_alpha(def), ns)
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end

function synthesis(excitation::AbstractVector,
                   param::SpectralParamState{MelCepstrum},
                   hopsize::Integer)
    def = paramdef(param)
    f = MLSADF(param_order(def), allpass_alpha(def))
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end

function synthesis(excitation::AbstractVector,
                   param::SpectralParamState{LinearCepstrum},
                   hopsize::Integer)
    f = LMADF(param_order(paramdef(param)))
    synthesis!(f, excitation, to_filtcoef(f, param), hopsize)
end
