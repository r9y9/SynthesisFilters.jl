module SynthesisFilters

using MelGeneralizedCepstrums

import Base: filt!

export
    # Abstract types
    SynthesisFilter,
    MelGeneralizedCepstrumSynthesisFilter, # LMADF or MLSADF or MGLSADF
    LinearPredictionVariantsSynthesisFilter, # AllPoleDF

    # Speech waveform synthesis filters
    AllPoleDF,            # All Pole Digital filter
    LMADF,                # Log Magnitude Approximation Digital Filter
    MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
    MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter

    # properties
    allpass_alpha,        # all-pass constant (alpha)
    glog_gamma,           # parameter of generalized log function
    nstage,

    # high-level interface for waveform synthesis
    synthesis_one_frame!, #
    synthesis!,           #

    # low-level interface
    filt!,                # filtering one sample
    to_filtcoef           # spectral parameter to filter coefficients

for fname in [
              "const",
              "common",
              "poledf",
              "lmadf",
              "mlsadf",
              "mglsadf"
    ]
    include(string(fname, ".jl"))
end

end # module
