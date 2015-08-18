module SynthesisFilters

using MelGeneralizedCepstrums

import SPTK
import Base: filt!
import MelGeneralizedCepstrums: allpass_alpha, glog_gamma

export
    # Synthesis filter types
    SynthesisFilter,
    MelGeneralizedCepstrumSynthesisFilter, # LMADF or MLSADF or MGLSADF
    LinearPredictionVariantsSynthesisFilter, # AllPoleDF

    # Speech waveform synthesis filters
    AllPoleDF,            # All Pole Digital filter
    AllPoleLatticeDF,     # All Pole Lattice Digital filter
    LSPDF,                # Line Spectral Pair  Digital Filter
    LMADF,                # Log Magnitude Approximation Digital Filter
    MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
    MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter

    # high-level interface for waveform synthesis
    synthesis!,           #
    synthesis_one_frame!, #

    to_filtcoef,          # spectral parameter to filter coefficients
    filt!,                # filtering one sample

    allpass_alpha,        # all-pass constant (alpha)
    glog_gamma,           # parameter of generalized log function

    nstage                # -1/γ; number of stages for MGLSADF

for fname in [
              "const",
              "common",
              "poledf",
              "ltcdf",
              "lspdf",
              "lmadf",
              "mlsadf",
              "mglsadf"
    ]
    include(string(fname, ".jl"))
end

end # module
