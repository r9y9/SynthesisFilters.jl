__precompile__()

module SynthesisFilters

using Compat
import Compat: view
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
    AllPoleDF,            # All-pole digital filter for synthesis from LPC
    AllPoleLatticeDF,     # All-pole lattice digital filter for synthesis from PARCOR
    LSPDF,                # LSP digital filter for synthesis from LSP
    LMADF,                # Log magnitude approximation digital filter for synthesis from cepstrum
    MLSADF,               # Mel-log spectrum approximation digital filter for synthesis from mel-cepstrum
    MGLSADF,              # Mel generalized-log spectrum approximation digital filter for synthesis from mel-generalized cepstrum

    # High-level interface for waveform synthesis
    synthesis,
    synthesis!,
    synthesis_one_frame!,

    filt!,                # filtering one sample
    to_filtcoef,          # spectral parameter to filter coefficients

    # Mel-generalized cepstrum synthesis filter properties
    allpass_alpha,        # all-pass constant (α)
    glog_gamma,           # parameter of generalized log function (γ)

    # MGLSADF property
    nstage                # -1/γ; number of stages for MGLSADF

for fname in [
              "const",
              "common",
              "poledf",
              "ltcdf",
              "lspdf",
              "lmadf",
              "mlsadf",
              "mglsadf",
              "helper"
    ]
    include(string(fname, ".jl"))
end

end # module
