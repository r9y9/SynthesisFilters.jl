module SynthesisFilters

import MelGeneralizedCepstrums: mc2b, mgc2b

export
  SynthesisFilter,
  MelGeneralizedCepstrumSynthesisFilter, # MLSADF or MGLSADF
  MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
  MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter
  allpass_alpha,        # all-pass constant (alpha)
  glog_gamma,           # parameter of generalized log function
  nstage,
  synthesis_one_frame!, #
  synthesis!,           #
  filter!               # filtering one sample (low-level)

abstract Filter
abstract SynthesisFilter <: Filter
abstract MelGeneralizedCepstrumSynthesisFilter <: SynthesisFilter

allpass_alpha(f::MelGeneralizedCepstrumSynthesisFilter) = error("not implemented")
glog_gamma(f::MelGeneralizedCepstrumSynthesisFilter) = error("not implemented")

for fname in [
              "mlsadf",
              "mglsadf",
              "synthesis",
    ]
    include(string(fname, ".jl"))
end

end # module
