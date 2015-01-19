module SynthesisFilters

import MelGeneralizedCepstrums: mc2b, mgc2b

export
  MelGeneralizedCepstrumSynthesisFilter, # LMADF or MLSADF or MGLSADF

  # Speech waveform sysnthesis filters
  LMADF,                # Log Magnitude Approximation Digital Filter
  MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
  MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter

  allpass_alpha,        # all-pass constant (alpha)
  glog_gamma,           # parameter of generalized log function
  nstage,

  # high-level interface for waveform synthesis
  synthesis_one_frame!, #
  synthesis!,           #

  # low-level interface
  filter!,              # filtering one sample
  mgc2b                 # mel-generalized cepstrum to filter coef.

for fname in [
              "common",
              "lmadf",
              "mlsadf",
              "mglsadf",
              "synthesis",
    ]
    include(string(fname, ".jl"))
end

end # module
