module SynthesisFilters

export
  SynthesisFilter,
  MelGeneralizedSynthesisFilter, # MLSADF or MGLSADF
  MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
  MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter
  synthesis_one_frame!, #
  synthesis!,           #
  filter!               # filtering one sample (low-level)

abstract Filter
abstract SynthesisFilter <: Filter
abstract MelGeneralizedSynthesisFilter <: SynthesisFilter

for fname in [
              "mlsadf",
              "mglsadf",
              "synthesis",
    ]
    include(string(fname, ".jl"))
end

end # module
