module SynthesisFilters

export
  SynthesisFilter,
  MelGeneralizedSynthesisFilter, # MLSADF or MGLSADF
  MLSADF,               # Mel-Log Spectrum Approximation Digital Filter
  MGLSADF,              # Mel Generalized-Log Spectrum Approximation Digital Filter
  synthesis_one_frame!, #
  synthesis!,           #
  filter!,              # filtering one sample (low-level)

  # Conversions (should be moved to other package?)
  gnorm,
  mc2b

abstract Filter
abstract SynthesisFilter <: Filter
abstract MelGeneralizedSynthesisFilter <: SynthesisFilter

function filtercoef_from_mgc(f::SynthesisFilter, mgc::Vector{Float64})
    error("not implemented")
end

for fname in [
              "gnorm",
              "mc2b",
              "mlsadf",
              "mglsadf",
              "synthesis",
    ]
    include(string(fname, ".jl"))
end

end # module
