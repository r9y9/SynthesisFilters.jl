import MelGeneralizedCepstrums: allpass_alpha, glog_gamma

abstract type Filter end
abstract type SynthesisFilter <: Filter end

### Generic interfaces ###

# Note that f has state and filt! should update state value in this function
# TODO(ryuichi): separate state and filt definition?
function filt!(f::SynthesisFilter, x::Real, coef::AbstractVector)
    error("should provide sample-by-sample filt!")
end

# A synthesis filter must provide conversion from spectral parameter
# (e.g. mel-cepstrum) to filter coefficients.
function to_filtcoef(f::SynthesisFilter, param::AbstractArray)
    error("cannot determine how to convert paramters")
end

### Interfaces for mel-generalized ceptrum synthesis filters ###

abstract type MelGeneralizedCepstrumSynthesisFilter <: SynthesisFilter end

function allpass_alpha(f::MelGeneralizedCepstrumSynthesisFilter)
    error("should provide get access to allpass constant")
end

function glog_gamma(f::MelGeneralizedCepstrumSynthesisFilter)
    error("should provide get access to glog gamma")
end

# LPC, PARCOR and LSP
abstract type LinearPredictionVariantsSynthesisFilter <: SynthesisFilter end

### Waveform synthesis implementation ###

# synthesis_one_frame! generates a one frame speech waveform given a excitation
# signal and successive filter coefficients of a synthesis filter.
function synthesis_one_frame!(f::SynthesisFilter,
                              y::AbstractVector, # result will be stored
                              excitation::AbstractVector,
                              bᵗ⁻¹::AbstractVector,
                              bᵗ::AbstractVector)
    @assert length(y) == length(excitation)
    slope = (bᵗ - bᵗ⁻¹) / length(excitation)

    interpolated_coef = copy(bᵗ⁻¹)

    for i in 1:length(excitation)
        scaled_excitation = excitation[i] * exp(interpolated_coef[1])
        y[i] = filt!(f, scaled_excitation, interpolated_coef)
        for j=1:length(slope)
            interpolated_coef[j] += slope[j]
        end
    end

    y
end

function synthesis_one_frame!(f::SynthesisFilter,
                              excitation::AbstractVector,
                              bᵗ⁻¹::AbstractVector,
                              bᵗ::AbstractVector)
    synthesis_one_frame!(f, similar(excitation), excitation, bᵗ⁻¹, bᵗ)
end

# synthesis! generates a speech waveform given a excitation signal and
# a sequence of spectral envelope parameter.
function synthesis!(f::SynthesisFilter, # filter states are will be updated
                    excitation::AbstractVector,
                    b::AbstractMatrix,
                    hopsize::Integer)
    synthesized = similar(excitation)
    fill!(synthesized, zero(eltype(synthesized)))

    bᵗ⁻¹ = b[:,1]
    bᵗ = similar(bᵗ⁻¹)
    buf = Vector{eltype(synthesized)}(hopsize)

    for i in 1:size(b, 2)
        if i > 1
            bᵗ⁻¹ = b[:,i-1]
        end
        bᵗ = b[:,i]

        s, e = (i-1)*hopsize+1, i*hopsize
        e > length(excitation) && break

        synthesis_one_frame!(f, buf, excitation[s:e], bᵗ⁻¹, bᵗ)
        copy!(synthesized, s, buf, 1, length(buf))
    end

    synthesized
end

### Dispatch on specific filter types ###

function synthesis!{T<:MelGeneralizedCepstrum}(f::MelGeneralizedCepstrumSynthesisFilter,
                    excitation::AbstractVector,
                    state::SpectralParamState{T},
                    hopsize::Integer)
    synthesis!(f, excitation, to_filtcoef(f, state), hopsize)
end

function synthesis!{T<:LinearPredictionCoefVariants}(f::LinearPredictionVariantsSynthesisFilter,
                    excitation::AbstractVector,
                    state::SpectralParamState{T},
                    hopsize::Integer)
    synthesis!(f, excitation, to_filtcoef(f, state), hopsize)
end
