abstract Filter
abstract SynthesisFilter <: Filter

# A synthesis filter must provide conversion from spectral parameter
# (e.g. mel-cepstrum) to coefficients of a synthesis filter.
function to_filtcoef(f::SynthesisFilter, param::AbstractArray)
    error("no conversion implementation found")
end

abstract MelGeneralizedCepstrumSynthesisFilter <: SynthesisFilter

allpass_alpha(f::MelGeneralizedCepstrumSynthesisFilter) = error("not implemented")
glog_gamma(f::MelGeneralizedCepstrumSynthesisFilter) = error("not implemented")

abstract MelLinearPredictionSynthesisFilter <: SynthesisFilter

allpass_alpha(f::MelLinearPredictionSynthesisFilter) = error("not implemented")

# synthesis_one_frame! generates a one frame speech waveform given a excitation
# signal and successive filter coefficients of a synthesis filter.
function synthesis_one_frame!(f::SynthesisFilter,
                              excitation::AbstractVector,
                              bᵗ⁻¹::Vector,
                              bᵗ::Vector)
    slope = (bᵗ - bᵗ⁻¹) / length(excitation)

    y = Array(eltype(excitation), length(excitation))
    interpolated_coef = copy(bᵗ⁻¹)

    for i=1:length(excitation)
        scaled_excitation = excitation[i] * exp(interpolated_coef[1])
        y[i] = filter!(f, scaled_excitation, interpolated_coef)
        for j=1:length(slope)
            interpolated_coef[j] += slope[j]
        end
    end

    y
end

# synthesis! generates a speech waveform given a excitation signal and
# a sequence of spectral envelope paramter.
function synthesis!(f::SynthesisFilter,
                    excitation::AbstractVector,
                    b::Matrix,
                    hopsize::Integer)
    synthesized = similar(excitation)
    fill!(synthesized, zero(eltype(synthesized)))

    bᵗ⁻¹ = b[:,1]
    for i=1:size(b, 2)
        if i > 1
            @inbounds bᵗ⁻¹ = b[:,i-1]
        end
        @inbounds bᵗ = b[:,i]

        s, e = (i-1)*hopsize+1, i*hopsize
        if e > length(excitation)
            break
        end

        x = synthesis_one_frame!(f, excitation[s:e], bᵗ⁻¹, bᵗ)
        copy!(synthesized, s, x, 1, length(x))
    end

    synthesized
end

function synthesis!(f::MelGeneralizedCepstrumSynthesisFilter,
                    excitation::AbstractVector,
                    mgc::MelGeneralizedCepstrum,
                    hopsize::Integer)
    synthesis!(f, excitation, to_filtcoef(f, mgc), hopsize)
end

function synthesis!(f::MelLinearPredictionSynthesisFilter,
                    excitation::AbstractVector,
                    mlpc::MelLinearPredictionCoef,
                    hopsize::Integer)
    synthesis!(f, excitation, to_filtcoef(f, mlpc), hopsize)
end
