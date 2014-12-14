function mgc2b(f::MLSADF, mc::Vector{Float64})
    α = alpha(f)
    mc2b(mc, α)
end

function mgc2b(f::MGLSADF, mgc::Vector{Float64})
    α = alpha(f)
    γ = gamma(f)
    mgc2b(mc, α, γ)
end

# synthesis_one_frame! generates speech waveform for one frame speech signal
# given a excitation signal and successive two mel generalized cepstrum.
function synthesis_one_frame!(f::MelGeneralizedCepstrumSynthesisFilter,
                              excite::Vector{Float64},
                              previous_mgc::Vector{Float64},
                              current_mgc::Vector{Float64})
    previous_coef = mgc2b(f, previous_mgc)
    current_coef = mgc2b(f, current_mgc)

    slope = (current_coef - previous_coef) / float(length(excite))

    part_of_speech = Array(eltype(excite), length(excite))
    interpolated_coef = copy(previous_coef)

    for i=1:endof(excite)
        scaled_excitation = excite[i] * exp(interpolated_coef[1])
        part_of_speech[i] = filter!(f, scaled_excitation, interpolated_coef)
        interpolated_coef += slope
    end

    part_of_speech
end

# synthesis! generates a speech waveform given a excitation signal and
# a sequence of mel generalized cepstrum.
function synthesis!(f::MelGeneralizedCepstrumSynthesisFilter,
                    excite::Vector{Float64},
                    mgc_sequence::Matrix{Float64},
                    hopsize::Int)
    const T = length(excite)
    synthesized = zeros(T)

    previous_mgc = mgc_sequence[:,1]
    for i=1:size(mgc_sequence, 2)
        if i > 1
            previous_mgc = mgc_sequence[:,i-1]
        end
        current_mgc = mgc_sequence[:,i]

        const s, e = (i-1)*hopsize+1, i*hopsize
        if e > T
            break
        end

        part_of_speech = synthesis_one_frame!(f, excite[s:e],
                                              previous_mgc,
                                              current_mgc)
        synthesized[s:e] = part_of_speech
    end

    synthesized
end
