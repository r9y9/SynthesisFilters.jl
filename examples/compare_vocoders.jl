# let's compare vocoders that are implemented in SynthesisFilters.jl
#
# After you've run the script, you'll get the following synthesized wav files:
# 1. test16k_poledf.wav
# 2. test16k_lmadf.wav
# 3. test16k_mlsadf.wav
# 4. test16k_mglsadf.wav
#
# The original wav file is found in data/test16k.wav.
#
# For more high quality speech waveform syntheis, you can check the `WORLD.jl`
# https://github.com/r9y9/WORLD.jl

using MelGeneralizedCepstrums
using SynthesisFilters
using WAV
using DSP

# Note about excitation
# fs: 16000
# frame period: 5.0 ms
# F0 analysis: esimated  by WORLD.dio and WORLD.stonemask
# Excitation genereration: perioic pulse for voiced segments and gaussian random
# value for un-voiced segments
base_excitation = vec(readdlm(joinpath(dirname(@__FILE__), "data",
                                       "test16k_excitation.txt")))

x, fs = wavread(joinpath(dirname(@__FILE__), "data", "test16k.wav"))
x = vec(x)
framelen = 512
hopsize = 80 # 5.0 ms for fs 16000
noverlap = framelen - hopsize

# Note that mgcep analysis basically assumes power-normalized window
# ŵ(n) = w(n)/|Σₙ w(n)²|¹/²
win = DSP.blackman(framelen) ./ sqrt(sumabs2(DSP.blackman(framelen)))

# create windowed signal matrix that each column represents a windowed time slice
as = arraysplit(x, framelen, noverlap)
xw = Array(Float64, framelen, length(as))
for t=1:length(as)
    xw[:,t] = as[t]
end
# col-wise windowing
xw .*= win

function test_poledf_synthesis(; order=25, savepath="test16k_poledf.wav")
    println("testing: poledf_synthesis")
    l = lpc(xw, order, use_mgcep=true)

    f = AllPoleDF(order)
    y = synthesis!(f, base_excitation, l, hopsize)
    wavwrite(y, savepath; Fs=fs)
    println("Dumped to $savepath")
end

function test_lmadf_synthesis(; order=25, savepath="test16k_lmadf.wav")
    println("testing: lmadf_synthesis")
    c = mgcep(xw, order, 0.0, 0.0)

    f = LMADF(order)
    y = synthesis!(f, base_excitation, c, hopsize)
    wavwrite(y, savepath; Fs=fs)
    println("Dumped to $savepath")
end

function test_mlsadf_synthesis(; order=25, savepath="test16k_mlsadf.wav")
    println("testing: mlsadf_synthesis")
    α = mcepalpha(fs) # automatic α selection
    mc = mcep(xw, order, α)

    f = MLSADF(order, α)
    y = synthesis!(f, base_excitation, mc, hopsize)
    wavwrite(y, savepath; Fs=fs)
    println("Dumped to $savepath")
end

function test_mglsadf_synthesis(; order=25, nstage=6,
                                savepath="test16k_mglsadf.wav")
    println("testing: mglsadf_synthesis")
    α = mcepalpha(fs)
    γ = -1./nstage
    mgc = mgcep(xw, order, α, γ)

    # waveform synthesis
    f = MGLSADF(order, α, nstage)
    y = synthesis!(f, base_excitation, mgc, hopsize)
    wavwrite(y, savepath; Fs=fs)
    println("Dumped to $savepath")
end

@time test_poledf_synthesis()
@time test_lmadf_synthesis()
@time test_mlsadf_synthesis()
@time test_mglsadf_synthesis()
