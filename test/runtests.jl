using MelGeneralizedCepstrums
using SynthesisFilters
using Base.Test

import SPTK

function test_poledf()
    srand(98765)
    x = rand(100)
    c = rand(21)
    order = length(c)-1

    f = AllPoleDF(order)

    delay = SPTK.poledf_delay(order)

    for i=1:length(x)
        y = SPTK.poledf(x[i], c, delay)
        ŷ = filter!(f, x[i], c)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_lmadf(pade::Int)
    srand(98765)
    x = rand(100)
    c = rand(21)
    order = length(c)-1

    f = LMADF(order, pade=pade)

    # setup for SPTK mlsadf
    delay = SPTK.lmadf_delay(order, pade)

    for i=1:length(x)
        y = SPTK.lmadf(x[i], c, pade, delay)
        ŷ = filter!(f, x[i], c)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_mlsadf(α::Float64, pade::Int)
    srand(98765)
    x = rand(100)
    mc = rand(21)
    order = length(mc)-1

    f = MLSADF(order, α, pade=pade)
    @test allpass_alpha(f) == α
    @test glog_gamma(f) == 0.0

    # setup for SPTK mlsadf
    delay = SPTK.mlsadf_delay(order, pade)

    for i=1:length(x)
        y = SPTK.mlsadf(x[i], mc, α, pade, delay)
        ŷ = filter!(f, x[i], mc)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_mglsadf(α::Float64=0.41, ns::Int=10)
    srand(98765)
    x = rand(100)
    mc = rand(21)
    order = length(mc)-1

    f = MGLSADF(order, α, ns)
    @test allpass_alpha(f) == α
    @test nstage(f) == ns
    @test glog_gamma(f) == -1.0/ns

    # setup for SPTK mlsadf
    delay = SPTK.mglsadf_delay(order, ns)

    for i=1:length(x)
        y = SPTK.mglsadf(x[i], mc, α, ns, delay)
        ŷ = filter!(f, x[i], mc)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_poledf_synthesis_one_frame(order::Int)
    srand(98765)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    f = AllPoleDF(order)
    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test any(!isnan(r))
end

function test_poledf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    T = 1024
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))
    l = MelLinearPredictionCoef(0.0, c, false)

    f = AllPoleDF(order)
    r = synthesis!(f, excite, l, hopsize)
    @test any(!isnan(r))
end

function test_poledf_exception()
    srand(98765)
    T = 1024
    order = 25
    hopsize = 80
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))

    l = MelLinearPredictionCoef(0.0, c, false)
    f = AllPoleDF(order)
    try
        synthesis!(f, excite, l, hopsize)
    catch
        @test false
    end
    l = MelLinearPredictionCoef(0.41, c, false)
    @test_throws ArgumentError synthesis!(f, excite, l, hopsize)
end

function test_lmadf_synthesis_one_frame(order::Int, pade::Int)
    srand(98765)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    f = LMADF(order; pade=pade)
    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test any(!isnan(r))
end

function test_lmadf_synthesis(order::Int, pade::Int, hopsize::Int)
    srand(98765)
    T = 1024
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))
    mgc = MelGeneralizedCepstrum(0.0, 0.0, c)

    f = LMADF(order; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

function test_lmadf_exception()
    srand(98765)
    T = 1024
    order = 25
    hopsize = 80
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))

    mgc = MelGeneralizedCepstrum(0.0, 0.0, c)
    f = LMADF(order)
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end
    mgc = MelGeneralizedCepstrum(0.41, 0.0, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
    mgc = MelGeneralizedCepstrum(0.41, -0.01, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
    mgc = MelGeneralizedCepstrum(0.0, -0.01, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
    mgc = MelGeneralizedCepstrum(0.0, -1.0, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
end

function test_mlsadf_synthesis_one_frame(order::Int, α::Float64, pade::Int)
    srand(98765)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    f = MLSADF(order, α; pade=pade)
    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test any(!isnan(r))
end

function test_mlsadf_synthesis(order::Int, α::Float64, pade::Int, hopsize::Int)
    srand(98765)
    T = 1024
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))
    mgc = MelGeneralizedCepstrum(α, 0.0, c)

    f = MLSADF(order, α; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

function test_mlsadf_exception()
    srand(98765)
    T = 1024
    order = 25
    hopsize = 80
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))

    mgc = MelGeneralizedCepstrum(0.0, 0.0, c)
    f = MLSADF(order, 0.0)
    synthesis!(f, excite, mgc, hopsize)
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end
    mgc = MelGeneralizedCepstrum(0.41, 0.0, c)
    f = MLSADF(order, 0.41)
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end
    mgc = MelGeneralizedCepstrum(0.41, -0.01, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
    mgc = MelGeneralizedCepstrum(0.0, -0.01, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
    mgc = MelGeneralizedCepstrum(0.0, -1.0, c)
    @test_throws ArgumentError synthesis!(f, excite, mgc, hopsize)
end

function test_mglsadf_synthesis_one_frame(order::Int, α::Float64, ns::Int)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    f = MGLSADF(order, α, ns)
    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test any(!isnan(r))
end

function test_mglsadf_synthesis(order::Int, α::Float64, ns::Int, hopsize::Int)
    T = 1024
    excite = rand(T)
    c = rand(order+1, div(T, hopsize))
    mgc = MelGeneralizedCepstrum(α, -1/ns, c)

    f = MGLSADF(order, α, ns)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

## Consistency tests

println("poledf: testing")
test_poledf()

for pade in [4, 5]
    println("lmadf: testing with pade=$pade")
    test_lmadf(pade)
end

for α in [0.35, 0.41, 0.544]
    for pade in [4, 5]
        println("mlsadf: testing with α=$α, pade=$pade")
        test_mlsadf(α, pade)
    end
end

for α in [0.35, 0.41, 0.544]
    for ns in 1:15
        println("mglsadf: testing with α=$α, nstage=$ns, γ=$(-1.0/ns)")
        test_mglsadf(α, ns)
    end
end

## Synthesis with AllPoleDF

test_poledf_exception()

for order in 20:5:40
    println("poledf_synthesis_one_frame: testing with order=$order")
    test_poledf_synthesis_one_frame(order)
end

for order in 20:5:40
    for hopsize in [80, 160]
        println("poledf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_poledf_synthesis(order, hopsize)
    end
end

## Synthesis with LMADF

test_lmadf_exception()

for order in 20:5:40
    for pade in [4, 5]
        println("lamdf_synthesis_one_frame: testing with order=$order, pade=$pade")
        test_lmadf_synthesis_one_frame(order, pade)
    end
end

for order in 20:5:40
    for pade in [4, 5]
        for hopsize in [80, 160]
            println("lmadf_synthesis: testing with order=$order, pade=$pade, hopsize=$hopsize")
            test_lmadf_synthesis(order, pade, hopsize)
        end
    end
end

## Synthesis with MLSADF

test_mlsadf_exception()

for order in 20:5:40
    for α in [0.35, 0.41, 0.544]
        for pade in [4, 5]
            println("mlsadf_synthesis_one_frame: testing with order=$order, α=$α, pade=$pade")
            test_mlsadf_synthesis_one_frame(order, α, pade)
        end
    end
end

for order in 20:5:40
    for α in [0.35, 0.41, 0.544]
        for pade in [4, 5]
            for hopsize in [80, 160]
                println("mlsadf_synthesis: testing with order=$order, α=$α, pade=$pade, hopsize=$hopsize")
                test_mlsadf_synthesis(order, α, pade, hopsize)
            end
        end
    end
end

## Synthesis with MGLSADF

for order in 20:5:40
    for α in [0.35, 0.41, 0.544]
        for ns in 2:10
            println("mglsadf_synthesis_one_frame: testing with order=$order, α=$α, nstage=$ns, γ=$(-1.0/ns)")
            test_mglsadf_synthesis_one_frame(order, α, ns)
        end
    end
end

for order in 20:5:40
    for α in [0.35, 0.41, 0.544]
        for ns in 2:10
            for hopsize in [80, 160]
                println("mglsadf_synthesis: testing with order=$order, α=$α, nstage=$ns, γ=$(-1.0/ns), hopsize=$hopsize")
                test_mglsadf_synthesis(order, α, ns, hopsize)
            end
        end
    end
end
