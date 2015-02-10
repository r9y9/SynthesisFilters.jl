using MelGeneralizedCepstrums
using SynthesisFilters
using Base.Test

@unix_only include("sptk.jl")

type TestSynthesisFilter <: SynthesisFilter
end
type TestMGCSF <: MelGeneralizedCepstrumSynthesisFilter
end
type TestMLPSF <: MelLinearPredictionSynthesisFilter
end

let
    @test_throws Exception to_filtercoef(TestSynthesisFilter(), rand(10))
    f = TestMGCSF()
    @test_throws Exception allpass_alpha(f)
    @test_throws Exception glog_gamma(f)
    f = TestMGCSF()
    @test_throws Exception allpass_alpha(f)
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
