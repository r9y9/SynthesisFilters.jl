using MelGeneralizedCepstrums
using SynthesisFilters
using Base.Test

@static is_linux() ? include("sptk.jl") : nothing

type TestSynthesisFilter <: SynthesisFilter
end

let
    @test_throws ErrorException to_filtcoef(TestSynthesisFilter(), rand(10))
    @test_throws ErrorException filt!(TestSynthesisFilter(), 1.0, rand(10))
end

function test_synthesis_one_frame(f::SynthesisFilter, order)
    srand(98765)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test all(isfinite.(r))
end

function test_poledf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    l = estimate(LinearPredictionCoef(order), x)

    f = AllPoleDF(order)
    r = synthesis!(f, excite, l, hopsize)
    @test all(isfinite.(r))
end

function test_synthesis(def::SpectralParam, hopsize)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    order = param_order(def)

    if isa(def, LineSpectralPair)
        state = lpc2lsp(estimate(LinearPredictionCoef(order), x))
    elseif isa(def, PartialAutoCorrelation)
        state = lpc2par(estimate(LinearPredictionCoef(order), x))
    else
        state = estimate(def, x)
    end

    r = synthesis(excite, state, hopsize)
    @test all(isfinite.(r))
end

function test_poledf_exception()
    srand(98765)
    N = 512
    order = 25
    hopsize = 80
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    l = estimate(LinearPredictionCoef(order), x)
    f = AllPoleDF(order)
    synthesis!(f, excite, l, hopsize)

    for def in [LinearCepstrum(order),
                GeneralizedCepstrum(order, -0.1),
                MelCepstrum(order, 0.41),
                MelGeneralizedCepstrum(order, 0.41, -0.1)
                ]
        l = estimate(def, x)
        @test_throws ErrorException synthesis!(f, excite, l, hopsize)
    end
end

function test_ltcdf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    N = 1024
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    par = lpc2par(estimate(LinearPredictionCoef(order), x))

    f = AllPoleLatticeDF(order)
    r = synthesis!(f, excite, par, hopsize)
    @test all(isfinite.(r))
end

function test_lspdf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    N = 1024
    excite = rand(N)
    x = rand(order+1, div(N, hopsize))
    lsp = lpc2lsp(estimate(LinearPredictionCoef(order), x))

    f = LSPDF(order)
    r = synthesis!(f, excite, lsp, hopsize)
    @test all(isfinite.(r))
end

function test_lmadf_synthesis(order::Int, pade::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(LinearCepstrum(order), x)

    f = LMADF(order; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test all(isfinite.(r))
end

function test_lmadf_exception()
    srand(98765)
    N = 512
    order = 25
    hopsize = 80
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    mgc = estimate(LinearCepstrum(order), x)
    f = LMADF(order)
    synthesis!(f, excite, mgc, hopsize)

    for def in [
                MelCepstrum(order, 0.41),
                GeneralizedCepstrum(order, -0.01),
                MelGeneralizedCepstrum(order, 0.41, -0.01),
                AllPoleCepstrum(order)
                ]
        mgc = estimate(def, x)
        @test_throws ErrorException synthesis!(f, excite, mgc, hopsize)
    end

    ## pade approximation
    @test_throws ErrorException LMADF(order, pade=3)
    try LMADF(order; pade=4); catch @test false; end
    try LMADF(order; pade=5); catch @test false; end
    @test_throws ErrorException LMADF(order, pade=6)
end

function test_mlsadf_synthesis(order::Int, α::Float64, pade::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(MelCepstrum(order, α), x)

    f = MLSADF(order, α; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test all(isfinite.(r))
end

function test_mlsadf_exception()
    srand(98765)
    N = 512
    order = 25
    hopsize = 80
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    # should accept LinearCepstrum
    mgc = estimate(LinearCepstrum(order), x)
    f = MLSADF(order, 0.0)
    synthesis!(f, excite, mgc, hopsize)

    mgc = estimate(MelCepstrum(order, 0.41), x)
    f = MLSADF(order, 0.41)
    synthesis!(f, excite, mgc, hopsize)

    for def in [
                GeneralizedCepstrum(order, -0.01),
                MelGeneralizedCepstrum(order, 0.41, -0.01),
                AllPoleCepstrum(order)
                ]
        mgc = estimate(def, x)
        @test_throws ErrorException synthesis!(f, excite, mgc, hopsize)
    end

    ## pade approximation
    @test_throws ErrorException MLSADF(order, 0.41; pade=3)
    try MLSADF(order, 0.41; pade=4); catch @test false; end
    try MLSADF(order, 0.41; pade=5); catch @test false; end
    @test_throws ErrorException MLSADF(order, 0.41; pade=6)
end

function test_mglsadf_synthesis(order::Int, α::Float64, ns::Int, hopsize::Int)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(MelGeneralizedCepstrum(order, α, -1/ns), x)

    f = MGLSADF(order, α, ns)
    r = synthesis!(f, excite, mgc, hopsize)
    @test all(isfinite.(r))
end

function test_mglsadf_exceptions()
    N = 512
    order = 25
    hopsize = 80
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    ns = 10.5
    α = 0.41
    mgc = estimate(MelGeneralizedCepstrum(order, α, -1/ns), x)
    # should raise exception for non-integer ns
    @test_throws ErrorException synthesis(excite, mgc, hopsize)
end

function test_synthesis()
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(MelGeneralizedCepstrum(order, α, -1/ns), x)
end

### Synthesis with AllPoleDF ###

test_poledf_exception()

for order in 20:5:30
    println("poledf_synthesis_one_frame: testing with order=$order")
    f = AllPoleDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:30
    for hopsize in [80, 160]
        println("poledf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_poledf_synthesis(order, hopsize)
    end
end

### Synthesis with AllPoleLatticeDF ###

for order in 20:5:30
    println("ltcdf_synthesis_one_frame: testing with order=$order")
    f = AllPoleLatticeDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:30
    for hopsize in [80, 160]
        println("ltcdf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_ltcdf_synthesis(order, hopsize)
    end
end

### Synthesis with LSPDF ###

for order in 20:5:30
    println("lspdf_synthesis_one_frame: testing with order=$order")
    f = LSPDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:30
    for hopsize in [80, 160]
        println("lspdf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_lspdf_synthesis(order, hopsize)
    end
end

### Synthesis with LMADF ###

let
    f = LMADF(20)
    @test SynthesisFilters.allpass_alpha(f) == 0.0
    @test SynthesisFilters.glog_gamma(f) == 0.0
end

test_lmadf_exception()

for order in 20:5:30
    for pade in [4, 5]
        println("lamdf_synthesis_one_frame: testing with order=$order, pade=$pade")
        f = LMADF(order, pade=pade)
        test_synthesis_one_frame(f, order)
    end
end

for order in 20:5:30
    for pade in [4, 5]
        for hopsize in [80, 160]
            println("lmadf_synthesis: testing with order=$order, pade=$pade, hopsize=$hopsize")
            test_lmadf_synthesis(order, pade, hopsize)
        end
    end
end

### Synthesis with MLSADF ###

let
    f = MLSADF(20, 0.41)
    @test SynthesisFilters.allpass_alpha(f) == 0.41
    @test SynthesisFilters.glog_gamma(f) == 0.0
end

test_mlsadf_exception()

for order in 20:5:30
    for α in [0.35, 0.544]
        for pade in [4, 5]
            println("mlsadf_synthesis_one_frame: testing with order=$order, α=$α, pade=$pade")
            f = MLSADF(order, α, pade=pade)
            test_synthesis_one_frame(f, order)
        end
    end
end

for order in 20:5:30
    for α in [0.35, 0.544]
        for pade in [4, 5]
            for hopsize in [80, 160]
                println("mlsadf_synthesis: testing with order=$order, α=$α, pade=$pade, hopsize=$hopsize")
                test_mlsadf_synthesis(order, α, pade, hopsize)
            end
        end
    end
end

### Synthesis with MGLSADF ###

let
    nstage = 10
    γ = -1/nstage
    f = MGLSADF(20, 0.41, nstage)
    @test SynthesisFilters.allpass_alpha(f) == 0.41
    @test SynthesisFilters.glog_gamma(f) == γ
end

for order in 20:5:30
    for α in [0.35, 0.544]
        for ns in 2:10
            println("mglsadf_synthesis_one_frame: testing with order=$order, α=$α, nstage=$ns, γ=$(-1.0/ns)")
            f = MGLSADF(order, α, ns)
            test_synthesis_one_frame(f, order)
        end
    end
end

for order in 20:5:30
    for α in [0.35, 0.544]
        for ns in 2:10
            for hopsize in [80, 160]
                println("mglsadf_synthesis: testing with order=$order, α=$α, nstage=$ns, γ=$(-1.0/ns), hopsize=$hopsize")
                test_mglsadf_synthesis(order, α, ns, hopsize)
            end
        end
    end
end

test_mglsadf_exceptions()

### Test Syntehsis ###

for def in [
            LinearPredictionCoef(20),
            LineSpectralPair(20),
            PartialAutoCorrelation(20),
            LinearCepstrum(20),
            MelCepstrum(20, 0.41),
            GeneralizedCepstrum(20, -1/10),
            MelGeneralizedCepstrum(20, 0.41, -1/10)
            ]
    println("test_synthesis: testing with $(def)")
    test_synthesis(def, 80)
end
