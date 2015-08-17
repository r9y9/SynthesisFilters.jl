using MelGeneralizedCepstrums
using SynthesisFilters
using Base.Test

@unix_only include("sptk.jl")

type TestSynthesisFilter <: SynthesisFilter
end
type TestMGCSF <: MelGeneralizedCepstrumSynthesisFilter
end
type TestMLPSF <: LinearPredictionVariantsSynthesisFilter
end

let
    @test_throws Exception to_filtercoef(TestSynthesisFilter(), rand(10))
    f = TestMGCSF()
    @test_throws Exception SynthesisFilters.allpass_alpha(f)
    @test_throws Exception SynthesisFilters.glog_gamma(f)
    f = TestMLPSF()
    @test_throws Exception SynthesisFilters.allpass_alpha(f)
end

function test_synthesis_one_frame(f::SynthesisFilter, order)
    srand(98765)
    excite = rand(1024)
    previous_b = rand(order+1)
    current_b = rand(order+1)

    r = synthesis_one_frame!(f, excite, previous_b, current_b)
    @test any(!isnan(r))
end

function test_poledf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    l = estimate(LinearPredictionCoef(order), x)

    f = AllPoleDF(order)
    r = synthesis!(f, excite, l, hopsize)
    @test any(!isnan(r))
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
    try
        synthesis!(f, excite, l, hopsize)
    catch
        @test false
    end

    for def in [LinearCepstrum(order),
                GeneralizedCepstrum(order, -0.1),
                MelCepstrum(order, 0.41),
                MelGeneralizedCepstrum(order, 0.41, -0.1)
                ]
        l = estimate(def, x)
        @test_throws Exception synthesis!(f, excite, l, hopsize)
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
    @test any(!isnan(r))
end

function test_lspdf_synthesis(order::Int, hopsize::Int)
    srand(98765)
    N = 1024
    excite = rand(N)
    x = rand(order+1, div(N, hopsize))
    lsp = lpc2lsp(estimate(LinearPredictionCoef(order), x))

    f = LSPDF(order)
    r = synthesis!(f, excite, lsp, hopsize)
    @test any(!isnan(r))
end

function test_lmadf_synthesis(order::Int, pade::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(LinearCepstrum(order), x)

    f = LMADF(order; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
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
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end

    for def in [
                MelCepstrum(order, 0.41),
                GeneralizedCepstrum(order, -0.01),
                MelGeneralizedCepstrum(order, 0.41, -0.01),
                AllPoleCepstrum(order)
                ]
        mgc = estimate(def, x)
        @test_throws Exception synthesis!(f, excite, mgc, hopsize)
    end

    ## pade approximation
    @test_throws Exception LMADF(order, pade=3)
    try LMADF(order; pade=4); catch @test false; end
    try LMADF(order; pade=5); catch @test false; end
    @test_throws Exception LMADF(order, pade=6)
end

function test_mlsadf_synthesis(order::Int, α::Float64, pade::Int, hopsize::Int)
    srand(98765)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(MelCepstrum(order, α), x)

    f = MLSADF(order, α; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

function test_mlsadf_exception()
    srand(98765)
    N = 512
    order = 25
    hopsize = 80
    excite = rand(N)
    x = rand(N, div(N, hopsize))

    mgc = estimate(LinearCepstrum(order), x)
    f = MLSADF(order, 0.0)
    synthesis!(f, excite, mgc, hopsize)
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end
    mgc = estimate(MelCepstrum(order, 0.41), x)
    f = MLSADF(order, 0.41)
    try
        synthesis!(f, excite, mgc, hopsize)
    catch
        @test false
    end


    for def in [
                GeneralizedCepstrum(order, -0.01),
                MelGeneralizedCepstrum(order, 0.41, -0.01),
                AllPoleCepstrum(order)
                ]
        mgc = estimate(def, x)
        @test_throws Exception synthesis!(f, excite, mgc, hopsize)
    end

    ## pade approximation
    @test_throws Exception MLSADF(order, 0.41; pade=3)
    try MLSADF(order, 0.41; pade=4); catch @test false; end
    try MLSADF(order, 0.41; pade=5); catch @test false; end
    @test_throws Exception MLSADF(order, 0.41; pade=6)
end

function test_mglsadf_synthesis(order::Int, α::Float64, ns::Int, hopsize::Int)
    N = 512
    excite = rand(N)
    x = rand(N, div(N, hopsize))
    mgc = estimate(MelGeneralizedCepstrum(order, α, -1/ns), x)

    f = MGLSADF(order, α, ns)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

### Synthesis with AllPoleDF ###

test_poledf_exception()

for order in 20:5:40
    println("poledf_synthesis_one_frame: testing with order=$order")
    f = AllPoleDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:40
    for hopsize in [80, 160]
        println("poledf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_poledf_synthesis(order, hopsize)
    end
end

### Synthesis with AllPoleLatticeDF ###

for order in 20:5:40
    println("ltcdf_synthesis_one_frame: testing with order=$order")
    f = AllPoleLatticeDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:40
    for hopsize in [80, 160]
        println("ltcdf_synthesis: testing with order=$order, hopsize=$hopsize")
        test_ltcdf_synthesis(order, hopsize)
    end
end

### Synthesis with LSPDF ###

for order in 20:5:40
    println("lspdf_synthesis_one_frame: testing with order=$order")
    f = LSPDF(order)
    test_synthesis_one_frame(f, order)
end

for order in 20:5:40
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
