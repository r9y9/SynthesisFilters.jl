using SynthesisFilters
using Base.Test

import SPTK

function test_mlsadf(α::Float64, pade::Int)
    srand(98765)
    x = rand(100)
    mc = rand(21)
    order = length(mc)-1

    f = MLSADF(order, α, pade=pade)
    @test alpha(f) == α
    @test gamma(f) == 0.0

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
    @test alpha(f) == α
    @test nstage(f) == ns
    @test gamma(f) == -1.0/ns

    # setup for SPTK mlsadf
    delay = SPTK.mglsadf_delay(order, ns)

    for i=1:length(x)
        y = SPTK.mglsadf(x[i], mc, α, ns, delay)
        ŷ = filter!(f, x[i], mc)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_mlsadf_synthesis_one_frame(order::Int, α::Float64, pade::Int)
    excite = rand(1024)
    previous_mgc = rand(order+1)
    current_mgc = rand(order+1)

    f = MLSADF(order, α; pade=pade)
    r = synthesis_one_frame!(f, excite, previous_mgc, current_mgc)
    @test any(!isnan(r))
end

function test_mlsadf_synthesis(order::Int, α::Float64, pade::Int, hopsize::Int)
    T = 1024
    excite = rand(T)
    mgc = rand(order+1, div(T, hopsize))

    f = MLSADF(order, α; pade=pade)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
end

function test_mglsadf_synthesis_one_frame(order::Int, α::Float64, ns::Int)
    excite = rand(1024)
    previous_mgc = rand(order+1)
    current_mgc = rand(order+1)

    f = MGLSADF(order, α, ns)
    r = synthesis_one_frame!(f, excite, previous_mgc, current_mgc)
    @test any(!isnan(r))
end

function test_mglsadf_synthesis(order::Int, α::Float64, ns::Int, hopsize::Int)
    T = 1024
    excite = rand(T)
    mgc = rand(order+1, div(T, hopsize))

    f = MGLSADF(order, α, ns)
    r = synthesis!(f, excite, mgc, hopsize)
    @test any(!isnan(r))
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

## Synthesis with MLSADF

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
