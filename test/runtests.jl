using SynthesisFilters
using Base.Test

import SPTK

function test_mlsadf(α::Float64=0.41, pade::Int=5)
    srand(98765)
    x = rand(100)
    mc = rand(21)
    order = length(mc)-1

    f = MLSADF(order, α, pade=pade)

    # setup for SPTK mlsadf
    delay = SPTK.mlsadf_delay(order, pade)

    for i=1:length(x)
        y = SPTK.mlsadf(x[i], mc, α, pade, delay)
        ŷ = filter!(f, x[i], mc)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_mglsadf(α::Float64=0.41, nstage::Int=10)
    srand(98765)
    x = rand(100)
    mc = rand(21)
    order = length(mc)-1

    f = MGLSADF(order, α, nstage)

    # setup for SPTK mlsadf
    delay = SPTK.mglsadf_delay(order, nstage)

    for i=1:length(x)
        y = SPTK.mglsadf(x[i], mc, α, nstage, delay)
        ŷ = filter!(f, x[i], mc)
        @test !isnan(ŷ)
        @test_approx_eq y ŷ
    end
end

function test_gnorm(γ::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    g = SPTK.gnorm(mc, γ)
    ĝ = gnorm(mc, γ)
    @test_approx_eq g ĝ
end

function test_mc2b(α::Float64)
    srand(98765)
    x = rand(100)
    mc = rand(21)

    g = SPTK.mc2b(mc, α)
    ĝ = mc2b(mc, α)
    @test_approx_eq g ĝ
end

for α in [0.35, 0.41, 0.544]
    for pade in [4, 5]
        println("mlsadf: testing with α=$α, pade=$pade")
        test_mlsadf(α, pade)
    end
end

for α in [0.35, 0.41, 0.544]
    for nstage in 1:15
        println("mglsadf: testing with α=$α, nstage=$nstage, γ=$(-1.0/nstage)")
        test_mglsadf(α, nstage)
    end
end

for nstage in 1:15
    γ=-1.0/nstage
    println("gnorm: testing with γ=$γ")
    test_gnorm(γ)
end

for α in [0.35, 0.41, 0.544]
    println("mc2b: testing with α=$α")
    test_mc2b(α)
end
