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

for (α, pade) = [(0.35, 5), (0.41, 5), (0.35, 4), (0.41, 4)]
    println("mlsadf: testing with α=$α, pade=$pade")
    test_mlsadf(α, pade)
end
