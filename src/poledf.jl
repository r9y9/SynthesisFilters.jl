# All pole digital filter for speech synthesis from LPC

poledf_delay(order::Int) = zeros(order)

type AllPoleDF <: MelGeneralizedCepstrumSynthesisFilter
    delay::Vector{Float64}

    function AllPoleDF(order::Int)
        new(poledf_delay(order))
    end
end

function filter!(f::AllPoleDF, x::Float64, a::Vector{Float64})
    d = delay(f)
    order = length(a) - 1

    y = x
    @inbounds for m=order:-1:2
        y -= a[m+1] * d[m]
        d[m] = d[m-1]
    end

    y -= a[2] * d[1]
    d[1] = y

    return y
end
