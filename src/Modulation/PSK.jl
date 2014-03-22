module Modulation

export pskmod

function pskmod( data, M::Integer, encoding::String )
    m = [ 0 : M-1 ]
    if M == 4
        Φ = pi/M
    else
        Φ = 0.0
    end
    ideal = exp(2*pi*m*im/M + Φ*im)

    outputVec = Array(Complex128, length(data))

    for i in 1:length(data)
        outputVec[i] = ideal[ data[i] + 1 ]
    end
    return outputVec
end

function pskmod( M::Integer, length::Integer )
    data = rand( 0:M-1, length )
    pskmod( data, M, "" )
end

end # module Modulation