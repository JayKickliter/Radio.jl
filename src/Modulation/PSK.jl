#==============================================================================#
#                         Phase Shift Keying Modulation                        #
#==============================================================================#

function pskmod( data, M::Integer, encoding::String, SPS::Integer = 1, ISIFilter::Vector = [] )
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
    
    if SPS == 1
        return outputVec
    end

    if ISIFilter == []        
        ISIFilter = rcos( 0.3, 10, SPS )
    end    
    
    outputVec = upsample( outputVec, SPS )
    outputVec = filt( complex128(ISIFilter), complex128(1.0), outputVec )
end

function pskmod( symbols::Integer, M::Integer, SPS::Integer = 1, ISIFilter::Vector = [] )
    data = rand( 0:M-1, symbols )
    pskmod( data, M, "", SPS, ISIFilter )
end