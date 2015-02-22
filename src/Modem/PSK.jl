#==============================================================================#
#                         Phase Shift Keying Modulation                        #
#==============================================================================#
# data  = Contains the integer data you want modulated.
#           It is up to you to ensure that 
# 
# 
# 
# 

function pskmod( data, modOrder::Integer, encoding::String = "", samplesPerSymbol::Integer = 1, ISIFilter::Vector = [] )
    m = [ 0 : modOrder-1 ]
    if modOrder == 4
        Φ = pi/modOrder
    else
        Φ = 0.0
    end
    ideal = exp(2*pi*m*im/modOrder .+ Φ*im)

    outputVec = Array(Complex128, length(data))

    for i in 1:length(data)
        outputVec[i] = ideal[ data[i] + 1 ]
    end
    
    if samplesPerSymbol == 1
        return outputVec
    end

    if ISIFilter == []        
        ISIFilter = rcos( 0.3, 10, samplesPerSymbol )
    end    
    
    outputVec = interpolate( outputVec, samplesPerSymbol, ISIFilter )
end

function pskmod( symbols::Integer, modOrder::Integer, samplesPerSymbol::Integer = 1, ISIFilter::Vector = [] )
    data = rand( 0:modOrder-1, symbols )
    pskmod( data, modOrder, "", samplesPerSymbol, ISIFilter )
end