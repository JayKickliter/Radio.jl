################################################################################
#                                ____ ____ ____                                #
#                                |__| | __ |                                   #
#                                |  | |__] |___                                #
################################################################################                                     

type AGC
    gain
    bandwidth
    alpha
    y2_prime
    initialized
end


function AGC( gain = 1.0, bandwidth = 0.01, alpha = bandwidth  )
    AGC( gain, bandwidth, alpha, 1.0, false )
end

function exec( agc::AGC, x::Number )
    y            = x * agc.gain
    y2           = abs2( y )
    agc.y2_prime = (1 - agc.alpha) * agc.y2_prime + agc.alpha*y2
    agc.gain *= exp( -0.5 * agc.alpha * log(agc.y2_prime))
    
    return y
end

function exec{T}( agc::AGC, x::Vector{T} )
    if ! agc.initialized
        for i in 1:min(1000, length(x))
           exec( agc, x[i]) 
        end
    end
    T[ exec( agc, xx) for xx in x ]
end

#=
include( "Modem/QAM.jl" )

agc        = AGC(1e-1)
M          = 16
modem      = QAMModem( M )
data       = rand( 0:M-1, 1000 )
symbols    = modulate( modem, data )
agcSymbols = exec( agc, symbols )
demodData  = demodulate( modem, agcSymbols )

Test.@test_approx_eq data demodData

using PyPlot
clf()
subplot( 2, 1, 1 )
plot( real(symbols), "b" )
plot( imag(symbols), "g" )

subplot( 2, 1, 2 )
plot( real(agcSymbols), "b" )
plot( imag(agcSymbols), "g" )
=#