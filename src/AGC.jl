################################################################################
#                                ____ ____ ____                                #
#                                |__| | __ |                                   #
#                                |  | |__] |___                                #
################################################################################                                     

type AGC
    gain
    bandwidth
    alpha
    y2
    initialized
end


function AGC( gain = 1.0, bandwidth = 0.01, alpha = bandwidth  )
    AGC( gain, bandwidth, alpha, 1.0, false )
end

function exec( agc::AGC, x::Number )
    y = x * agc.gain
    Ey = abs2(y)
    agc.y2 = (1-agc.alpha)*agc.y2 + Ey
    agc.gain *= exp( -0.5*agc.alpha*log(agc.y2) )
    return y
end

function exec{T}( agc::AGC, x::Vector{T} )
    T[ exec( agc, xx) for xx in x ]
end

using PyPlot

agc = AGC()
x = [cos(2*pi*0.05*t) for t in 0:10000]
y = exec( agc, x )

clf()
subplot(2,1,1)
plot(x)
subplot(2,1,2)
plot(y)
