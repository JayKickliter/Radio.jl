include( "../coding.jl" )

abstract Modem

################################################################################
#                                ____ ____ _  _                                #
#                                |  | |__| |\/|                                #
#                                |_\| |  | |  |                                #
################################################################################                                     

type QAMModem
    M::Int                             # Number of symbols in the constellation
    m::Int                             # Bits per symbol
    α::FloatingPoint                   # Scaling factor
    constellation::Vector{Complex}     # Symbol constellation
    grayConstellation::Vector{Complex} # Symbol constellation with gray code
    bitsMap::Vector{Int}               # bitsMap[i] maps to constellation[i]
end

function QAMModem( M::Integer )
    m                 = sqrt( M )
    isinteger( m )    || error( "sqrt(M) must be an integer value" )
    m = int( m )
    m2 = int(m/2)
    m22 = m2/2
    constellation = Array( Complex64, M )
    z = [-3, -1, 1, 3]
    
    # http://www.gaussianwaves.com/2012/10/constructing-a-rectangular-constellation-for-16-qam/
    
    for index in 0:M-1
        inPhase = index >> m2
        quadrature = index & (m2^2-1)
        println( inPhase, " ", quadrature )
        println( (z[decode(Gray, inPhase) + 1], z[decode( Gray, quadrature ) + 1]) )
        constellation[index + 1] = Complex( decode(Gray, inPhase)*2-m+1, decode(Gray, quadrature)*2-m+1)
    end
    constellation
    
    
    bitsMap           = [ encode( Gray, i ) for i in 0:M-1 ]
    constellation     = [  for i in 0:M-1 ]
    grayConstellation = [ exp(Δϕ*im*decode( Gray, i)) for i in 0:M-1 ]    
    m                 = log2( M )
    
    QAMModem( M, m, constellation, grayConstellation, bitsMap )
end




################################################################################
#                           ____ _    ___  _  _ ____                           #
#                           |__| |    |__] |__| |__|                           #
#                           |  | |___ |    |  | |  |                           #
################################################################################                                     

function α( M::Integer )
    1/sqrt( (M-1) * 2/3 )
end




################################################################################
#                                _  _ ____ ___                                 #
#                                |\/| |  | |  \                                #
#                                |  | |__| |__/                                #
################################################################################                                     

function modulate( modem::QAMModem, bits::Integer )
    modem.grayConstellation[ bits+1 ]
end


function modulate( modem, data::AbstractVector )
    [ modulate( modem, datum ) for datum in data ]
end




################################################################################
#                           ___  ____ _  _ ____ ___                            #
#                           |  \ |___ |\/| |  | |  \                           #
#                           |__/ |___ |  | |__| |__/                           #
################################################################################     

function demodulate( modem::QAMModem, symbol::Complex )
end


function demodulate( modem, symbols::AbstractVector{Complex} )
    Int[ demodulate( modem, symbol ) for symbol in symbols ]
end




################################################################################
#                      ____ _  _ ____ _  _ ___  _    ____                      #
#                      |___  \/  |__| |\/| |__] |    |___                      #
#                      |___ _/\_ |  | |  | |    |___ |___                      #
################################################################################     

#=

using PyPlot

function PyPlot.plot( modem::QAMModem )
    fig = figure()
    scatter( real(modem.constellation), imag(modem.constellation) )
    plot( real(modem.grayConstellation), imag(modem.grayConstellation) )
    
    for i in 0:modem.M-1
        symbol = modem.constellation[ i + 1 ]
        x      = real( symbol )
        y      = imag( symbol )
        bits   = modem.bitsMap[i+1]

        annotate( bin(bits, modem.m ), (x, y+0.05), ha="center" )
    end
    
    setp( axes(), aspect = 1  )
    xlabel( "In Phase" )
    ylabel( "Quadrature" )
    
    return fig
end

modem   = QAMModem( 4 )

plot( modem )

modData   = 0:modem.M-1 #rand(0:8-1, 10)
symbols   = modulate( modem, modData)
demodData = demodulate( modem, symbols )

Test.@test_approx_eq modData demodData

=# 
