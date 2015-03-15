include( "../coding.jl" )

abstract Modem

################################################################################
#                                ___  ____ _  _                                #
#                                |__] [__  |_/                                 #
#                                |    ___] | \_                                #
################################################################################                                     

type PSKModem
    M::Int                             # Number of symbols in the constellation
    m::Int                             # Bits per symbol
    constellation::Vector{Complex}     # Symbol constellation
    grayConstellation::Vector{Complex} # Symbol constellation with gray code
    bitsMap::Vector{Int}               # bitsMap[i] maps to constellation[i]
end

function PSKModem( M::Integer )
    ispow2( M ) || error( "M must be a power of 2" )
    Δϕ                = 2π/M
    bitsMap           = [ encode( Gray, i ) for i in 0:M-1 ]
    constellation     = [ exp(Δϕ*im*i) for i in 0:M-1 ]
    grayConstellation = [ exp(Δϕ*im*decode( Gray, i)) for i in 0:M-1 ]    
    m                 = log2( M )
    PSKModem( M, m, constellation, grayConstellation, bitsMap )
end




################################################################################
#                                _  _ ____ ___                                 #
#                                |\/| |  | |  \                                #
#                                |  | |__| |__/                                #
################################################################################                                     

function modulate( modem::PSKModem, bits::Integer )
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

function demodulate( modem::PSKModem, symbol::Complex )
    ϕ = angle( symbol )
    ϕ = ϕ < 0 ? ϕ += 2π : ϕ
    
    index = int( ϕ*modem.M / 2π ) + 1
    modem.bitsMap[index]
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

function PyPlot.plot( modem::PSKModem )
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

modem   = PSKModem( 4 )

plot( modem )

modData   = 0:modem.M-1 #rand(0:8-1, 10)
symbols   = modulate( modem, modData)
demodData = demodulate( modem, symbols )

Test.@test_approx_eq modData demodData

=# 
