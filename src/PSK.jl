################################################################################
#                                ___  ____ _  _                                #
#                                |__] [__  |_/                                 #
#                                |    ___] | \_                                #
################################################################################                                     

type PSK <: Modulation
    M::Int                    # Number of symbols in the constellation
    m::Int                    # Bits per symbol
    constellation::Vector     # Symbol constellation
    grayConstellation::Vector # Symbol constellation with gray code
    bitsMap::Vector           # bitsMap[i] maps to constellation[i]
end

function PSK( M::Integer )
    ispow2( M ) || error( "M must be a power of 2" )
    Δϕ                = 2π/M
    bitsMap           = [ encode( Gray, i ) for i in 0:M-1 ]
    # constellation     = Complex128[ exp(Δϕ*im*i) for i in 0:M-1 ]
    # grayConstellation = Complex128[ exp(Δϕ*im*decode( Gray, i)) for i in 0:M-1 ]
    constellation     = Array( Complex128, M )
    grayConstellation = Array( Complex128, M )   
    for i in 0:M-1
       constellation[i+1]     = exp(Δϕ*im*i)
       grayConstellation[i+1] = exp(Δϕ*im*decode( Gray, i))
    end
    m                 = log2( M )
    PSK( M, m, constellation, grayConstellation, bitsMap )
end




################################################################################
#                                _  _ ____ ___                                 #
#                                |\/| |  | |  \                                #
#                                |  | |__| |__/                                #
################################################################################                                     

function modulate( psk::PSK, bits::Integer )
    psk.grayConstellation[ bits+1 ]
end


function modulate( psk::PSK, data::AbstractVector )
    T = eltype(psk.grayConstellation)
    T[ modulate( psk, datum ) for datum in data ]
end




################################################################################
#                           ___  ____ _  _ ____ ___                            #
#                           |  \ |___ |\/| |  | |  \                           #
#                           |__/ |___ |  | |__| |__/                           #
################################################################################     

function demodulate( psk::PSK, symbol::Complex )
    ϕ = angle( symbol )
    ϕ = ϕ < 0 ? ϕ += 2π : ϕ
    
    index = int( ϕ*psk.M / 2π ) + 1
    index = mod1( index, psk.M )
    psk.bitsMap[index]
end


function demodulate( psk::PSK, symbols::AbstractVector )
    Int[ demodulate( psk, symbol ) for symbol in symbols ]
end




################################################################################
#                      ____ _  _ ____ _  _ ___  _    ____                      #
#                      |___  \/  |__| |\/| |__] |    |___                      #
#                      |___ _/\_ |  | |  | |    |___ |___                      #
################################################################################     

#=

using PyPlot

function PyPlot.plot( psk::PSK )
    fig = figure()
    scatter( real(psk.constellation), imag(psk.constellation) )
    plot( real(psk.grayConstellation), imag(psk.grayConstellation) )
    
    for i in 0:psk.M-1
        symbol = psk.constellation[ i + 1 ]
        x      = real( symbol )
        y      = imag( symbol )
        bits   = psk.bitsMap[i+1]

        annotate( bin(bits, psk.m ), (x, y+0.05), ha="center" )
    end
    
    setp( axes(), aspect = 1  )
    xlabel( "In Phase" )
    ylabel( "Quadrature" )
    
    return fig
end

psk   = PSK( 4 )

plot( psk )

modData   = 0:psk.M-1 #rand(0:8-1, 10)
symbols   = modulate( psk, modData)
demodData = demodulate( psk, symbols )

Test.@test_approx_eq modData demodData

=# 
