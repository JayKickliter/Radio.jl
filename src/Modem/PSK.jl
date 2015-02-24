abstract Modem
abstract CodingScheme


type Gray <: CodingScheme end


function encode( ::Type{Gray}, n::Integer )
    n $ (n >> 1)
end


function decode( ::Type{Gray}, n::Integer )
    p = n
     while (n >>= 1) != 0
         p $= n
     end
     return p
end

encode{T<:CodingScheme}( ::Type{T}, N::AbstractVector ) = [ encode( T, n) for n in N ]



type PSKModem
    M::Integer                     # Modulation order, or bits per symbol. The constellation has M^2 symbols
    constellation::Vector{Complex} # ideal symbol constellation
end

function PSKModem( M::Integer )
    ispow2( M ) || error( "M must be a power of 2" )
    Δ𝜙            = 2π/M
    constellation = [ exp(Δ𝜙*im*decode( Gray, i)) for i in 0: M-1 ]
    PSKModem( M, constellation )
end


function symbol2index( symbol::Complex, constellationSize::Integer )
    θ = angle( symbol )

    if θ < 0
        θ += 2*pi
    end
    
    α = (constellationSize)/(2*pi)
    
    index = int(α * θ) + 1
    
    index = index > constellationSize ? 0 : index
    return index
end

function modulate( modem::PSKModem, bits::Integer )
    modem.constellation[bits+1]
end

function modulate( modem, data::AbstractVector )
    [ modulate( modem, datum ) for datum in data ]
end


function demodulate( modem::PSKModem, symbol::Complex )
    θ = angle( symbol )
    θ = θ < 0 ? θ += 2π : θ
    
    bits = int( θ*modem.M / 2π )
    encode( Gray, bits )
end

function demodulate( modem, symbols::AbstractVector{Complex} )
    [ demodulate( modem, symbol ) for symbol in symbols ]
end




#=

modem   = PSKModem( 4 )
data    = [0:modem.M-1]
symbols = modulate( modem, data)
plot(symbols)
demodulate( modem, symbols )

=# 
