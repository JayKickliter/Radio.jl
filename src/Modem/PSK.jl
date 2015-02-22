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
    modOrder::Integer                   # Modulation order, or bits per symbol. The constellation has modOrder^2 symbols
    constellation::Vector{Complex} # ideal symbol constellation
end

function PSKModem( symbols::Integer )
    ispow2( symbols ) || error( "symbols must be a power of 2" )
    
    modOrder      = log2( symbols )
    Δ𝜙            = 2π/symbols
    constellation = [ exp(Δ𝜙*im*decode( Gray, i)) for i in 0: symbols-1 ]
    
    return PSKModem( modOrder, constellation )
end




function modulate( modem::PSKModem, datum::Integer )
    modem.constellation[datum+1]
end

function modulate( modem, data::AbstractVector )
    [ modulate( modem, datum ) for datum in data ]
end

function demodulate( modem::PSKModem, datum::Integer )
    
end

function demodulate( modem, data::AbstractVector{Integer} )
    [ demodulate( modem, datum ) for datum in data ]
end


#=

modem   = PSKModem( 16 )
symbols = modulate( modem, [0:2^modem.modOrder-1])
plot(symbols)

=#

