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
    # constellation = [ exp(Δ𝜙*im*decode( Gray, i)) for i in 0: symbols-1 ]
    constellation = [ exp(Δ𝜙*im*i) for i in 0: symbols-1 ]
    
    return PSKModem( modOrder, constellation )
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

function modulate( modem::PSKModem, datum::Integer )
    modem.constellation[decode(Gray, datum)+1]
end

function modulate( modem, data::AbstractVector )
    [ modulate( modem, datum ) for datum in data ]
end

# https://www.rocq.inria.fr/scicos/ScicosModNum/modnum_web/src/modnum_421/interf/scicos/help/eng/htm/DEMODPSK_f.htm
function demodulate( modem::PSKModem, symbol::Complex )
    θ = angle( symbol )
    θ = θ < 0 ? θ += 2π : θ
    
    bits = int( ( θ*modem.modOrder^2 )/2π )
    encode( Gray, bits )
end

function demodulate( modem, symbols::AbstractVector{Complex} )
    [ demodulate( modem, symbol ) for symbol in symbols ]
end




#=

modem   = PSKModem( 4 )
data    = [0:2^modem.modOrder-1]
symbols = modulate( modem, data)
plot(symbols)
demodulate( modem, symbols )

=#
