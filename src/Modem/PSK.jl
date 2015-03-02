abstract Modem
abstract CodingScheme



################################################################################
#               ____ _  _ ____ ____ ___  _ _  _ ____                           #
#               |___ |\ | |    |  | |  \ | |\ | | __                           #
#               |___ | \| |___ |__| |__/ | | \| |__]                           #
################################################################################                                     
                                     
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
decode{T<:CodingScheme}( ::Type{T}, N::AbstractVector ) = [ decode( T, n) for n in N ]




################################################################################
#                                ___  ____ _  _                                #
#                                |__] [__  |_/                                 #
#                                |    ___] | \_                                #
################################################################################                                     

type PSKModem
    M::Integer                     # Number of symbols in the constellation
    constellation::Vector{Complex} # Symbol constellation
end

function PSKModem( M::Integer )
    ispow2( M ) || error( "M must be a power of 2" )
    Δ𝜙            = 2π/M
    constellation = [ exp(Δ𝜙*im*i) for i in 0: M-1 ]
    PSKModem( M, constellation )
end




################################################################################
#                                _  _ ____ ___                                 #
#                                |\/| |  | |  \                                #
#                                |  | |__| |__/                                #
################################################################################                                     

function modulate( modem::PSKModem, bits::Integer )
    modem.constellation[decode( Gray, bits )+1]
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
    
    bits = int( ϕ*modem.M / 2π )
    encode( Gray, bits )
end


function demodulate( modem, symbols::AbstractVector{Complex} )
    [ demodulate( modem, symbol ) for symbol in symbols ]
end

# function Winston.plot( modem::PSKModem )
#     p             = scatter( modem.constellation )
#     M             = modem.M
#     constellation = modem.constellation
#     for i in 0:M-1
#         symbol = constellation[ i + 1 ]
#         x      = real( symbol )
#         y      = imag( symbol )
#         bits   = decode( Gray, i )
#
#         println( x, " ", y, " ", bin(bits) )
#
#         add( p, PlotLabel(x*0.05, y*0.05, bin(bits)) )
#     end
#
#     return p
# end



#=
using Winston

modem   = PSKModem( 8 )
data    = [0:modem.M-1]
symbols = modulate( modem, data)
plot( symbols, "o-", xlabel="I", ylabel="Q" )
demodulate( modem, symbols )



=# 
