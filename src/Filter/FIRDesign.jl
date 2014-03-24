#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# M      = Filter order; Returns M+1 samples
# F_t    = Normalized transition frequency. CutoffFreq/SampleRate
# window = Window function to apply to ideal truncated sinc response
function firdes( M::Integer, F_t::FloatingPoint, windowFunction::Function )
    if F_t < 0 && F_t > 0.5
        error("F_t must be > 0.0 and < 0.5")
    end
    
    [ 2*F_t*sinc(2*F_t*(n-M/2)) * windowFunction(n+1, M+1) for n = 0:M ]
end

#==============================================================================#
#                                Kaiser Order                                  #
#==============================================================================#
# δ  = ripple in passband and stopband
# Δω = transition width
# 
# returns:
#   M = approxomate filter order,
#     = length of filter - 1
function kaiserord( δ::Real, Δω )    
    A = -20*log10( δ )
    if A > 50
        α = 0.1102*( A - 8.7 )
    elseif A >= 21
        α = 0.5842*( A - 21 )^( 0.4 ) + 0.07886*( A - 21 )
    else
        α = 0
    end
    M = ( A - 7.95 )/( 2.285 * Δω )
    M = int(ceil(M))
    
    return ( M, α )
end
