#==============================================================================#
#                               Constants/Enums                                #
#==============================================================================#
module FIR_TYPE
const LOW_PASS  = 0
const BAND_PASS = 1
const HIGH_PASS = 2
const BAND_STOP = 3
end


#==============================================================================#
#                             Prototype FIR Filter                             #
#==============================================================================#
# M       = Filter order; Returns M+1 samples
# F       = Cutoff frequency. Real for high-pass & lowpass, Vector{Real} for
#           band-pass & band-reject
# FIRType = Type of filter: FIR_TYPE.LOW_PASS, FIR_TYPE.BAND_PASS,
#           FIR_TYPE.BAND_REJECT
function firprototype( M::Integer, F::Union(Real, Vector), FIRType::Integer )
    if     FIRType == FIR_TYPE.LOW_PASS
        return [ 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]
    elseif FIRType == FIR_TYPE.BAND_PASS
        return [ 2*(F[1]*sinc(2*F[1]*(n-M/2)) - F[2]*sinc(2*F[2]*(n-M/2))) for n = 0:M ]        
    elseif FIRType == FIR_TYPE.HIGH_PASS
        return [ sinc(n-M/2) - 2*F*sinc(2*F*(n-M/2)) for n = 0:M ]        
    elseif FIRType == FIR_TYPE.BAND_STOP
        return [ 2*(F[2]*sinc(2*F[2]*(n-M/2)) - F[1]*sinc(2*F[1]*(n-M/2))) for n = 0:M ]
        
    else
        error("Not a valid FIR_TYPE")
    end
end



#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# M      = Filter order; Returns M+1 samples
# F_c    = Transition frequency.
# window = Window function to apply to ideal truncated sinc response
# F_s    = Sample rate. Defualts to 2.
function firdes( M::Integer, F_c::Real, windowFunction::Function, F_s::Real=2.0 )
    # TODO: add argument valdiation
    F_t = F_c/F_s
    if F_t < 0 || F_t > 0.5
        error("F_c/F_s must be > 0.0 and < 0.5")
    end
    
    [ 2*F_t*sinc(2*F_t*(n-M/2)) * windowFunction(n, M) for n = 0:M ]
end

# F_c    = Transition frequency.
# window = A vector containging pre-calculated window coefficients
function firdes( F_c::Real, window::Vector, F_s=2.0 )
    # TODO: add argument valdiation
    F_t = F_c/F_s
    if F_t < 0 || F_t > 0.5
        error("F_c/F_s must be > 0.0 and < 0.5")
    end
    
    M = length( window ) - 1
    
    [ 2*F_t*sinc(2*F_t*(n-M/2)) * window[n+1] for n = 0:M ]
end


#==============================================================================#
#                                Kaiser Order                                  #
#==============================================================================#
# δ  = ripple in passband and stopband
# Δω = ω_p - ω_s
#    = transition width in radians 
# 
# returns:
#   M = approxomate filter order,
#     = length of filter - 1
function kaiserord( δ::Real, Δω::Real )
    # TODO: add argument valdiation
    A = -20*log10( δ )
    if A > 50
        β = 0.1102*( A - 8.7 )
    elseif A >= 21
        β = 0.5842*( A - 21 )^( 0.4 ) + 0.07886*( A - 21 )
    else
        β = 0.0
    end
    M = int(ceil( ( A - 7.95 )/( 2.285 * Δω )))

    return ( M, β )
end
