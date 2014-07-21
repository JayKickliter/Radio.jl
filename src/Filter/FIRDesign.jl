import DSP: firfilt

#==============================================================================#
#                                    Types                                     #
#==============================================================================#
abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter taps
type FIRKernelSingleRate <: FIRKernel
    taps::Vector
end

# Multirate FIR kernels, holds a polyphase filter bank
type FIRKernelInterpolator <: FIRKernel
    PFB::Matrix
    interploation::Int
end

type FIRFilter <: Filter
    kernel::FIRKernel
    state::Vector
end


function FIRFilter{Tx}( ::Type{Tx}, taps::Vector )
    Nt     = length( taps )
    state  = zeros( Tx, Nt-1 )
    kernel = FIRKernelSingleRate( taps )
    FIRFilter( kernel, state )
end

FIRFilter{Tt}( taps::Vector{Tt} ) = FIRFilter( Tt, taps )

function FIRFilter{Tx}( ::Type{Tx}, taps::Vector, interploation::Integer )
    PFB    = polyize( taps, interploation )
    Ns     = size( PFB )[1]
    state  = zeros( Tx, Ns )
    kernel = FIRKernelInterpolator( PFB, interploation )
    FIRFilter( kernel, state )
end

FIRFilter{Tt}( taps::Vector{Tt}, interploation::Integer ) = FIRFilter( Tt, taps, interploation )

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
# N       = Filter order; Returns M (N+1) samples
#         = If FIRType is HIGH_PASS, may return M+1 (N+2) samples to make it
#           a type 1 filter.
# F       = Cutoff frequency. Real for high-pass & lowpass, Vector{Real} for
#           band-pass & band-reject
# FIRType = Type of filter: FIR_TYPE.LOW_PASS, FIR_TYPE.BAND_PASS,
#           FIR_TYPE.BAND_REJECT
# 
# Reference(s): [3]
function firprototype( N::Integer, F::Union(Real, Vector), FIRType::Integer )
    if     FIRType == FIR_TYPE.LOW_PASS
        return [ 2*F*sinc(2*F*(n-N/2)) for n = 0:N ]
    elseif FIRType == FIR_TYPE.BAND_PASS
        return [ 2*(F[1]*sinc(2*F[1]*(n-N/2)) - F[2]*sinc(2*F[2]*(n-N/2))) for n = 0:N ]        
    elseif FIRType == FIR_TYPE.HIGH_PASS
        N = isodd( N ) ? N+1 : N 
        return [ sinc(n-N/2) - 2*F*sinc(2*F*(n-N/2)) for n = 0:N ]        
    elseif FIRType == FIR_TYPE.BAND_STOP
        return [ 2*(F[2]*sinc(2*F[2]*(n-N/2)) - F[1]*sinc(2*F[1]*(n-N/2))) for n = 0:N ]
    else
        error("Not a valid FIR_TYPE")
    end
end



#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# N      = Filter order; Returns N+1 samples
# F_c    = Transition frequency.
# window = Window function to apply to ideal truncated sinc response
# F_s    = Sample rate. Defualts to 2.
function firdes( N::Integer, F_c::Real, windowFunction::Function, F_s::Real=2.0 )
    # TODO: add argument valdiation
    F_t = F_c/F_s
    if F_t < 0 || F_t > 0.5
        error("F_c/F_s must be > 0.0 and < 0.5")
    end
    
    [ 2*F_t*sinc(2*F_t*(n-N/2)) * windowFunction(n, N) for n = 0:N ]
end

# F_c    = Transition frequency.
# window = A vector containging pre-calculated window coefficients
function firdes( F_c::Real, window::Vector, F_s=2.0 )
    # TODO: add argument valdiation
    F_t = F_c/F_s
    if F_t < 0 || F_t > 0.5
        error("F_c/F_s must be > 0.0 and < 0.5")
    end
    
    N = length( window ) - 1
    
    [ 2*F_t*sinc(2*F_t*(n-N/2)) * window[n+1] for n = 0:N ]
end


#==============================================================================#
#                                Kaiser Order                                  #
#==============================================================================#
# δ  = ripple in passband and stopband
# Δω = ω_p - ω_s
#    = transition width in radians 
# 
# returns:
#   N = approxomate filter order,
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
    N = int(ceil( ( A - 7.95 )/( 2.285 * Δω )))

    return ( N, β )
end