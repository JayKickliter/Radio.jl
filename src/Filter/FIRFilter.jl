#==============================================================================#
#                                    Types                                     #
#==============================================================================#
abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter taps
type FIRKernelSingleRate <: FIRKernel
    taps::Vector
end

# Interpolator FIR kernel
type FIRKernel⬆︎ <: FIRKernel
    PFB::Matrix
    interploation::Int
    leftovers::Vector
end

# Decimator FIR kernel
type FIRKernel⬇︎ <: FIRKernel
    taps::Vector
    decimation::Int
    xLeftover::Vector
end

# Rational resampler FIR kernel
type FIRKernel⬆︎⬇︎  <: FIRKernel
    PFB::Matrix
    interploation::Int
    decimation::Int
    xLeftover::Vector
end

type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    state::Vector
end


function FIRFilter{Tx}( ::Type{Tx}, taps::Vector )
    Nt     = length( taps )
    state  = zeros( Tx, Nt-1 )
    kernel = FIRKernelSingleRate( taps )
    FIRFilter( kernel, state )
end

FIRFilter{Tt}( taps::Vector{Tt} ) = FIRFilter( Tt, taps )

function FIRFilter{Tx}( ::Type{Tx}, taps::Vector, resampleRatio::Rational )
    
    interploation = num( resampleRatio )
    decimation    = den( resampleRatio )
    
    if resampleRatio == 1
        return FIRFilter( Tx, taps )        
    elseif interploation == 1
        PFB    = taps
        kernel = FIRKernel⬇︎( PFB, decimation, Tx[] )
    elseif decimation == 1
        PFB    = polyize( taps, interploation )
        kernel = FIRKernel⬆︎( PFB, interploation, Tx[] )
    else
        PFB    = polyize( taps, interploation )
        kernel = FIRKernel⬆︎⬇︎( PFB, interploation, decimation, Tx[] )
    end 
    
    hLen  = size(PFB)[1]
    state = zeros( Tx, hLen - 1 )         
        
    FIRFilter( kernel, state )
end

FIRFilter{Tt}( taps::Vector{Tt}, resampleRatio::Rational ) = FIRFilter( Tt, taps, resampleRatio )

#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function filt{T}( taps::Vector{T}, signal::Vector{T} )
    @assert length(signal) > length(taps)
    signal_len     = length( signal )
    taps_len       = length( taps )
    output_buffer  = zeros( T, signal_len )
    for n = 1:taps_len-1
        for m = 1:n
            @inbounds output_buffer[n] += taps[m] * signal[n-m+1]
        end
    end
    for n = taps_len:signal_len
        @simd for m = 1:taps_len
            @inbounds output_buffer[n] += taps[m] * signal[n-m+1]
        end
    end
    output_buffer
end




#==============================================================================#
#                      _  _    ___ ____    ___  ____ ___                       #
#                      |__|     |  |  |    |__] |___ |__]                      #
#                      |  |     |  |__|    |    |    |__]                      #
#==============================================================================#

function polyize{T}( h::Vector{T}, interpolation )
    hLen         = length( h )
    tapsPerPhase = int( ceil( hLen/interpolation ))
    pfbSize = tapsPerPhase * interpolation
    # check that the vector is an integer multiple of interpolation
    if hLen != pfbSize 
        hExtended             = similar( h, pfbSize )
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end
    nFilters     = interpolation
    hLen         = length( h )
    tapsPerPhase = int( hLen/nFilters )
    pfb          = reshape( h, nFilters, tapsPerPhase )'
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function interpolate{T}( PFB::Array{T, 2}, x::Vector{T} )
    (hLen, Nφ) = size( PFB )           # each column is a phase of the PFB, the rows hold the individual taps
    xLen        = length( x )           # number of input items
    yLen        = xLen * Nφ
    y           = similar( x, xLen * Nφ ) # yLen = xLen * Nφ
    
    for Xn = 1:hLen-1, φ = 1:Nφ        # until Xn == hLen, x[Xn-m+1] would reach out of bounds 
                                        # this first loop limits Xn-m+1 to a minimum of 1
        Yn          = Nφ*(Xn-1)+φ
        accumulator = zero(T)
        for Tn = 1:Xn                   # for each tap in phase[n]
            @inbounds accumulator += PFB[Tn, φ] * x[Xn-Tn+1]
        end
        y[Yn] = accumulator
    end
    
    PFB = flipud(PFB)
    
    for Xn = hLen:xLen, φ = 1:Nφ         # no longer in danger of stepping out of bounds
        
        XnBase      = Xn-hLen               
        Yn          = Nφ*(Xn-1)+φ
        accumulator = zero(T)        
                
        @simd for Tn = 1:hLen
            @inbounds accumulator += PFB[Tn, φ] * x[XnBase + Tn]
        end
        
        y[Yn] = accumulator
    end

    return y
end

interpolate( h, x, interpolation ) = interpolate( polyize( h, interpolation ), x )

function filt( self::FIRFilter{FIRKernel⬆︎}, x )
   interpolate( self.kernel.PFB, x ) 
end

#=
# Short interpolate test
factor = 4;
h      = rand( 56 );
x      = rand( 1000 );
xx     = zeros( length(x) * factor );
for n = 0:length(x)-1;
    xx[ n*factor+1 ] = x[ n+1 ];
end

nativeResult = interpolate( h, x, factor );
baseResult   = Base.filt( h, 1.0, xx );
[ baseResult nativeResult ]
areApprox( nativeResult, baseResult )
=#



#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function resample{T}( PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )
    
    (hLen, Nφ)    = size( PFB ) # each column is a phase of the PFB, the rows hold the individual taps
    interpolation = num(ratio)
    decimation    = den( ratio )
    xLen          = length( x )        # number of input items    
    
    xLen * interpolation % decimation == 0 || error("signal length * interpolation mod decimation must be 0")    
    
    yLen = int(xLen*interpolation/decimation)
    y    = zeros( T, yLen )
        
    for m = 0:yLen-1
        
        φ    = mod( m*decimation, interpolation)
        nm   = int( floor( m*decimation / interpolation ))
        kMax = nm < hLen ? nm+1 : hLen
        acc  = zero(T)        
        
        for k = 0:kMax-1
            acc += PFB[ k+1, φ+1 ] * x[ nm+1-k ]
        end
                
        y[m+1] = acc
    
    end

    return y
end

resample{T}( h::Vector{T}, x::Vector{T}, ratio::Rational ) = resample( polyize(h, num(ratio)), x, ratio )

#=
# Short rational resample test
upfactor   = 3;
downfactor = 4;
h          = rand( 56 );
x          = rand( 1000 );
xx         = zeros( length(x) * upfactor );
baseResult = similar( x, int( length(x) * upfactor / downfactor ))

for n = 0:length(x)-1;
    xx[ n*upfactor+1 ] = x[ n+1 ];
end

nativeResult           = resample( h, x, upfactor//downfactor );
baseResultInterpolated = Base.filt( h, 1.0, xx );
baseResult = [ baseResultInterpolated[n] for n = 1:downfactor:length( baseResultInterpolated ) ]

[ baseResult nativeResult ]
areApprox( nativeResult, baseResult )
=#




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function decimate!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, decimation::Integer )
    (xLen   = length( x )) % decimation      == 0    || error( "signal length % decimation must be 0" )
    (bufLen = length( buffer )) * decimation >= xLen || error( "buffer lenght must be >= signal length * decimation" )        
    hLen    = length( h )
    yLen    = int(xLen / decimation)    
    
    
    criticalYidx = int(ceil(hLen / decimation)) # The index of y where our taps would overlap
    xIdx = 1
    
    for yIdx = 1:criticalYidx
            
        acc  = zero(T)
        kMax = xIdx < hLen ? xIdx : hLen
                
        for k = 1:kMax
            @inbounds acc += h[ k ] * x[ xIdx+1-k ]
        end
        
        @inbounds buffer[yIdx] = acc
        xIdx += decimation
    end
    
    h = flipud(h)
    
    xIdx -= hLen
    
    for yIdx = criticalYidx+1:yLen
            
        acc  = zero(T)
                
        for k = 1:hLen
            @inbounds acc += h[ k ] * x[ xIdx+k ]
        end
        
        @inbounds buffer[yIdx] = acc
        xIdx += decimation
    end
        
    return buffer
end

decimate{T}( h::Vector{T}, x::Vector{T}, decimation::Integer ) = decimate!( similar(x, int(length(x)/decimation)), h::Vector{T}, x::Vector{T}, decimation::Integer )

#=
# Short decimation test
h            = rand( 56 )
x            = rand( 1000 )
factor       = 10
nativeResult = decimate( h, x, factor )
baseResult   = Base.filt( h, 1.0, x )
baseResult   = baseResult[1:factor:end]
areApprox( nativeResult, baseResult )
=#
