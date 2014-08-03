#==============================================================================#
#                                    Types                                     #
#==============================================================================#
abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter taps
type FIRStandard <: FIRKernel
    taps::Vector
end

# Interpolator FIR kernel
type FIRInterpolator <: FIRKernel
    PFB::Matrix
    interploation::Int
    leftovers::Vector
end

# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    taps::Vector
    decimation::Int
    xLeftover::Vector
end

# Rational resampler FIR kernel
type FIRRational  <: FIRKernel
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
    kernel = FIRStandard( taps )
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
        kernel = FIRDecimator( PFB, decimation, Tx[] )
    elseif decimation == 1
        PFB    = polyize( taps, interploation )
        kernel = FIRInterpolator( PFB, interploation, Tx[] )
    else
        PFB    = polyize( taps, interploation )
        kernel = FIRRational( PFB, interploation, decimation, Tx[] )
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

function filt{T}( h::Vector{T}, x::Vector{T} )
    xLen = length( x )
    hLen = length( h )
    @assert length(x) > length(h)

    buffer  = zeros( T, xLen )

    for n = 1:hLen-1
        for m = 1:n
            @inbounds buffer[n] += h[m] * x[n-m+1]
        end
    end
    for n = hLen:xLen
        @simd for m = 1:hLen
            @inbounds buffer[n] += h[m] * x[n-m+1]
        end
    end
    buffer
end

#=
# Short single-rate test
h          = rand( 56 );
x          = rand( 1000 );

nativeResult = filt( h, x );
baseResult   = Base.filt( h, 1.0, x );

[ baseResult nativeResult ]
areApprox( nativeResult, baseResult )
=#


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

function interpolate!{T}( buffer::Vector{T}, PFB::Array{T, 2}, x::Vector{T} )
    (hLen, Nφ)  = size( PFB )      # each column is a phase of the PFB, the rows hold the individual taps
    xLen        = length( x )      # number of input items
    bufLen      = length( buffer )
    
    bufLen >= xLen*Nφ || error( "buffer length must be >= signal length")
    
    yIdx = 1
    
    for xIdx = 1:hLen-1, φ = 1:Nφ  # until xIdx == hLen, x[xIdx-m+1] would reach out of bounds 
                                   # this first loop limits xIdx-m+1 to a minimum of 1
        accumulator = zero(T)
        for Tn = 1:xIdx            # for each tap in phase[n]
            @inbounds accumulator += PFB[Tn, φ] * x[xIdx-Tn+1]
        end
        buffer[yIdx] = accumulator
        yIdx += 1
    end
    
    PFB = flipud(PFB)
    
    for xIdx = hLen:xLen, φ = 1:Nφ # no longer in danger of stepping out of bounds
        
        xStartIdx      = xIdx-hLen               
        yIdx          = Nφ*(xIdx-1)+φ
        accumulator = zero(T)        
                
        @simd for Tn = 1:hLen
            @inbounds accumulator += PFB[Tn, φ] * x[xStartIdx + Tn]
        end
        
        buffer[yIdx] = accumulator
        yIdx += 1
    end

    return buffer
end

interpolate{T}( PFB::Array{T, 2}, x::Vector{T} ) = interpolate!( similar( x, length(x)*size(PFB)[2] ), PFB, x )
interpolate( h, x, interpolation )               = interpolate( polyize( h, interpolation ), x )

function filt( self::FIRFilter{FIRInterpolator}, x )
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

function resample!{T}( buffer::Vector{T}, PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )   
    (hLen, Nφ)    = size( PFB ) # each column is a phase of the PFB, the rows hold the individual taps
    interpolation = num(ratio)
    decimation    = den( ratio )
    xLen          = length( x ) # number of input items    
    bufLen        = length( buffer )

    
    xLen * interpolation % decimation == 0 || error("signal length * interpolation mod decimation must be 0")    
        
    for m = 0:bufLen-1
        
        φ    = mod( m*decimation, interpolation)
        nm   = int( floor( m*decimation / interpolation ))
        kMax = nm < hLen ? nm+1 : hLen
        accumulator  = zero(T)        
        
        for k = 0:kMax-1
            accumulator += PFB[ k+1, φ+1 ] * x[ nm+1-k ]
        end
                
        buffer[m+1] = accumulator
    
    end

    return buffer
end

function resample{T}( PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )
    xLen   = length( x )    
    bufLen = int( floor( xLen * ratio ))
    buffer = Array( T, bufLen )
    resample!( buffer, PFB, x, ratio )
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

function decimate!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, decimation::Integer, state::Vector{T} = T[] )
    xLen   = length( x )
    hLen   = length( h )
    outLen = floor(int(xLen / decimation))

    length( buffer ) * decimation >= xLen || error( "buffer lenght must be >= signal length * decimation" )        
    
    criticalYidx = int(ceil(hLen / decimation)) # The index of y where our taps would overlap
    xIdx         = 1    

    for yIdx = 1:criticalYidx
            
        accumulator  = zero(T)
        kMax = xIdx < hLen ? xIdx : hLen
                
        for k = 1:kMax
            @inbounds accumulator += h[ k ] * x[ xIdx+1-k ]
        end
        
        @inbounds buffer[yIdx] = accumulator
        xIdx += decimation
    end
    
    h = flipud(h)
    
    xIdx -= hLen
    
    for yIdx = criticalYidx+1:outLen
            
        accumulator  = zero(T)
                
        for k = 1:hLen
            @inbounds accumulator += h[ k ] * x[ xIdx+k ]
        end
        
        @inbounds buffer[yIdx] = accumulator
        xIdx += decimation
    end
        
    return buffer
end

decimate{T}( h::Vector{T}, x::Vector{T}, decimation::Integer ) = decimate!( similar(x, int(floor(length(x)/decimation))), h, x, decimation )

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
