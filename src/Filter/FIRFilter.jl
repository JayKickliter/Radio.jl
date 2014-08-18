module Multirate

export  FIRFilter,      polyize,
        filt!,          filt,
        decimate!,      decimate,
        interpolate!,   interpolate,
        resample!,      resample

#==============================================================================#
#                                    Types                                     #
#==============================================================================#

abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter h
type FIRStandard <: FIRKernel
    h::Vector
end

# Interpolator FIR kernel
type FIRInterpolator <: FIRKernel
    PFB::Matrix
    interpolation::Int
end

# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    h::Vector
    decimation::Int
end

# Rational resampler FIR kernel
type FIRRational  <: FIRKernel
    PFB::Matrix
    interploation::Int
    decimation::Int
end

type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    dlyLine::Vector
    reqDlyLineLen::Int
end

function FIRFilter( h::Vector, resampleRatio::Rational = 1//1 )
    interploation = num( resampleRatio )
    decimation    = den( resampleRatio )
    reqDlyLineLen = 0

    if resampleRatio == 1                                     # single-rate
        reqDlyLineLen = length( h ) - 1
        kernel        = FIRStandard( h )
    elseif interploation == 1                                 # decimate
        reqDlyLineLen = max( decimation-1, length(h)-1 )
        kernel        = FIRDecimator( h, decimation )
    elseif decimation == 1                                    # interpolate
        PFB           = polyize( h, interploation )
        reqDlyLineLen = size( PFB )[1] - 1
        kernel        = FIRInterpolator( PFB, interploation ) # TODO: this is a placeholder, needs to be figured out before stateful version can work
    else                                                      # rational
        PFB           = polyize( h, interploation )
        reqDlyLineLen = size(PFB)[1] - 1                      # TODO: this is a placeholder, needs to be figured out before stateful version can work
        kernel        = FIRRational( PFB, interploation, decimation )
    end

    dlyLine = zeros( eltype( h ), reqDlyLineLen )

    FIRFilter( kernel, dlyLine, reqDlyLineLen )
end




#==============================================================================#
#                           ____ _    _  _ ____ _  _                           #
#                           |___ |    |  | [__  |__|                           #
#                           |    |___ |__| ___] |  |                           #
#==============================================================================#

function flush( self::FIRFilter )
    dlyLine = self.dlyLine
    for i in 1:length( dlyLine )
        dlyLine[i] = 0
    end
end




#==============================================================================#
#                      _  _    ___ ____    ___  ____ ___                       #
#                      |__|     |  |  |    |__] |___ |__]                      #
#                      |  |     |  |__|    |    |    |__]                      #
#==============================================================================#

function polyize{T}( h::Vector{T}, interpolation )
    hLen      = length( h )
    hLenPerφ  = int( ceil( hLen/interpolation ))
    pfbSize   = hLenPerφ * interpolation
    # check that the vector is an integer multiple of interpolation
    if hLen != pfbSize
        hExtended             = similar( h, pfbSize )
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end
    nFilters  = interpolation
    hLen      = length( h )
    hLenPerφ  = int( hLen/nFilters )
    pfb       = reshape( h, nFilters, hLenPerφ )'
end




#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, dlyLine::Vector{T} = T[] )

    bufLen        = length( buffer )
    xLen          = length( x )
    hLen          = length( h )
    dlyLineLen    = length( dlyLine )
    reqDlyLineLen = hLen - 1

    bufLen >= xLen || error( "buffer length must be >= x length" )

    if dlyLineLen != reqDlyLineLen
        if dlyLineLen == 0
            dlyLine = zeros( T, reqDlyLineLen )
        elseif dlyLineLen < reqDlyLineLen
            dlyLine = [ zeros( T, reqDlyLineLen ), dlyLine ] # TODO: write the filtering logic to not depends on dlyLine being a certain length, as the current implementation allocates useless zeros
        else
            dlyLine = dlyLine[ end+1-reqDlyLineLen:end ]
        end
    end

    h = flipud( h )                                          # flip the h to make the multiplication more SIMD friendly

    for bufIdx in 1:hLen-1                                   # this first loop takes care of filter ramp up and previous dlyLine

        accumulator = zero(T)
        hIdx        = 1

        for dlyLineIdx in bufIdx:dlyLineLen                  # this loop takes care of previous dlyLine
            @inbounds accumulator += h[hIdx] * dlyLine[dlyLineIdx]
            hIdx += 1
        end

        hIdx = hLen-bufIdx+1
        for xIdx in 1:bufIdx                                 # this loop takes care of the first hlen-1 samples in x
            @inbounds accumulator += h[hIdx] * x[xIdx]
            hIdx += 1
        end

        @inbounds buffer[bufIdx] = accumulator
    end

    for bufIdx in hLen:xLen                                  # filter ramp is complete, normal filtering from this point on

        accumulator = zero(T)
        xIdx        = bufIdx-hLen+1

        @simd for hIdx in 1:hLen
            @inbounds accumulator += h[hIdx] * x[xIdx]
            xIdx += 1
        end

        @inbounds buffer[bufIdx] = accumulator
    end

    return buffer
end


filt{T}( h::Vector{T}, x::Vector{T}, dlyLine::Vector{T} = T[] ) = filt!( similar(x), h, x, dlyLine )

function filt( self::FIRFilter{FIRStandard}, x )
    h            = self.kernel.h
    hLen         = length( h )
    nextDlyLine  = x[end-hLen+2:end]
    y            = filt( self.kernel.h, x, self.dlyLine )
    self.dlyLine = nextDlyLine

    return y
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function interpolate!{T}( buffer::AbstractVector{T}, PFB::AbstractArray{T, 2}, x::AbstractVector{T}, dlyLine::AbstractVector{T} = T[] )

    (φLen, Nφ)    = size( PFB )                                                    # each column is a phase of the PFB, the rows hold the individual h
    xLen          = length( x )                                                    # number of input items
    bufLen        = length( buffer )
    reqDlyLineLen = φLen - 1
    dlyLineLen    = length( dlyLine )
    outLen        = xLen * Nφ
    criticalYidx  = min( reqDlyLineLen*Nφ, outLen )

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    if dlyLineLen != reqDlyLineLen
        if dlyLineLen == 0
            dlyLine = zeros( T, reqDlyLineLen )
        elseif dlyLineLen < reqDlyLineLen
            dlyLine = prepend!( dlyLine, zeros(  T, reqDlyLineLen - dlyLineLen ) ) # TODO: write the filtering logic to not depends on dlyLine being a certain length, as the current implementation allocates useless zeros
        else
            dlyLine = dlyLine[ end+1-reqDlyLineLen:end ]
        end
        dlyLineLen = length( dlyLine )
    end

    PFB      = flipud(PFB)
    inputIdx = 1
    φ        = 1

    for yIdx in 1:criticalYidx

        accumulator = zero(T)

        for k in 1:φLen-inputIdx
            @inbounds accumulator += PFB[k, φ] * dlyLine[k+inputIdx-1]
        end

        for k in 1:inputIdx
            @inbounds accumulator += PFB[φLen-inputIdx+k, φ] * x[k]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
    end

    for yIdx in criticalYidx+1:outLen

        accumulator = zero(T)

        for k in 1:φLen
            @inbounds accumulator += PFB[ k, φ ] * x[ inputIdx - φLen + k ]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
    end

    return buffer
end

interpolate{T}( PFB::Array{T, 2}, x::Vector{T}, dlyLine::AbstractVector{T} = T[] ) = interpolate!( similar( x, length(x)*size(PFB)[2] ), PFB, x, dlyLine )
interpolate( h, x, interpolation, dlyLine = eltype(x)[] )                          = interpolate( polyize( h, interpolation ), x, dlyLine )

function filt( self::FIRFilter{FIRInterpolator}, x )
   interpolate( self.kernel.PFB, x )
end

function filt( self::FIRFilter{FIRInterpolator}, x::AbstractVector )
    xLen          = length( x )
    interpolation = self.kernel.interpolation
    reqDlyLineLen = self.reqDlyLineLen
    
    y = interpolate( self.kernel.PFB, x, self.dlyLine )
    
    if xLen >= reqDlyLineLen
        self.dlyLine = x[ end-reqDlyLineLen+1 : end ]
    else
        self.dlyLine = [ self.dlyLine, x ][end-reqDlyLineLen+1:end]
    end
    return y
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function resample!{T}( buffer::Vector{T}, PFB::Array{T, 2}, x::Vector{T}, ratio::Rational; dlyLine::Vector{T} = T[], xStartIdx = 1 )

    (tapsPerφ, Nφ) = size( PFB )                                           
    interpolation  = num( ratio )
    decimation     = den( ratio )
    xLen           = length( x ) - xStartIdx + 1
    bufLen         = length( buffer )
    outLen         = int(ceil( xLen * interpolation / decimation ))
    criticalYidx   = int(floor( tapsPerφ * interpolation / decimation ))
    criticalYidx   = min( criticalYidx, outLen )

    1 <= xStartIdx <= length( x ) || error( "xStartIdx must be >= 1 and =< length( x ) " )
    bufLen >= outLen              || error( "length( buffer ) must be >= int(ceil( xLen * interpolation / decimation ))")

    PFB = flipud( PFB )
    
    for yIdx in 1:criticalYidx
        φIdx               = mod( (yIdx-1)*decimation, interpolation ) + 1
        inputIdx           = int( floor( (yIdx-1)*decimation / interpolation )) + xStartIdx
        accumulator        = zero(T)
                                        println(); println( "yIdx = $yIdx; φIdx = $φIdx; inputIdx = $inputIdx")
        for k in 1:inputIdx             
            @inbounds accumulator += PFB[ end-inputIdx+k, φIdx ] * x[ k ]
        end
        buffer[ yIdx ] = accumulator
    end
    
    for yIdx in criticalYidx+1:outLen
        φIdx        = mod( (yIdx-1)*decimation, interpolation ) + 1
        inputIdx    = int( floor( (yIdx-1)*decimation / interpolation )) + xStartIdx
        accumulator = zero(T)
        xFirstIdx   = inputIdx-tapsPerφ # this is actually ones less than the input of our first, with k added below it is the actuall index
                                        println(); println( "< yIdx = $yIdx; φIdx = $φIdx; inputIdx = $inputIdx >")
        for k in 1:tapsPerφ
            @inbounds accumulator += PFB[k, φIdx] * x[ xFirstIdx + k ]
        end
        
        buffer[ yIdx ] = accumulator
    end
    
    # TODO: figure out what the next inputIdx would be based off of th next yIdx.
    #       then figure out how many input samples are missing and return a new delay line 
    #       that is missing that number of samples

    return buffer
end

function resample{T}( PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )
    xLen   = length( x )
    bufLen = int( floor( xLen * ratio ))
    buffer = Array( T, bufLen )
    resample!( buffer, PFB, x, ratio )
end

resample{T}( h::Vector{T}, x::Vector{T}, ratio::Rational ) = resample( polyize(h, num(ratio)), x, ratio )




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function decimate!{T}( buffer::AbstractVector{T}, h::AbstractVector{T}, x::AbstractVector{T}, decimation::Integer, dlyLine::AbstractVector{T} = T[]; xStartIdx = 1 )

    xOffset       = xStartIdx - 1
    xLen          = length( x ) - xOffset
    hLen          = length( h )
    reqDlyLineLen = hLen-1
    dlyLineLen    = length( dlyLine )
    outLen        = int( ceil( xLen / decimation ))
    criticalYidx  = min( int(ceil(hLen / decimation)), outLen )

    1 <= xStartIdx <= length( x ) || error( "xStartIdx must be >= 1 and =< length( x ) " )
    length( buffer ) >= outLen    || error( "length(buffer) must be >= floor( ( length(x) ) / decimation)" )

    if dlyLineLen != reqDlyLineLen
        if dlyLineLen == 0
            dlyLine = zeros( T, reqDlyLineLen )
        elseif dlyLineLen < reqDlyLineLen
            dlyLine = prepend!( dlyLine, zeros(  T, reqDlyLineLen - dlyLineLen ) ) # TODO: write the filtering logic to not depends on dlyLine being a certain length, as the current implementation allocates useless zeros
        else
            dlyLine = dlyLine[ end+1-reqDlyLineLen:end ]
        end
        dlyLineLen = length( dlyLine )
    end

    h = flipud(h)                                                         # TODO: figure out a way to not always have to flip taps each time
    inputIdx = 1                                                          # inputIdx is is the current input, not the actual index of of the current x in the delay line
    for yIdx in 1:criticalYidx                                            # Filtering is broken up into two outer loops
                                                                          # This first loop takes care of filter ramp up.
        hIdx        = 1
        accumulator = zero(T)

        for dlyLineIdx in inputIdx:dlyLineLen                             # Inner Loop 1: Handles convolution of taps and delay line
            accumulator += h[hIdx] * dlyLine[dlyLineIdx]
            hIdx += 1
        end

        for k in 1:inputIdx                                               # Inner Loop 2: handles convolution of taps and x
            accumulator += h[ hIdx ] * x[ k + xOffset ]
            hIdx += 1
        end

        buffer[yIdx] = accumulator
        inputIdx += decimation
    end

    inputIdx -= hLen

    for yIdx in criticalYidx+1:outLen                                     # second outer loop, we are now in the clear to to convolve without x[inputIdx-]
        accumulator  = zero(T)

        for k in 1:hLen
            accumulator += h[ k ] * x[ inputIdx + k + xOffset ]
        end

        buffer[yIdx] = accumulator
        inputIdx += decimation
    end

    inputIdx += hLen - decimation
    xLeftoverLen = max( xLen - inputIdx, 0 )

    return buffer, xLeftoverLen                                           # TODO: return number of elements written to buffer
end


function decimate{T}( h::Vector{T}, x::AbstractVector{T}, decimation::Integer, dlyLine::Vector{T} = T[]; xStartIdx = 1 )
    xLen   = length( x ) - xStartIdx + 1
    buffer = similar( x, int(ceil( xLen / decimation )) )
    decimate!( buffer, h, x, decimation, dlyLine, xStartIdx = xStartIdx )
end

function filt{T}( self::FIRFilter{FIRDecimator}, x::Vector{T} )
    xLen          = length( x )
    h             = self.kernel.h
    decimation    = self.kernel.decimation
    dlyLine       = self.dlyLine
    dlyLineLen    = length( dlyLine )
    reqDlyLineLen = self.reqDlyLineLen
    combinedLen   = dlyLineLen + xLen
    y             = T[]

    if combinedLen > reqDlyLineLen
        xStartIdx         = reqDlyLineLen - dlyLineLen + 1
        self.dlyLine      = [ dlyLine, x[1:xStartIdx-1] ]
        (y, xLeftoverLen) = decimate( h, x, decimation, self.dlyLine; xStartIdx = xStartIdx )
        nextDlyLineLen    = reqDlyLineLen - decimation + 1 + xLeftoverLen
        self.dlyLine      = [ self.dlyLine, x[xStartIdx:end] ][end-nextDlyLineLen+1:end]
    else
        append!( dlyLine, x )
    end

    return y
end

end # module
