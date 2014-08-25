module Multirate

import Base: filt, filt!, reset

export  FIRFilter,      polyize,
        filt!,          filt,
        decimate!,      decimate,
        interpolate!,   interpolate,
        resample!,      resample,
        reset

#==============================================================================#
#                                    Types                                     #
#==============================================================================#

typealias PFB{T} Array{T,2}

abstract Filter
abstract FIRKernel

# Single rate FIR kernel, just hold filter h
type FIRStandard <: FIRKernel
    h::Vector
    hLen::Int
end

# Interpolator FIR kernel
type FIRInterpolator <: FIRKernel
    pfb::PFB
    interpolation::Int
    Nφ::Int
    tapsPerφ::Int
end

# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    h::Vector
    hLen::Int
    decimation::Int
end

# Rational resampler FIR kernel
type FIRRational  <: FIRKernel
    pfb::PFB
    ratio::Rational
    Nφ::Int
    tapsPerφ::Int
    φIdx::Int
end

type FIRFilter{Tk<:FIRKernel} <: Filter
    kernel::Tk
    dlyLine::Vector
    reqDlyLineLen::Int
end

function FIRFilter( h::Vector, resampleRatio::Rational = 1//1 )
    interpolation = num( resampleRatio )
    decimation    = den( resampleRatio )
    reqDlyLineLen = 0
    hLen          = length( h )

    if resampleRatio == 1                                     # single-rate
        reqDlyLineLen = length( h ) - 1
        h             = flipud( h )
        kernel        = FIRStandard( h, hLen )
    elseif interpolation == 1                                 # decimate
        reqDlyLineLen = max( decimation-1, length(h)-1 )
        h             = flipud( h )
        kernel        = FIRDecimator( h, hLen, decimation )
    elseif decimation == 1                                    # interpolate
        pfb              = flipud(polyize( h, interpolation ))
        ( tapsPerφ, Nφ ) = size( pfb )
        reqDlyLineLen    = tapsPerφ - 1
        kernel           = FIRInterpolator( pfb, interpolation, Nφ, tapsPerφ )
    else                                                      # rational
        pfb              = polyize( h, interpolation )
        ( tapsPerφ, Nφ ) = size( pfb )
        reqDlyLineLen    = tapsPerφ - 1
        kernel           = FIRRational( pfb, resampleRatio, Nφ, tapsPerφ, 1 )
    end

    dlyLine = zeros( eltype( h ), reqDlyLineLen )

    FIRFilter( kernel, dlyLine, reqDlyLineLen )
end




#==============================================================================#
#                            ____ ____ ____ ____ ___                           #
#                            |__/ |___ [__  |___  |                            #
#                            |  \ |___ ___] |___  |                            #
#==============================================================================#

# Resets filter and its kernel to an initial state

# Does nothing for non-rational kernels
reset( self::FIRKernel ) = self

# For rational kernel, set φIdx back to 1
reset( self::FIRRational ) = self.φIdx = 1

# For FIRFilter, set delay line to zeros of same tyoe and required length
function reset( self::FIRFilter )
    self.dlyLine = zeros( eltype( self.dlyLine ), self.reqDlyLineLen )
    reset( self.kernel )
    return self
end




#==============================================================================#
#                      _  _    ___ ____    ___  ____ ___                       #
#                      |__|     |  |  |    |__] |___ |__]                      #
#                      |  |     |  |__|    |    |    |__]                      #
#==============================================================================#

# Converts a vector of coefficients to a matrix. Each column is a filter.
# Appends zeros if necessary.
# Example:
#   julia> polyize( [1:9], 4 )
#   3x4 Array{Int64,2}:
#    1  2  3  4
#    5  6  7  8
#    9  0  0  0

function polyize{T}( h::Vector{T}, numFilters::Integer )
    hLen      = length( h )
    hLenPerφ  = iceil(  hLen/numFilters  )
    pfbSize   = hLenPerφ * numFilters

    if hLen != pfbSize                                # check that the vector is an integer multiple of numFilters
        hExtended             = similar( h, pfbSize ) # No? extend and zero pad
        hExtended[1:hLen]     = h
        hExtended[hLen+1:end] = 0
        h                     = hExtended
    end

    hLen      = length( h )
    hLenPerφ  = int( hLen/numFilters )
    pfb       = reshape( h, numFilters, hLenPerφ )'
end

#==============================================================================#
#               ____ _  _ ___ ___  _  _ ___    _    ____ _  _                  #
#               |  | |  |  |  |__] |  |  |     |    |___ |\ |                  #
#               |__| |__|  |  |    |__|  |     |___ |___ | \|                  #
#==============================================================================#

# Calculates the resulting length of a multirate filtering operation, given
#   resampling ratio, input length, and initial phase of the filter bank.
#
# It's hard to explain how this works without a diagram.

function outputlength( ratio::Rational, inputLen, φ = 1 )
    interpolation = num( ratio )
    decimation    = den( ratio )
    outLen        = (( inputLen * interpolation ) - φ + 1 ) / decimation
    iceil(  outLen  )
end


#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

# Stateful single-rate filtering with pre-allocated buffer
function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRStandard}, x::Vector{T} )
    dlyLine::Vector{T} = self.dlyLine
    h::Vector{T}       = self.kernel.h
    hLen               = self.kernel.hLen
    reqDlyLineLen      = self.reqDlyLineLen
    bufLen             = length( buffer )
    xLen               = length( x )
    outLen             = xLen
    criticalYidx       = min( hLen, outLen )

    bufLen >= xLen || error( "buffer length must be >= x length" )

    for yIdx in 1:criticalYidx                                   # this first loop takes care of filter ramp up and previous dlyLine
        accumulator = zero(T)

        for k in 1:hLen-yIdx
            @inbounds accumulator += h[k] * dlyLine[k+yIdx-1]
        end

        for k in 1:yIdx
            @inbounds accumulator += h[hLen-yIdx+k] * x[k]
        end

        @inbounds buffer[yIdx] = accumulator
    end

    for yIdx in criticalYidx+1:xLen
        accumulator = zero(T)

        for k in 1:hLen
            @inbounds accumulator += h[k] * x[yIdx-hLen+k]
        end

        @inbounds buffer[yIdx] = accumulator
    end

    if xLen >= reqDlyLineLen
        self.dlyLine = x[end-reqDlyLineLen+1:end]
    else
        self.dlyLine = [ dlyLine, x ][end-reqDlyLineLen+1:end]
    end

    return buffer
end

function filt{T}( self::FIRFilter{FIRStandard}, x::Vector{T} )
    buffer = zeros( eltype(x), length(x) )
    filt!( buffer, self, x )
end


#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRInterpolator}, x::Vector{T} )
    pfb::PFB{T}        = self.kernel.pfb
    dlyLine::Vector{T} = self.dlyLine
    interpolation      = self.kernel.interpolation
    Nφ                 = self.kernel.Nφ
    tapsPerφ           = self.kernel.tapsPerφ
    xLen               = length( x )
    bufLen             = length( buffer )
    reqDlyLineLen      = self.reqDlyLineLen
    outLen             = xLen * interpolation
    criticalYidx       = min( reqDlyLineLen*interpolation, outLen )

    bufLen >= outLen || error( "length( buffer ) must be >= interpolation * length(x)")

    inputIdx = 1
    φ        = 1

    for yIdx in 1:criticalYidx

        accumulator = zero(T)

        for k in 1:tapsPerφ-inputIdx
            @inbounds accumulator += pfb[k, φ] * dlyLine[k+inputIdx-1]
        end

        for k in 1:inputIdx
            @inbounds accumulator += pfb[tapsPerφ-inputIdx+k, φ] * x[k]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
    end

    for yIdx in criticalYidx+1:outLen

        accumulator = zero(T)

        for k in 1:tapsPerφ
            @inbounds accumulator += pfb[ k, φ ] * x[ inputIdx - tapsPerφ + k ]
        end

        @inbounds buffer[yIdx]  = accumulator
        (φ, inputIdx) = φ == Nφ ? ( 1, inputIdx+1 ) : ( φ+1, inputIdx )
    end

    if xLen >= reqDlyLineLen
        self.dlyLine = x[end-reqDlyLineLen+1:end]
    else
        self.dlyLine = [ dlyLine, x ][end-reqDlyLineLen+1:end]
    end

    return buffer
end

function filt( self::FIRFilter{FIRInterpolator}, x::Vector )
    xLen   = length( x )
    bufLen = xLen * self.kernel.interpolation
    buffer = similar( x, bufLen )
    filt!( buffer, self, x )
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRRational}, x::Vector{T} )
    xLen          = length( x )
    pfb           = self.kernel.pfb
    ratio         = self.kernel.ratio
    dlyLine       = self.dlyLine
    dlyLineLen    = length( dlyLine )
    φIdx          = self.kernel.φIdx
    reqDlyLineLen = self.reqDlyLineLen
    combinedLen   = dlyLineLen + xLen

    if combinedLen <= reqDlyLineLen
        append!( dlyLine, x )
        return T[]
    end

    xStartIdx     = reqDlyLineLen - dlyLineLen + 1
    xOffset       = xStartIdx - 1
    xLen          = length( x ) - xOffset
    outLen        = outputlength( ratio, xLen, φIdx )
    dlyLine       = [ dlyLine, x[1:xStartIdx-1] ]
    criticalYidx  = min( ifloor( self.kernel.tapsPerφ * ratio ), outLen )
    interpolation = num( ratio )
    decimation    = den( ratio )
    φIdxStepSize  = mod( decimation, interpolation )
    Nφ            = self.kernel.Nφ
    tapsPerφ      = self.kernel.tapsPerφ
    criticalφIdx  = Nφ - φIdxStepSize

    yIdx         = 0
    inputIdx     = 1
    φIdxLast     = 0
    inputIdxLast = 0

    for yIdx in 1:criticalYidx
        accumulator = zero(T)

        for k in 1:tapsPerφ-inputIdx
            accumulator += pfb[ k, φIdx ] * dlyLine[ k+inputIdx-1]
        end

        for k in 1:inputIdx
            accumulator += pfb[ end-inputIdx+k, φIdx ] * x[ k+xOffset ]
        end

        inputIdxLast = inputIdx                                                             # TODO: get rid of need to store last
        φIdxLast     = φIdx                                                                 # TODO: get rid of need to store last
        inputIdx    += int( floor( ( φIdx + decimation - 1 ) / interpolation ) )
        φIdx         = φIdx > criticalφIdx ? φIdx + φIdxStepSize - Nφ : φIdx + φIdxStepSize #

        buffer[ yIdx ] = accumulator
    end

    for yIdx in criticalYidx+1:outLen
        accumulator = zero(T)
        xFirstIdx   = inputIdx-tapsPerφ                                                     # this is actually ones less than the input of our first, with k added below it is the actuall index

        for k in 1:tapsPerφ
            accumulator += pfb[k, φIdx] * x[ k + xFirstIdx + xOffset ]
        end

        inputIdxLast   = inputIdx                                                           # TODO: get rid of need to store last
        φIdxLast       = φIdx                                                               # TODO: get rid of need to store last
        inputIdx      += int( floor( ( φIdx + decimation - 1 ) / interpolation ) )
        φIdx           = φIdx > criticalφIdx ? φIdx + φIdxStepSize - Nφ : φIdx + φIdxStepSize
        buffer[ yIdx ] = accumulator
    end

    xLeftoverLen     = length( x ) - xOffset - inputIdxLast
    sampleDeficit    = inputIdx - inputIdxLast - xLeftoverLen
    dlyLineLen       = reqDlyLineLen - sampleDeficit + 1    
    self.dlyLine     = [ dlyLine, x[xStartIdx:end] ][end-dlyLineLen+1:end]
    self.kernel.φIdx = φIdx
    
    return buffer
end

function filt{T}( self::FIRFilter{FIRRational}, x::Vector{T} )
    xLen   = max( length( x ) - self.reqDlyLineLen + length( self.dlyLine ), 0 )
    outLen = outputlength( self.kernel.ratio, xLen, self.kernel.φIdx )
    buffer = similar( x, outLen )
    filt!( buffer, self, x )        
end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

# Stateful decimation
function filt!{T}( buffer::Vector{T}, self::FIRFilter{FIRDecimator}, x::Vector{T} )
    h::Vector{T}       = self.kernel.h
    dlyLine::Vector{T} = self.dlyLine
    xLen               = length( x )
    dlyLineLen         = length( dlyLine )
    reqDlyLineLen      = self.reqDlyLineLen

    if xLen + dlyLineLen  < reqDlyLineLen
        append!( dlyLine, x )
        return T[]
    end

    decimation   = self.kernel.decimation
    xStartIdx    = reqDlyLineLen - dlyLineLen + 1
    xOffset      = xStartIdx - 1
    xLen         = length( x ) - xOffset
    hLen         = self.kernel.hLen
    outLen       = iceil(  xLen / decimation  )
    # buffer       = similar( x, outLen )
    criticalYidx = min( iceil( hLen / decimation ), outLen ) #
    dlyLine      = [ dlyLine, x[1:xStartIdx-1] ]

    inputIdx = 1                                               # inputIdx is is the current input, not the actual index of of the current x in the delay line
    for yIdx in 1:criticalYidx                                 # Filtering is broken up into two outer loops
        hIdx        = 1                                        # This first loop takes care of filter ramp up.
        accumulator = zero(T)

        for dlyLineIdx in inputIdx:reqDlyLineLen               # Inner Loop 1: Handles convolution of taps and delay line
            @inbounds accumulator += h[hIdx] * dlyLine[dlyLineIdx]
            hIdx += 1
        end

        for k in 1:inputIdx                                    # Inner Loop 2: handles convolution of taps and x
            @inbounds accumulator += h[ hIdx ] * x[ k + xOffset ]
            hIdx += 1
        end

        @inbounds buffer[yIdx] = accumulator
        inputIdx += decimation
    end

    inputIdx -= hLen

    for yIdx in criticalYidx+1:outLen                          # second outer loop, we are now in the clear to to convolve without going out of bounds
        accumulator  = zero(T)

        for k in 1:hLen
            @inbounds accumulator += h[ k ] * x[ inputIdx + k + xOffset ]
        end

        @inbounds buffer[yIdx] = accumulator
        inputIdx += decimation
    end
    inputIdx += hLen - decimation

    xLeftoverLen = max( xLen - inputIdx, 0 )
    dlyLineLen   = reqDlyLineLen - decimation + xLeftoverLen + 1
    self.dlyLine = [ dlyLine, x[xStartIdx:end]][end-dlyLineLen+1:end]

    return buffer
end

function filt{T}( self::FIRFilter{FIRDecimator}, x::Vector{T} )
    xLen       = length( x ) - self.reqDlyLineLen + length( self.dlyLine )
    outLen     = iceil(  xLen / self.kernel.decimation  )
    buffer     = similar( x, outLen )
    filt!( buffer, self, x )
end


end # module
