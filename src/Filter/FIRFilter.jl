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
    interploation::Int
    xLeftover::Vector
end

# Decimator FIR kernel
type FIRDecimator <: FIRKernel
    h::Vector
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
    dlyLine::Vector
end


function FIRFilter{Tx}( ::Type{Tx}, h::Vector )
    hLen     = length( h )
    dlyLine  = zeros( Tx, hLen-1 )
    kernel = FIRStandard( h )
    FIRFilter( kernel, dlyLine )
end

FIRFilter{Tt}( h::Vector{Tt} ) = FIRFilter( Tt, h )

function FIRFilter{Tx}( ::Type{Tx}, h::Vector, resampleRatio::Rational )

    interploation = num( resampleRatio )
    decimation    = den( resampleRatio )

    if resampleRatio == 1
        return FIRFilter( Tx, h )
    elseif interploation == 1
        PFB    = h
        kernel = FIRDecimator( PFB, decimation, zeros( Tx, decimation-1 ) )
    elseif decimation == 1
        PFB    = polyize( h, interploation )
        kernel = FIRInterpolator( PFB, interploation, Tx[] )
    else
        PFB    = polyize( h, interploation )
        kernel = FIRRational( PFB, interploation, decimation, Tx[] )
    end

    hLen  = size(PFB)[1]
    dlyLine = zeros( Tx, hLen - 1 )

    FIRFilter( kernel, dlyLine )
end

FIRFilter{Tt}( h::Vector{Tt}, resampleRatio::Rational ) = FIRFilter( Tt, h, resampleRatio )




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

    bufLen      = length( buffer )
    xLen        = length( x )
    hLen        = length( h )
    dlyLineLen    = length( dlyLine )
    reqDlyLineLen = hLen - 1

    bufLen >= xLen || error( "buffer length must be >= x length")

    if dlyLineLen != reqDlyLineLen          # TODO: write the filtering logic to not depends on dlyLine being a certain length, as the current implementation allocates useless zeros
        if dlyLineLen == 0
            dlyLine = zeros( T, reqDlyLineLen )
        elseif dlyLineLen < reqDlyLineLen
            dlyLine = [ zeros( T, reqDlyLineLen ), dlyLine ]
        else
            dlyLine = dlyLine[ end+1-reqDlyLineLen:end ]
        end
    end

    h = flipud( h )                     # flip the h to make the multiplication more SIMD friendly

    for bufIdx in 1:hLen-1              # this first loop takes care of filter ramp up and previous dlyLine

        accumulator = zero(T)
        hIdx        = 1

        for dlyLineIdx in bufIdx:dlyLineLen # this loop takes care of previous dlyLine
            @inbounds accumulator += h[hIdx] * dlyLine[dlyLineIdx]
            hIdx += 1
        end

        hIdx = hLen-bufIdx+1
        for xIdx in 1:bufIdx            # this loop takes care of the first hlen-1 samples in x
            @inbounds accumulator += h[hIdx] * x[xIdx]
            hIdx += 1
        end

        @inbounds buffer[bufIdx] = accumulator
    end

    for bufIdx in hLen:xLen             # filter ramp is complete, normal filtering from this point on

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
    h          = self.kernel.h
    hLen       = length( h )
    nextDlyLine  = x[end-hLen+2:end]
    y          = filt( self.kernel.h, x, self.dlyLine )
    self.dlyLine = nextDlyLine

    return y
end


function short_singlerate_test( h, x, dlyLine )
    @printf( "Radio's Single-rate filt\n\t")
    @time nativeResult = filt( h, x, dlyLine )
    @printf( "Base Single-rate filt\n\t")
    @time baseResult   = Base.filt( h, 1.0, x )

    self = FIRFilter( h )
    
    @printf( "dlyLineful Single-rate filt\n\t")
    @time y = [ filt( self, x[1:250] ) , filt( self, x[251:end] ) ]

    # [ baseResult nativeResult y  ]
    areApprox( nativeResult, baseResult ) && areApprox( y, baseResult )    
end

#=============================
h = rand( 56 )
x = rand( 1_000_000 )
dlyLine = zeros( length(h) - 1 )

short_singlerate_test( h, x )
=============================#




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function interpolate!{T}( buffer::Vector{T}, PFB::Array{T, 2}, x::Vector{T} )
    (hLen, Nφ)  = size( PFB )      # each column is a phase of the PFB, the rows hold the individual h
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

        xStartIdx   = xIdx-hLen
        yIdx        = Nφ*(xIdx-1)+φ
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

function short_interpolate_test( h, x, factor )    
    @printf( "Radio's interpolate\n\t")
    @time nativeResult = interpolate( h, x, factor )
    @printf( "Naive interpolate\n\t")
    @time begin
        xx = zeros( length(x) * factor )
        for n = 0:length(x)-1;
            xx[ n*factor+1 ] = x[ n+1 ]
        end
        baseResult = Base.filt( h, 1.0, xx )        
    end
        
    areApprox( nativeResult, baseResult )
end

#=============================
h = rand( 56 )
x = rand( 1_000_000 )
factor = 4

short_interpolate_test( h, x, factor )
=============================#




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function resample!{T}( buffer::Vector{T}, PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )
    (hLen, Nφ)    = size( PFB ) # each column is a phase of the PFB, the rows hold the individual h
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

function short_rational_test( h, x, ratio )
    upfactor   = num( ratio )
    downfactor = den( ratio )
    
    @printf( "Radio's rational resampling\n\t")
    @time nativeResult = resample( h, x, upfactor//downfactor );
    
    @printf( "Naive resampling\n\t")
    @time begin
        xx         = zeros( length(x) * upfactor );
        baseResult = similar( x, int( length(x) * upfactor / downfactor ))
        
        for n = 0:length(x)-1;
            xx[ n*upfactor+1 ] = x[ n+1 ];
        end
        
        baseResultInterpolated = Base.filt( h, 1.0, xx );
        baseResult = [ baseResultInterpolated[n] for n = 1:downfactor:length( baseResultInterpolated ) ]
    end
    
    [ baseResult nativeResult ]
    areApprox( nativeResult, baseResult )
end

#=========================
ratio = 3//4
h     = rand( 56 )
x     = rand( 1000000 )
short_rational_test( h, x, ratio )
=========================#




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function decimate!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, decimation::Integer, dlyLine::Vector{T} = T[] )

    xLen          = length( x )
    hLen          = length( h )
    dlyLineLen    = length( dlyLine )
    reqDlyLineLen = hLen-1
    outLen        = int( ceil( xLen / decimation ))

    length( buffer ) >= outLen || error( "length(buffer) must be >= floor( length(x) / decimation)" )

    if dlyLineLen != reqDlyLineLen                                      # TODO: write the filtering logic to not depends on dlyLine being a certain length, as the current implementation allocates useless zeros
        println( "dlyLineLen, $dlyLineLen, is not the correct size")
        if dlyLineLen == 0
            dlyLine = zeros( T, reqDlyLineLen )
        elseif dlyLineLen < reqDlyLineLen
            dlyLine = [ zeros( T, reqDlyLineLen - dlyLineLen ), dlyLine ]
        else
            dlyLine = dlyLine[ end+1-reqDlyLineLen:end ]
        end
        dlyLineLen = length( dlyLine )
    end
    
    criticalYidx = int(ceil(hLen / decimation))                     # The maxximum index of y where our h*x would would rech out of bounds
    xIdx         = 1
    
    h = flipud(h)                                                   # TODO: figure out a way to not always have to flip taps each time

    for yIdx = 1:criticalYidx
        
        hIdx        = 1
        accumulator = zero(T)

        for dlyLineIdx = xIdx:dlyLineLen                                # this loop takes care of previous dlyLine
            @inbounds accumulator += h[hIdx] * dlyLine[dlyLineIdx]
            hIdx += 1
        end
        
        for k = 1:xIdx
            @inbounds accumulator += h[ hIdx ] * x[ k ]
            hIdx += 1
        end

        @inbounds buffer[yIdx] = accumulator
        xIdx += decimation
    end

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

decimate{T}( h::Vector{T}, x::Vector{T}, decimation::Integer, dlyLine::Vector{T} = T[] ) = decimate!( similar(x, int(ceil( length( x ) / decimation )) ), h, x, decimation, dlyLine )

function filt( self::FIRFilter{FIRDecimator}, x::Vector )
    xLen          = length( x )
    h             = self.kernel.h
    hLen          = length( h )
    reqDlyLineLen = hLen - 1
    dlyLine       = self.dlyLine
    dlyLineLen    = length( dlyLine )
    decimation    = self.kernel.decimation
    xLeftover     = self.kernel.xLeftover
    xLeftoverLen  = length( xLeftover )
    
    emptyResult = similar( x, 0 )
    
    # if there are no leftovers, the next x index an integer multiple of decimation
    if xLeftoverLen == 0
        nextLeftoverLen       = mod( xLen, decimation )
        nextLeftover          = x[end-nextLeftoverLen+1:end]
        self.kernel.xLeftover = nextLeftover
        return decimate( h, x, decimation, dlyLine )
    end
    
    if xLeftoverLen + xLen >= decimation
        dlyLine         = [dlyLine, xLeftover, x[1:decimation-1-xLeftoverLen]][1:reqDlyLineLen]
        self.dlyLine    = dlyLine
        xStart          = decimation-xLeftoverLen
        x               = x[xStart:end]
        xLen            = length( x )
        nextLeftoverLen = mod( xLen, decimation )
        nextLeftover    = x[end-nextLeftoverLen+1:end]
        self.kernel.xLeftover = nextLeftover
        return decimate( h, x, decimation, dlyLine )
    end
    
    append!( xLeftover, x )
    
    return emptyResult
end


function short_decimate_test( h, x, factor )
    @printf( "\nRadio's decimation\n\t")
    @time nativeResult = decimate( h, x, factor )
    display( nativeResult )
    
    self = FIRFilter( h, 1//factor )
    @printf( "\nDlyLineful decimation\n\t")
    @time begin
        y1   = filt( self, x[1:8])
        y2   = filt( self, x[9:100] )
    end
    dlyLinefulResult = [y1, y2]
    display( dlyLinefulResult )
    
    @printf( "\nNaive resampling\n\t")
    @time begin
        baseResult = Base.filt( h, 1.0, x )[1:factor:end]
    end
    display( baseResult )
    
    areApprox( dlyLinefulResult, baseResult ) & areApprox( nativeResult, baseResult ) ? println( "Tests passed" ) : println( "1 or more tests failed")
end

#===================================

h      = ones(10)./10;
x      = [1.0:100];
factor = 6
self   = FIRFilter( h, 1//factor )

y = similar(x, 0)
for i = 1:length(x)
    yNext = filt( self, x[i:i] )
    display( yNext )
    y = [y, yNext]
end
y


short_decimate_test( h, x, factor )

===================================#
