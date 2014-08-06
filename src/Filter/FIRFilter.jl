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
    state::Vector
end


function FIRFilter{Tx}( ::Type{Tx}, h::Vector )
    hLen     = length( h )
    state  = zeros( Tx, hLen-1 )
    kernel = FIRStandard( h )
    FIRFilter( kernel, state )
end

FIRFilter{Tt}( h::Vector{Tt} ) = FIRFilter( Tt, h )

function FIRFilter{Tx}( ::Type{Tx}, h::Vector, resampleRatio::Rational )

    interploation = num( resampleRatio )
    decimation    = den( resampleRatio )

    if resampleRatio == 1
        return FIRFilter( Tx, h )
    elseif interploation == 1
        PFB    = h
        kernel = FIRDecimator( PFB, decimation, Tx[] )
    elseif decimation == 1
        PFB    = polyize( h, interploation )
        kernel = FIRInterpolator( PFB, interploation, Tx[] )
    else
        PFB    = polyize( h, interploation )
        kernel = FIRRational( PFB, interploation, decimation, Tx[] )
    end

    hLen  = size(PFB)[1]
    state = zeros( Tx, hLen - 1 )

    FIRFilter( kernel, state )
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

function filt!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, state::Vector{T} = T[] )

    bufLen      = length( buffer )
    xLen        = length( x )
    hLen        = length( h )
    stateLen    = length( state )
    reqStateLen = hLen - 1

    bufLen >= xLen || error( "buffer length must be >= x length")

    if stateLen != reqStateLen          # TODO: write the filtering logic to not depends on state being a certain length, as the current implementation allocates useless zeros
        if stateLen == 0
            state = zeros( T, reqStateLen )
        elseif stateLen < reqStateLen
            state = [ zeros( T, reqStateLen ), state ]
        else
            state = state[ end+1-reqStateLen:end ]
        end
    end

    h = flipud( h )                     # flip the h to make the multiplication more SIMD friendly

    for bufIdx in 1:hLen-1              # this first loop takes care of filter ramp up and previous state

        accumulator = zero(T)
        hIdx        = 1

        for stateIdx in bufIdx:stateLen # this loop takes care of previous state
            @inbounds accumulator += h[hIdx] * state[stateIdx]
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


filt{T}( h::Vector{T}, x::Vector{T}, state::Vector{T} = T[] ) = filt!( similar(x), h, x, state )

function filt( self::FIRFilter{FIRStandard}, x )
    h       = self.kernel.h
    hLen       = length( h )
    nextState  = x[end-hLen+2:end]
    y          = filt( self.kernel.h, x, self.state )
    self.state = nextState

    return y
end


function short_singlerate_test( h, x, state )
    @printf( "Radio's Single-rate filt\n\t")
    @time nativeResult = filt( h, x, state )
    @printf( "Base Single-rate filt\n\t")
    @time baseResult   = Base.filt( h, 1.0, x )

    self = FIRFilter( h )
    
    @printf( "Stateful Single-rate filt\n\t")
    @time y = [ filt( self, x[1:250] ) , filt( self, x[251:end] ) ]

    # [ baseResult nativeResult y  ]
    areApprox( nativeResult, baseResult ) && areApprox( y, baseResult )    
end

#=============================
h = rand( 56 )
x = rand( 1_000_000 )
state = zeros( length(h) - 1 )

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

function decimate!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, decimation::Integer, state::Vector{T} = T[] )
    xLen   = length( x )
    hLen   = length( h )
    outLen = floor(int(xLen / decimation))

    length( buffer ) * decimation >= xLen || error( "buffer lenght must be >= signal length * decimation" )

    criticalYidx = int(ceil(hLen / decimation)) # The index of y where our h would overlap
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

function short_decimate_test( h, x, factor )
    @printf( "Radio's decimation\n\t")
    @time nativeResult = decimate( h, x, factor );
    
    @printf( "Naive resampling\n\t")
    @time begin
        baseResult   = Base.filt( h, 1.0, x );
        baseResult   = baseResult[1:factor:end];
    end
    areApprox( nativeResult, baseResult )
end

#===================================
h      = rand( 128 );
x      = rand( 10_000_000 );
factor = 100
short_decimate_test( h, x, factor )
===================================#