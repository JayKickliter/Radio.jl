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
    (Ntaps, Nφ) = size( PFB )           # each column is a phase of the PFB, the rows hold the individual taps
    Nx          = length( x )           # number of input items
    Ny          = Nx * Nφ
    y           = similar( x, Nx * Nφ ) # Ny = Nx * Nφ
    
    for Xn = 1:Ntaps-1, φ = 1:Nφ        # until Xn == Ntaps, x[Xn-m+1] would reach out of bounds 
                                        # this first loop limits Xn-m+1 to a minimum of 1
        Yn          = Nφ*(Xn-1)+φ
        accumulator = zero(T)
        for Tn = 1:Xn                   # for each tap in phase[n]
            @inbounds accumulator += PFB[Tn, φ] * x[Xn-Tn+1]
        end
        y[Yn] = accumulator
    end
    
    PFB = flipud(PFB)
    
    for Xn = Ntaps:Nx, φ = 1:Nφ         # no longer in danger of stepping out of bounds
        
        XnBase      = Xn-Ntaps               
        Yn          = Nφ*(Xn-1)+φ
        accumulator = zero(T)        
                
        @simd for Tn = 1:Ntaps
            @inbounds accumulator += PFB[Tn, φ] * x[XnBase + Tn]
        end
        
        y[Yn] = accumulator
    end

    return y
end

interpolate( h, x, interpolation ) = polyinterpolate( polyize( h, interpolation ), x )




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function resample{T}( PFB::Array{T, 2}, x::Vector{T}, ratio::Rational )
    
    (Ntaps, Nφ) = size( PFB ) # each column is a phase of the PFB, the rows hold the individual taps
    L    = num(ratio)
    M    = den( ratio )    
    xLen = length( x )        # number of input items    
    
    xLen * L % M == 0 || error("signal length * interpolation mod decimation must be 0")    
    
    yLen = int(xLen*L/M)
    y    = zeros( T, yLen )
        
    for m = 0:yLen-1
        
        φ    = mod( m*M, L)
        nm   = int( floor( m*M / L ))
        kMax = nm < Ntaps ? nm+1 : Ntaps
        acc  = zero(T)

        # @printf( "\n\nφ = %d", φ+1 )
        # @printf( "\tnm = %d", nm+1 )
        # @printf( "\tkMax = %d\n", kMax)
        # @printf( "y[%d] =", m+1 )
        
        
        for k = 0:kMax-1
            # @printf( "\tPFB[%d, %d] * x[%d]", k+1, φ+1, nm+1-k )
            # @printf( "\t%f * %f", PFB[k+1, φ+1], x[nm+1-k] )
            acc += PFB[ k+1, φ+1 ] * x[ nm+1-k ]
        end
        
        # @printf( "\ny[%d] = %f", m+1, acc )
        
        y[m+1] = acc
    
    end

    return y
end

resample{T}( h::Vector{T}, x::Vector{T}, ratio::Rational ) = resample( polyize(h, num(ratio)), x, ratio )



#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function decimate{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, M::Integer )
    xLen   = length( x )    
    hLen   = length( h )
    bufLen = length( bufLen )
    yLen   = int(xLen / M)    
    
    xLen   % M == 0     || error( "signal length % decimation must be 0" )
    bufLen * M >= xLen  || error( "buffer lenght must be >= signal length * decimation" )
    b
    
    y    = Array( T, yLen ) # yLen = xLen * Nφ

    criticalYidx = int(ceil(hLen / M))
        
    xIdx = 1
    
    for yIdx = 1:criticalYidx
            
        acc  = zero(T)
        
        kMax = xIdx < hLen ? xIdx : hLen
                
        for k = 1:kMax
            @inbounds acc += h[ k ] * x[ xIdx+1-k ]
        end
        
        @inbounds y[yIdx] = acc
        xIdx += M
    end
    
    h = flipud(h)
    
    xIdx -= hLen
    
    for yIdx = criticalYidx+1:yLen
            
        acc  = zero(T)
                
        for k = 1:hLen
            @inbounds acc += h[ k ] * x[ xIdx+k ]
        end
        
        @inbounds y[yIdx] = acc
        xIdx += M
    end
        
    return y
end