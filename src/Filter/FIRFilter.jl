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

