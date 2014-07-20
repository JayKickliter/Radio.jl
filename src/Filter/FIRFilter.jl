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