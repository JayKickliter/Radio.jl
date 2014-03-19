export firdes

#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# M      = Filter order: number of impulse response samples to return
# F_t    = Normalized transition frequency. CutoffFreq/SampleRate
# window = Window function to apply to ideal truncated sinc response
function firdes( M::Integer, F_t::Real, windowFunction::Function )
    if F_t < 0 && F_t > 1
        error("F_t must be greater than 0 and less than 1")
    end
    if !method_exists( windowFunction, ())
        error("Specified window doesn't exist. See Window.jl for valid windows")
    end
    
    H_t = windowFunction( M )
    
    for n in 1:M
        H_t[n] = sinc(2*π*F_t*(n-M/2)) * H_t[n]
    end
    
    return H_t
end