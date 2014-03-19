export firdes

#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# M      = Filter order; Returns M+1 samples
# F_t    = Normalized transition frequency. CutoffFreq/SampleRate
# window = Window function to apply to ideal truncated sinc response
function firdes( M::Integer, F_t::Real, windowFunction::Function )
    if F_t < 0 && F_t > 1
        error("F_t must be greater than 0 and less than 1")
    end
    
    H_t = windowFunction( M+1 )
    
    for n in 0:M
        H_t[n+1] = 2*F_t*sinc(2*F_t*(n-M/2)) * H_t[n+1]
    end
    
    return H_t
end