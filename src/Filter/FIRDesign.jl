module Filter

export firdes

#==============================================================================#
#                                Windowed Filter                               #
#==============================================================================#
# M      = Filter order; Returns M+1 samples
# F_t    = Normalized transition frequency. CutoffFreq/SampleRate
# window = Window function to apply to ideal truncated sinc response
function firdes( M::Integer, F_t::FloatingPoint, windowFunction::Function )
    if F_t < 0 && F_t > 0.5
        error("F_t must be > 0.0 and < 0.5")
    end
    
    [ 2*F_t*sinc(2*F_t*(n-M/2)) * windowFunction(n+1, M+1) for n = 0:M ]
end

end # module Filter