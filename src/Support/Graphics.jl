#==============================================================================#
#                           Plot Frequency Response                            #
#==============================================================================#
# x = vector containing a filter's impulse response

function plot_response( x )
    if !method_exists( plot, ())
        error( "To use plot_response, you must load the Winston package")
    end
    
    X = fft(x)
    X = abs(X)
    X = 20*log10( X )
    f = linspace( -0.5, 0.5, length(X) )
    display( plot(f, X) )
end
