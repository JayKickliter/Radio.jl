export plot_response

#==============================================================================#
#                           Plot Frequency Response                            #
#==============================================================================#
# x = vector containing a filter's impulse response

function plot_response( x )
    if !method_exists( plot, ())
        error( "To use plot_response, you must load the Winston package")
    end
    
    truncatedLength = length(x) >> 1
    X = fft(x)
    X = X[1:truncatedLength]
    X = 20*log10( abs(X) )
    f = linspace( 0.0, 1.0, truncatedLength )    
    display( plot(f, X) )
end
