#==============================================================================#
#                           Plot Frequency Response                            #
#==============================================================================#
# x = vector containing a filter's impulse response
# 
# Example:
#    plot_response( firdes(0.25, kaiser(37, 5.653) ))
# See example 7.8 in DTSP
    
function plot_response( x )
    if !method_exists( plot, ())
        error( "To use plot_response, you must load the Winston package")
    end
    M = [0:length(x) - 1]
    n = int(ceil(length( x )/2))
    X = fft(x)[1:n]
    X = abs(X)
    X = 20*log10( X )
    f = linspace( 0, 1, n )
    
    impulse = FramedPlot( 
                            title  = "Impulse Response",
                            xlabel = "Sample #",
                            ylabel = "Amplitude"
                        )
    add(impulse, Curve(M, x))                    
    add(impulse, Points(M, x, kind="filled circle"))
    
    freq = FramedPlot( 
                            title  = "Frequency Response",
                            xlabel = "Normalized Frequency",
                            ylabel = "Decibels"
                        )
    add(freq, Curve(f, X))                    
    add(freq, Points(f, X, kind="filled circle"))
    
    t = Table(2, 1)
    t[1,1] = impulse
    t[2,1] = freq
    display(t)
end
