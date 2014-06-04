import Winston
import DSP

#==============================================================================#
#                 Plot Filter's Impulse & Frequency Response                   #
#==============================================================================#
# x = vector containing a filter's impulse response
# 
# Example:
#    plot_response( firdes(0.5, kaiser(37, 5.653) ))
# See example 7.8 in DTSP
    
function freqz( coefficients::Vector )
    if !method_exists( Winston.plot, ())
        error( "To use plot_response, you must load the Winston package")
    end
    x  = coefficients
    xx = [x, zeros(1024 - length(x))]
    M  = [0:length(x) - 1]
    n  = int(ceil(length( x )/2))
    X  = fft(xx)[1:512]
    X  = abs(X)
    X  = 20*log10( X )
    f  = linspace( 0, 1, 512 )
    
    impulse = Winston.stem( M,
                            x,
                            title  = "Impulse Response",
                            xlabel = "Sample #",
                            ylabel = "Amplitude"
                          )
    Winston.setattr(impulse.frame1, draw_grid=true)
    
    freq = Winston.FramedPlot( 
                            title  = "Frequency Response",
                            xlabel = "f/f_{Nyquist}",
                            ylabel = "dB"
                        )
    Winston.add(freq, Winston.Curve(f, X))
    Winston.setattr(freq.frame1, draw_grid=true)                    
    
    t = Winston.Table(2, 1)
    t[1,1] = impulse
    t[2,1] = freq

    return t
end

function freqz( FIR::FIRFilter )
   freqz( FIR.coefficients ) 
end


#==============================================================================#
#                               Plot Constellation                             #
#==============================================================================#

function plot_constellation( symbols::Vector )
    p = Winston.scatter( symbols,
                         ".",
                         xlabel = "In Phase",
                         ylabel = "Quadrature"
                       )

    return p
end



#==============================================================================#
#                                Plot Spectrum                                 #
#==============================================================================#

function plot_spectrum( signal::Vector, sample_rate::Real = 1.0 )
    spectrum = DSP.periodogram( signal )
    spectrum = 10*log10( spectrum.^2 ) 
    spectrum = fftshift( spectrum )
    
    frequencies = linspace( -sample_rate/2, sample_rate/2, length(spectrum) )
    
    spectrum_plot = Winston.FramedPlot( 
                            title  = "Spectrum",
                            xlabel = "frequency",
                            ylabel = "dB"
                        )
    Winston.add(spectrum_plot, Winston.Curve( frequencies, spectrum ))
    
    return spectrum_plot
end