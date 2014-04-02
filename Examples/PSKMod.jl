using Radio
using Winston
import DSP: welch_pgram

symbols  = pskmod( 100000, 4, 4 )
noise    = wgn(length(symbols), 5, "dBm", 50, true)
signal   = symbols .+ noise
spectrum = welch_pgram( signal, 100, 50 )
spectrum = fftshift( spectrum )
spectrum = 10*log10( spectrum )

@printf( "Oversampled symbol RMS: %f", rms(symbols) )


display(plot( spectrum ))