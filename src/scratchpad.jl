import Radio
import DSP

M          = 8
modulation = Radio.PSK(M)
modem      = Radio.Modem( modulation, samplesPerSymbol=2 )
txData     = rand( 0:M-1, 1000 )
symbols    = Radio.modulate( modem, txData )
rxData     = Radio.demodulate( modem, symbols )

spectrum =  DSP.periodogram( symbols )



using Winston
p = plot(fftshift(log10(abs2(fft(symbols)))))
display( p )
