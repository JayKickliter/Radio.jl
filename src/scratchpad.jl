using Winston
using Radio
using DSP

M          = 8
modulation = Radio.PSK(M)
modem      = Radio.Modem( modulation, samplesPerSymbol=2 )
data       = rand( 0:M-1, 1000 )
txSymbols  = Radio.modulate( modem, data ) 

spectrum =  DSP.periodogram( txSymbols )

p = plot(fftshift(log10(abs2(fft(txSymbols)))))
display( p )
