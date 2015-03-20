import Radio
import DSP

M          = 4
modulation = Radio.PSK(M)
modem      = Radio.Modem( modulation, samplesPerSymbol=2 )
txData     = rand( 0:M-1, 1000 )
symbols    = Radio.modulate( modem, txData )
rxData     = Radio.demodulate( modem, symbols )
taps       = modem.txFilter.h

p = DSP.welch_pgram( symbols, window=DSP.hanning )
freqs = fftshift(DSP.freq( p ))
power = log10( fftshift(DSP.power( p )) )


using PyPlot

clf()
subplot( 2,1,1 )
scatter( real(symbols), imag(symbols) )

subplot( 2,1,2 )
plot( freqs, power )
