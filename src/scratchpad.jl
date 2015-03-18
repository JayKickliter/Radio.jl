import Radio

M          = 8
modulation = Radio.PSK(M)
modem      = Radio.Modem( modulation )
data       = rand( 0:M-1, 1000 )
txSymbols  = Radio.mod( modem, data ) 
txSymbols  = Radio.mod( modem, data ) 

p = plot( real(txSymbols) )

p = plot(fftshift(abs2(fft(txSymbols))))
display( p )
