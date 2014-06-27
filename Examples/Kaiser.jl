import Radio
const rad = Radio

( N, β )        = rad.kaiserord( 0.001, 0.2*π )
window          = rad.kaiser( N, β )
impulseResponse = rad.firdes( 0.5, window )
p               = rad.freqz( impulseResponse )

display( p )
file( p, "Kaiser.png" )