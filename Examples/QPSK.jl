using Radio
using Winston

symbols       = pskmod( 10000, 4 )
noise         = wgn( length( symbols ), 10, "dBm", 1.0, true )
signal        = symbols .+ noise
constellation = scatter( real(signal), imag(signal), "." )

display( constellation )
file( constellation,  "QPSK.png" )
