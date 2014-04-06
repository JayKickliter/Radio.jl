using Radio, Winston

symbols       = pskmod( 10000, 4 )
noise         = wgn( length( symbols ), 10, "dBm", 1.0, true )
signal        = symbols .+ noise

constellation = scatter( real(signal), imag(signal), "." )
setattr( constellation,
            title = "QPSK Constellation",
            xlabel = "In Phase",
            ylabel = "Quadrature"
)                

display( constellation )
file( constellation,  "QPSK.png" )
