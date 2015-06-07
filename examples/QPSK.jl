using Radio, Winston

# generate 10,000 random QPSK symbols
symbols = pskmod( 10000, 4 )
# create some gaussian noise and add it to the symbols
noise  = wgn( length( symbols ), 10, "dBm", 1.0, true )
signal = symbols .+ noise

constellation = plot_constellation( signal )
setattr( constellation, title = "QPSK Modulation" )

display( constellation )
file( constellation, "QPSK.png" )