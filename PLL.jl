# parameters and simulation options
phase_in     = 3.00  # carrier phase offset
frequency_in = -0.20 # carrier frequency offset
alpha        = 0.002 # PLL bandwidth
n            = 400   # number of samples

# initialize states
beta          = sqrt(alpha) # phase adjustment factor
phase_out     = 0.0         # output signal phase
frequency_out = 0.0         # output signal frequency

x = Array(Complex128, n)
y = Array(Complex128, n)

# print line legend to standard output
# printf("# %6s %12s %12s %12s %12s %12s\n", "index", "real(in)", "imag(in)", "real(out)", "imag(out)", "error")

# run basic simulation
for i in 0:n-1
    # compute input and output signals
    signal_in  = exp(im * phase_in)
    signal_out = exp(im * phase_out)
    
    println(i)
    x[i+1] = signal_in
    y[i+1] = signal_out
    
    # compute phase error estimate
    phase_error = angle( signal_in * conj(signal_out) )

    # apply loop filter and correct output phase and frequency
    frequency_out += alpha * phase_error    # adjust phase
    phase_out     +=  beta * phase_error    # adjust frequency

    # increment input and output phase values
    phase_in  += frequency_in
    phase_out += frequency_out
end


using Winston
plot( t, real(x), t, real(y), "p")