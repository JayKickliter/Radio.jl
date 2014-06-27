import Radio
import Winston

mod_order    = 4
interp_ratio = 2
symbol_count = 1000
plot_filename = "CZT Resampling.png"

# generate random data, make it non-square so it's easier to see of we downsampled correctly
data = rand( 0:mod_order-2, symbol_count )

# generate 1e3 random QPSK symbols
symbols = Radio.pskmod( data, 4 )

# generate root-cosine nyquist filter coefficients
nyq_coef = Radio.rcos(0.35, 4 .* interp_ratio, interp_ratio)
nyq_coef = nyq_coef ./ sum( nyq_coef )

# interpolate symbols using the nyquest filter we just created
symbols_interp = Radio.interpolate( symbols, interp_ratio, nyq_coef )

# create some gaussian noise and add it to interpolated symbols
noise  = Radio.wgn( length( symbols_interp ), 0, "dBm", 1.0, true )
signal = symbols_interp .+ noise

constellation_base = Radio.plot_constellation( symbols )
Winston.setattr( constellation_base, title = "Baseband QPSK" )

spectrum_base = Radio.plot_spectrum( symbols, 1.0 )
Winston.setattr( spectrum_base, title = "Baseband Spectrum" )

constellation_interp = Radio.plot_constellation( signal )
Winston.setattr( constellation_interp, title = @sprintf("L=%d Interpolated QPSK Constellation", interp_ratio ) )

spectrum_interp = Radio.plot_spectrum( signal, interp_ratio )
Winston.setattr( spectrum_interp, title = @sprintf("L=%d Interpolated QPSK Spectrum", interp_ratio) )


# For spectral zoom, set:
#   W = exp(-j*2*pi*(f2-f1)/(m*fs));
#   A = exp(j*2*pi*f1/fs);
# where:
#   f1 = start freq
#   f2 = end freq
#   m  = length of x
#   fs = sampling rate

m  = symbol_count
f1 = -0.5
f2 =  0.5
fs = interp_ratio
W  = exp(-im*2*pi*(f2-f1)/(m*fs))
A  = exp(im*2*pi*f1/fs)

spectrum_czt = Radio.czt( signal, m, W, A )
# spectrum_czt = Radio.czt( signal, length(signal)) 
signal_decim = ifft( fftshift(spectrum_czt) )[interp_ratio*4:end]


constellation_decim = Radio.plot_constellation( signal_decim )
Winston.setattr( constellation_decim, title = "Chrip-z Decimated QPSK Constellation" )

spectrum_decim = Radio.plot_spectrum( signal_decim, 1.0 )
Winston.setattr( spectrum_decim, title = "Chrip-z Decimated QPSK Spectrum" )

t = Winston.Table( 3, 2 )
t[1,1] = constellation_base
t[1,2] = spectrum_base
t[2,1] = constellation_interp
t[2,2] = spectrum_interp
t[3,1] = constellation_decim
t[3,2] = spectrum_decim

Winston.display( t )
Winston.file( t, joinpath( dirname(Base.source_path()), plot_filename ) )

