# Reference: Chirp Z-Transform Spectral Zoom Optimization with MATLAB
# http://prod.sandia.gov/techlib/access-control.cgi/2005/057084.pdf

# For spectral zoom, set:
#   W = exp(-j*2*pi*(f2-f1)/(m*fs));
#   A = exp(j*2*pi*f1/fs);
# where:
#   f1 = start freq
#   f2 = end freq
#   m  = length of x
#   fs = sampling rate
function czt(x::Vector{Complex128}, m::Int, w::Complex128, a::Complex128)
    n = length(x)
    N = [0:n-1]+n
    NM = [-(n-1):(m-1)]+n
    M = [0:m-1]+n

    nfft = nextpow2(n+m-1)
    W2 = w.^(([-(n-1):max(m-1,n-1)].^2)/2)

    fg = zeros(Complex128, nfft)
    fg[1:n] = x.*(a.^-(N-n)).*W2[N]
    fg = fft(fg)

    fw = zeros(Complex128, nfft)
    fw[1:length(NM)] = 1./W2[NM]
    fw = fft(fw)
    gg = ifft(fg.*fw)

    return gg[M].*W2[M]
end

function czt(x::Vector{Complex128}, m::Int)
    czt( x, m, exp(-im*2*pi/m), 1.0+0.0im)
end