ƒ = -0.2
ϕ = 0
n = 400
t = 0:n-1
x = exp(im*( ƒ*t + ϕ) )
y = exec( pll, x)

pll = PLL()

using Winston
plotTable = Table(3,1)
plotTable[1,1] = plot( t, real(x), t, real(y), "b", title = "Real" )
plotTable[2,1] = plot( t, imag(x), t, imag(y), "b", title = "Imaginary" )
display( plotTable )