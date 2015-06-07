pll   = PLL()
ƒ     = -0.1
ϕ     = 3.0
n     = 300
t     = 0:n-1
x     = exp( im*( ƒ*t + ϕ) )
(y,e) = exec( pll, x )


using Winston
plotTable = Table(3,1)
plotTable[1,1] = plot( t, real(x), t, real(y), "r", title = "Real" )
plotTable[2,1] = plot( t, imag(x), t, imag(y), "r", title = "Imaginary" )
plotTable[3,1] = plot( t, e, title = "ϕ error" )
display( plotTable )
