type PLL
    α # Bandwidth
    β # Phase adjustor
    ƒ # Output frequency
    ϕ # Output phase
end

function PLL( α = 0.002, β = sqrt( α ) )
    PLL( α, β, 0.0, 0.0 )
end

function exec( pll::PLL, x::Number )
    y  = exp( pll.ϕ * im )
    Δϕ = angle( x * conj(y) )
    pll.ƒ += pll.α * Δϕ
    pll.ϕ += pll.β * Δϕ
    pll.ϕ += pll.ƒ
    
    return y, Δϕ
end

function exec{T}( pll::PLL, x::AbstractVector{T} )
    y = Array( T, length(x) )
    e = Array( typeof(real(x[1])), length(x) )

    for i in 1:length( x )
        ( yi, ei ) = exec( pll, x[i] )
        y[i] = yi
        e[i] = ei
    end
    
    return y, e
end

