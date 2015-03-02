type PLL
    α # Bandwidth
    β # Phase adjustor
    ƒ # Output frequency
    ϕ # Output phase
end

function PLL( α = 0.002, β = sqrt( α ) )
    PLL( α, β, 0.0, 0.0 )
end

function exec( pll::PLL, x::Complex )
    y  = exp( pll.ϕ * im )
    Δϕ = angle( x * conj(y) )
    pll.ƒ += pll.α * Δϕ
    pll.ϕ += pll.β * Δϕ
    pll.ϕ += pll.ƒ
    
    return y
end

function exec{T}( pll::PLL, x::AbstractVector{T} )
    T[ exec( pll, x ) for x in x ]
end
