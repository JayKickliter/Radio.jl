function Base.LinAlg.dot{T<:Union(Float32, Float64)}( n::Integer, x::Vector{T<:Union(Float32, Float64)}, xStart::Integer, y::Vector{T<:Union(Float32, Float64)}, yStart::Integer )
    xLen  = length( x )
    yLen  = length( y )    
    xxLen = xLen - xStart + n
    yyLen = yLen - yStart + n
    
    xPtr = Ptr( x ) + ( xStart - 1 ) * sizeof(T)
    yPtr = Ptr( y ) + ( yStart - 1 ) * sizeof(T) 
    
    Base.BLAS.dot( n, xPtr, 1, yPtr, 1 )
end

# dot(100, x, 150, y, 75 )