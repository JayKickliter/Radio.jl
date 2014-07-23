function areApprox( x1::Vector, x2::Vector )
    Nx1 = length( x1 )
    Nx2 = length( x2 )    
    
    if Nx1 != Nx2
        println("x1 & x2 are different lengths vectors")
        return false
    end
    
    for i = 1:Nx1
        if !isapprox( x1[i], x2[i] )
            println( "Something went wrong at index $i" )
            return false
        end
    end
    
    return true
end

include("filters.jl")