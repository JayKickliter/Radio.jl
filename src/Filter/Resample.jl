#==============================================================================#
#                                  Upsample                                    #
#==============================================================================#
function upsample( x::Vector, L::Int )
    # TODO: add argument valdiation
    M = length( x )
    K = M*L
    
    h = similarzeros( x, K )
    
    for i = 1:M
       h[i*L-L+1] = x[i] 
    end
    
    return h
end

#==============================================================================#
#                                  Interpolate                                 #
#==============================================================================#
function interpolate( x::Vector, L::Int, Filter::Vector )
    # TODO: add argument valdiation
    M = length( x )
    K = M*L
    
    h = similarzeros( x, K )
    
    for i = 1:M
       h[i*L-L+1] = x[i] 
    end
    
    h = filt( complex128( Filter ), complex128(1.0), h )
end