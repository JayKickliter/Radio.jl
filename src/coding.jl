abstract CodingScheme



################################################################################
#               ____ _  _ ____ ____ ___  _ _  _ ____                           #
#               |___ |\ | |    |  | |  \ | |\ | | __                           #
#               |___ | \| |___ |__| |__/ | | \| |__]                           #
################################################################################                                     
                                     
type Gray <: CodingScheme end


function encode( ::Type{Gray}, n::Integer )
    n $ (n >> 1)
end


function decode( ::Type{Gray}, n::Integer )
    p = n
     while (n >>= 1) != 0
         p $= n
     end
     return p
end


encode{T<:CodingScheme}( ::Type{T}, N::AbstractVector ) = [ encode( T, n) for n in N ]
decode{T<:CodingScheme}( ::Type{T}, N::AbstractVector ) = [ decode( T, n) for n in N ]