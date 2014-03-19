export blackaman, hamming, hann, kaiser, rectangular

#==============================================================================#
#                                Blackman Window                               #
#==============================================================================#

function blackaman( N::Integer )
    n = 0:N-1
    α = 0.42
    β = 0.5
    α - β*cos(2*π*n/(N-1)) + (β-α)*cos(4*π*n/(N-1))
end

#==============================================================================#
#                                Hamming Window                                #
#==============================================================================#

function hamming( N::Integer )
    n = 0:N-1
    α = 0.54
    β = 0.46 # 1 - α
    α - β*cos(2*π*n/(N-1))
end

#==============================================================================# 
#                                  Hann Window                                 #
#==============================================================================#

function hann( N::Integer )
    n = 0:N-1
    α = 0.5
    β = 0.5
    α - β*cos(2*π*n/(N-1))
end

#==============================================================================#
#                                Kaiser Window                                 #
#==============================================================================#

function kaiser( N::Integer, α::Real )
    β   = π*α
    n   = 0:N-1
    num = besseli(0, β*sqrt(1-(2*n/(N-1)-1).^2))
    dem = besseli(0, β)
    num/dem
end

#==============================================================================#
#                              Rectangular Window                              #
#==============================================================================#

function rectangular( N::Integer )
    ones( Float64, N )
end