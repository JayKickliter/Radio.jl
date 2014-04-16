#==============================================================================#
#                                Blackman Window                               #
#==============================================================================#

function blackman( n::Integer, N::Integer )
    # TODO: add argument valdiation
    α = 0.42
    β = 0.5
    α - β*cos(2*π*n/N) + (β-α)*cos(4*π*n/N)
end

function blackman( N::Integer )
    # TODO: add argument valdiation
    [ blackman(n, N) for n = 0:N ]
end

#==============================================================================#
#                                Hamming Window                                #
#==============================================================================#

function hamming( n::Integer, N::Integer )
    # TODO: add argument valdiation
    α = 0.54
    β = 0.46 # 1 - α
    α - β*cos(2*π*n/N)
end

function hamming( N::Integer )
    # TODO: add argument valdiation
    [ hamming(n, N) for n = 0:N ]
end

#==============================================================================# 
#                                  Hann Window                                 #
#==============================================================================#

function hann( n::Integer, N::Integer )
    # TODO: add argument valdiation
    α = 0.5
    β = 0.5
    α - β*cos(2*π*n/N)
end

function hann( N::Integer )
    # TODO: add argument valdiation
    [ hann(n, N) for n = 0:N ]
end

#==============================================================================#
#                                Kaiser Window                                 #
#==============================================================================#

function kaiser( n::Integer, N::Integer, β::Real )
    # TODO: add argument valdiation
    num = besseli(0, β*sqrt(1-(2*n/N-1).^2))
    dem = besseli(0, β)
    num/dem
end

function kaiser( N::Integer, β::Real )
    # TODO: add argument valdiation
    [ kaiser( n, N, β ) for n = 0:N ]
end

#==============================================================================#
#                              Rectangular Window                              #
#==============================================================================#

function rectangle( n::Integer, N::Integer )
    # TODO: add argument valdiation
    1.0
end

function rectangle( N::Integer )
    # TODO: add argument valdiation
    ones( typeof(1.0), N )
end