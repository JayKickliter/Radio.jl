#==============================================================================#
#                                Blackman Window                               #
#==============================================================================#

function blackman( n::Integer, M::Integer )
    # TODO: add argument valdiation
    α = 0.42
    β = 0.5
    α - β*cos(2*π*n/M) + (β-α)*cos(4*π*n/M)
end

function blackman( M::Integer )
    # TODO: add argument valdiation
    [ blackman(n, M) for n = 0:M ]
end

#==============================================================================#
#                                Hamming Window                                #
#==============================================================================#

function hamming( n::Integer, M::Integer )
    # TODO: add argument valdiation
    α = 0.54
    β = 0.46 # 1 - α
    α - β*cos(2*π*n/M)
end

function hamming( M::Integer )
    # TODO: add argument valdiation
    [ hamming(n, M) for n = 0:M ]
end

#==============================================================================# 
#                                  Hann Window                                 #
#==============================================================================#

function hann( n::Integer, M::Integer )
    # TODO: add argument valdiation
    α = 0.5
    β = 0.5
    α - β*cos(2*π*n/M)
end

function hann( M::Integer )
    # TODO: add argument valdiation
    [ hann(n, M) for n = 0:M ]
end

#==============================================================================#
#                                Kaiser Window                                 #
#==============================================================================#

function kaiser( n::Integer, M::Integer, β::Real )
    # TODO: add argument valdiation
    num = besseli(0, β*sqrt(1-(2*n/M-1).^2))
    dem = besseli(0, β)
    num/dem
end

function kaiser( M::Integer, β::Real )
    # TODO: add argument valdiation
    [ kaiser( n, M, β ) for n = 0:M ]
end

#==============================================================================#
#                              Rectangular Window                              #
#==============================================================================#

function rectangle( n::Integer, M::Integer )
    # TODO: add argument valdiation
    1.0
end

function rectangle( M::Integer )
    # TODO: add argument valdiation
    ones( typeof(1.0), M )
end