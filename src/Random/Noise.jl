function wgn(length::Integer, power::Real=1.0, impedence::Real=1.0, units::String = "linear", returnComplex::Bool=false)
    # TODO: add argument valdiation
    if units == "linear"
        np = power
    elseif units == "dBm"
        np = 10^((power - 30)/10)
    elseif units == "dBW"
        np = 10^(power/10)
    end

    if returnComplex
        noiseVec = sqrt(impedence*np/2) * ( randn(length) + randn(length)im )
    else
        noiseVec = sqrt(impedence*np/2) * randn(length)
    end

    return noiseVec
end