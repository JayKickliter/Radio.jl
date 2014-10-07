#==============================================================================#
#                                  Upsample                                    #
#==============================================================================#
# length        = desired length of WGN noise vector
# power         = desired power level
# units         = "unitless", "watts", "dBm", "dBW". Case insensitive.
# impedence     = resistance of the load
# returnComplex = true/false. Whether or not to return complex noise 
# 
# Returned vector is in volts. The function uses paramter 'units' and 'impedence' to determine the RMS voltage level
function wgn(length::Integer; power::Real=1.0, units::String = "unitless", impedence::Real=1.0, returnComplex::Bool=false)
    # TODO: add argument valdiation
    units = lowercase( units )
    if units == "unitless"
        np = power^2/impedence
    elseif units == "watts"
        np = power
    elseif units == "dbm"
        np = 10^((power - 30)/10)
    elseif units == "dbw"
        np = 10^(power/10)
    end

    if returnComplex
        noiseVec = sqrt(impedence*np) * ( randn(length) + randn(length)im )
    else
        noiseVec = sqrt(impedence*np) * randn(length)
    end

    return noiseVec
end