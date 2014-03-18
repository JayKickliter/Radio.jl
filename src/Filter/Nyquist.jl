#==============================================================================#
#                                Raised Cosine                                 #
#==============================================================================#
# α    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# sps  = Samples per symbol

function rcos( α, span, sps )
    hlen = sps*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h    = Array( Float64, hlen )
    
    for i = 1:hlen
        t = i - 1 - (hlen - 1)/2
        if abs(α*t/sps) == 0.5      # check for case where (1-(2*α*t/sps)^2 is zero, which would cause divide by zero
            h[i] = 0.5
        else
            h[i] = sinc(t/sps)*cos(π*α*t/sps)/(1-(2*α*t/sps)^2)
        end
    end 

    return h
end

#==============================================================================#
#                             Root Raised Cosine                               #
#==============================================================================#
# α    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# sps  = Samples per symbol

function rrcos( α, span, sps )
    hlen = sps*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h    = Array( Float64, hlen )
    
    for i = 1:hlen
        t = i - (hlen - 1)/2 - 1 
        if t == 0                   # full equation is undefined at t = 0
            h[i] = (1/sqrt(sps)) * (1+α*(4/π -1))  
        elseif abs(t) == sps/(4*α)  # another case when full equation is undefined
            h[i] = (α/sqrt(2*sps)) * ( (1+(2/π))*sin(π/(4*α)) + (1-(2/π))*cos(π/(4*α)) ) 
        else
            h[i] = (4*α)/(π*sqrt(sps)) * ( cos((1+α)*π*t/sps) + (sps/(4*α*t))*sin((1-α)*π*t/sps) ) / (1-(4*α*t/sps)^2)
        end
    end

    return h
end
