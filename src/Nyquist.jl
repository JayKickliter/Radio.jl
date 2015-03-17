#==============================================================================#
#                                Raised Cosine                                 #
#==============================================================================#
# β    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# sps  = Samples per symbol

function rcos( β, span, sps )
    # TODO: add argument valdiation
    hlen = sps*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h    = Array( Float64, hlen )
    
    for i = 1:hlen
        t = i - 1 - (hlen - 1)/2
        if abs(β*t/sps) == 0.5      # check for case where (1-(2*β*t/sps)^2 is zero, which would cause divide by zero
            h[i] = sinc(t/sps)
        else
            h[i] = sinc(t/sps)*cos(π*β*t/sps)/(1-(2*β*t/sps)^2)
        end
    end 

    return h
end

#==============================================================================#
#                             Root Raised Cosine                               #
#==============================================================================#
# β    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# sps  = Samples per symbol

function rrcos( β, span, sps )
    # TODO: add argument valdiation
    hlen = sps*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h    = Array( Float64, hlen )
    
    for i = 1:hlen
        t = i - (hlen - 1)/2 - 1 
        if t == 0                   # full equation is undefined at t = 0
            h[i] = (1/sqrt(sps)) * (1+β*(4/π -1))  
        elseif abs(t) == sps/(4*β)  # another case when full equation is undefined
            h[i] = (β/sqrt(2*sps)) * ( (1+(2/π))*sin(π/(4*β)) + (1-(2/π))*cos(π/(4*β)) ) 
        else
            h[i] = (4*β)/(π*sqrt(sps)) * ( cos((1+β)*π*t/sps) + (sps/(4*β*t))*sin((1-β)*π*t/sps) ) / (1-(4*β*t/sps)^2)
        end
    end
  
    return h
end