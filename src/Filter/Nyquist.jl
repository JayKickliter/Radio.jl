#==============================================================================#
#                                Raised Cosine                                 #
#==============================================================================#
# α    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# k    = Sample period; Samples per symbol

function rcos( α, span, k)
    hlen = k*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h = Array( Float64, hlen )
    for i = 1:hlen
        t = i - 1 - (hlen - 1)/2
        if abs(α*t/k) == 0.5 # check for case where (1-(2*α*t/k)^2 is zero, which would cause divide by zero
            h[i] = 0.5
        else
            h[i] = sinc(t/k)*cos(pi*α*t/k)/(1-(2*α*t/k)^2)
        end
    end 

    return h
end

#==============================================================================#
#                             Root Raised Cosine                               #
#==============================================================================#
# α    = Rolloff factor
# span = filter span. How many symbols the filter should cover
# k    = Sample period; Samples per symbol

function rrcos( α, span, k)
    hlen = k*span
    hlen = iseven(hlen) ? hlen+1 : hlen
    h    = Array( Float64, hlen )
    for i = 1:hlen
        t = i - (hlen - 1)/2 - 1 
        if t == 0
            h[i] = (1/sqrt(k)) * (1+α*(4/pi -1))
        elseif abs(t) == k/(4*α)
            h[i] = (α/sqrt(2*k)) * ( (1+(2/pi))*sin(pi/(4*α)) + (1-(2/pi))*cos(pi/(4*α)) ) 
        else
            h[i] = (4*α)/(pi*sqrt(k)) * ( cos((1+α)*pi*t/k) + (k/(4*α*t))*sin((1-α)*pi*t/k) ) / (1-(4*α*t/k)^2)
        end
    end

    return h
end