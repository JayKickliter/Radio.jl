function rms{T<:Real}(s::Array{T,1})
   m = zero(T)
   for si in s
       m += si * si
   end
   sqrt(m/length(s))
end

function rms{T<:Real}(s::Array{Complex{T},1})
   m = zero(T)
   for si in s
       m += abs2(si)
   end
   sqrt(m/length(s))
end