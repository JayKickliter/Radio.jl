function rms{T<:Real}(s::Array{T,1})
   m = zero(T)
   for si in s
       m += si * si
   end
   sqrt(m/length(s))
end