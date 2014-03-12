type IQ{T}
    centerFrequency::Real
    sampleRate::Real
    sampleCount::Int
    data::Vector{Complex{T}}
end  