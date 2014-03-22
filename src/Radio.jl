# /Volumes/HD1/Repositories/Radio.jl/src/Radio.jl

module Radio

include( "Filter/FIRDesign.jl"  )
include( "Filter/FIRFilter.jl"  )
include( "Filter/Nyquist.jl"    )  
include( "Filter/Window.jl"     )
include( "Math/CZT.jl"          )
include( "Math/Power.jl"        )
include( "Modulation/PSK.jl"    )
include( "Random/Noise.jl"      )
include( "Support/Types.jl"     )

using Radio.Math, Radio.Filter, Radio.Modulation, Radio.Random

export
    # Filter
    firdes, rcos, rrcos, blackaman, hamming, hann, kaiser, rectangle,
    
    # Math
    czt, rms,
    
    # Random
    wgn

end # Radio module