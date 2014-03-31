module Radio

# firfilt is needed by several Radio functions
# import DSP: firfilt

include( "Filter/FIRDesign.jl"  )
include( "Filter/FIRFilter.jl"  )
include( "Filter/Nyquist.jl"    )
include( "Filter/Resample.jl"   )
include( "Filter/Window.jl"     )
include( "Math/CZT.jl"          )
include( "Math/Misc.jl"         )
include( "Math/Power.jl"        )
include( "Modulation/PSK.jl"    )
include( "Random/Noise.jl"      )
include( "Support/Types.jl"     )
include( "Support/Graphics.jl"  )

export
    # Filter
    firdes, rcos, rrcos, blackman,
    hamming, hann, kaiser, rectangle,
    kaiserord, interpolate, upsample,
    
    # Math
    czt, rms, similarzeros,
    
    # Random
    wgn,
    
    # Modulation
    pskmod,
    
    # Support
    plot_response

end # Radio module