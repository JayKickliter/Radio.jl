using Base.Test

function test_decimate{Th, Tx}( ::Type{Th}, ::Type{Tx}, hLen, xLen, factor )
    h = rand( Th, hLen )
    x = rand( Tx, xLen )
    

    nativeResult = decimate( h, x, factor )
    baseResult   = Base.filt( h, one(Th), x )
    baseResult   = baseResult[1:factor:end]
    areApprox( nativeResult, baseResult )
end

function test_interpolate{Th, Tx}( ::Type{Th}, ::Type{Tx}, hLen, xLen, factor )
    h = rand( Th, hLen )
    h = isreal(h[1]) ? h : real(h).+im
    x = rand( Tx, xLen )
    xx = upsample( x, factor )
    
    nativeResult = interpolate( h, x, factor )
    baseResult   = Base.filt( h, one(Th), xx )
    areApprox( nativeResult, baseResult )
end

filter_types = [    Float32,
                    Float64,
                    Complex64,
                    Complex128 ]

tap_lengths = rand(1:128, 5)
factors     = rand(1:32,  5)

for theType in filter_types, tap_len in tap_lengths, factor in factors
    sig_len = factor*tap_len
    println( "Testing interpolate, type: $theType, nTaps: $tap_len, nX: $sig_len, factor: $factor  " )
    @test test_interpolate( theType, theType, tap_len, sig_len, factor )
end

for theType in filter_types, factor in factors
    @test test_decimate( theType, theType, 4*factor, 10*factor*tap_len, factor )
end