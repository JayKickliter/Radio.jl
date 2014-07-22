import IPPDSP

function test_decimate( Th, Tx, hLen, xLen, factor )
    h = rand( Th, hLen )
    x = rand( Tx, xLen )
    

    print("Native: ")
    @time nativeResult = decimate( h, x, factor )

    print("IPP: ")
    @time ippResult = IPPDSP.filt( h, x, 1, factor )
    
    areApprox( nativeResult, ippResult )
end


