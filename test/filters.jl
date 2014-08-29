reload(expanduser("~/.julia/v0.3/Radio/src/Filter/FIRFilter.jl"))

using Base.Test
import Multirate



#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function test_singlerate( h, x )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    println(); println()
    println( "Testing single-rate fitering. xLen = $xLen, hLen = $hLen")

    @printf( "\nBase.filt\n\t\t")
    @time naiveResult = Base.filt( h, 1.0, x )

    @printf( "\nMultirate.filt filt. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\nMultirate.filt filt. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    @time begin
        for i in 1:length(x1)
            y1[i] = Multirate.filt( self, x1[i:i] )[1]
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]


    areApprox( statefulResult, naiveResult ) && areApprox( piecewiseResult, naiveResult )
end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function test_decimate( h, x, decimation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    println(); println()
    println( "Testing decimation. xLen = $xLen, hLen = $hLen, decimation = $decimation")

    @printf( "\nNaive decimation with Base.filt\n\t\t")
    @time begin
        naiveResult   = Base.filt( h, one(eltype(h)), x )
        naiveResult   = naiveResult[1:decimation:end]
    end

    @printf( "\nMultirate.filt decimation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//decimation )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\nMultirate.filt decimation. Piecewise for first %d inputs.\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    areApprox( naiveResult, statefulResult ) && areApprox( piecewiseResult, naiveResult )
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function test_interpolate( h, x, interpolation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    println(); println()
    println( "Testing interpolation. xLen = $xLen, hLen = $hLen, interpolation = $interpolation")

    @printf( "\nNaive interpolation with Base.filt\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        naiveResult = Base.filt( h, one(eltype(h)), xZeroStuffed )
    end

    @printf( "\nMultirate.filt interpolation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, interpolation//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\nMultirate.filt interpolation. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    areApprox( naiveResult, statefulResult ) && areApprox( piecewiseResult, naiveResult )
end




#==============================================================================#
#           ____ ____ ___     ____ ____ ____ ____ _  _ ___  _    ____          #
#           |__/ |__|  |      |__/ |___ [__  |__| |\/| |__] |    |___          #
#           |  \ |  |  |  .   |  \ |___ ___] |  | |  | |    |___ |___          #
#==============================================================================#

function test_rational( h, x, ratio )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]
    upfactor   = num( ratio )
    downfactor = den( ratio )
    resultType = promote_type( eltype(h), eltype(x) )

   # println(); println()
   # println( "Testing rational resampling. xLen = $xLen, hLen = $hLen, ratio = $ratio")

    # @printf( "\nIPPDSP.filt\n\t\t")
    # # newXlen = xLen - mod( xLen, downfactor )
    # self = IPPDSP.FIRFilter( eltype(x), h, upfactor, downfactor )
    # @time ippResult = IPPDSP.filt( self, x )
    # display(ippResult)

    @printf( "\nNaive rational resampling with Base.filt\n\t\t")
    @time begin
        xx          = zeros( resultType, length(x) * upfactor )
        naiveResult = Array( resultType, int( ceil( length(x) * ratio )))

        for n = 0:length(x)-1;
            xx[ n*upfactor+1 ] = x[ n+1 ]
        end

        naiveResultInterpolated = Base.filt( h, one(eltype(h)), xx );
        naiveResult = [ naiveResultInterpolated[n] for n = 1:downfactor:length( naiveResultInterpolated ) ]
    end    
    # display( naiveResult )
    
   # println( "_ _  _ ___  _  _ ___    ___  ____ ____ ____ ____ ____ ____ ____ _ ____ _  _ " )
   # println( "| |\\ | |__] |  |  |     |__] |__/ |  | | __ |__/ |___ [__  [__  | |  | |\\ | " )
   # println( "| | \\| |    |__|  |     |    |  \\ |__| |__] |  \\ |___ ___] ___] | |__| | \\| " )
    
    z = [1:xLen]'
    z = repmat( z, upfactor, 1 )
    z = reshape( z, xLen*upfactor, 1)
    z = [ z[n] for n in 1:downfactor:length(z) ]
   # display( [1:length(z) z ] )
    
   # println( "____ _ _  _ ____ _    ____    ___  ____ ____ ____ " )
   # println( "[__  | |\\ | | __ |    |___    |__] |__| [__  [__  " )
   # println( "___] | | \\| |__] |___ |___    |    |  | ___] ___] " )
    
    self = Multirate.FIRFilter( h, ratio )
    @time singlepassResult = Multirate.filt( self, x )
    
    
   # println( "___ _ _ _ ____    ___  ____ ____ ___ " )
   # println( " |  | | | |  |    |__] |__| |__/  |  " )
   # println( " |  |_|_| |__|    |    |  | |  \\  |  " )
                                         
    
    @printf( "\nMultirate.filt rational resampling. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, ratio )
    @time begin
        s1 = Multirate.filt( self, x1 )
        s2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ s1, s2 ]
    # display(statefulResult)
    
   # println()
   # println()
   # println( "___  _ ____ ____ ____ _ _ _ ____ _ ____ ____ " )
   # println( "|__] | |___ |    |___ | | | |___ | [__  |___ " )
   # println( "|    | |___ |___ |___ |_|_| |___ | ___] |___ " )
                                                 
    @printf( "\nMultirate.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    self = Multirate.FIRFilter( h, ratio )
    # println( "pfb = $(self.kernel.pfb)")
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x)
            append!( y1, Multirate.filt( self, x[i:i] ) )
        end
    end
    piecewiseResult = y1
    # display( piecewiseResult )


    
    if areApprox( naiveResult, singlepassResult ) && areApprox( naiveResult, statefulResult ) && areApprox( naiveResult, piecewiseResult )
        return true
    end

    
    st1 = [ s1, zeros(eltype(s1), length(s2)) ]
    st2 = [ zeros(eltype(s2), length(s1)) , s2 ]
    
    commonLen = min( length(naiveResult), length( singlepassResult ), length(statefulResult), length(piecewiseResult) )
   # display( [ [1:commonLen] naiveResult[1:commonLen] singlepassResult[1:commonLen] statefulResult[1:commonLen] #=st1[1:commonLen] st2[1:commonLen]=# piecewiseResult[1:commonLen] ])

    return false
end





# Multirate.filt( Multirate.FIRFilter( h ), x[1:min(100, length(x))]);
# @test test_singlerate( h, x );
#
# Multirate.filt( Multirate.FIRFilter( h, 1//decimation ), x[1:min(100, length(x))]);
# @test test_decimate( h, x, decimation );
#
# Multirate.filt( Multirate.FIRFilter( h, interpolation//1 ), x[1:min(100, length(x))]);
# @test test_interpolate( h, x, interpolation )


function run_tests()
    for i in 1:100
        for j in 1:10; println(); end
        println( "Test number $i")
        interpolation = rand(2:32, 1)[1]
        decimation    = rand(2:32, 1)[1]
        # interpolation   = 7
        # decimation      = 11
        h             = rand(Float32, rand(1:128,1)[1] )
        x             = rand(Float32, int(1e3)+rand( 1:100, 1 )[1] )
        # x               = [ 1.0:40 ]
        # h = [1.0:40]
        ratio           = interpolation//decimation
        while ( num( ratio ) == 1 || den( ratio ) == 1 )
            ratio = rand(2:32, 1)[1]//rand(2:32, 1)[1]
        end
        
        
        # h           = [ 1.0, zeros(3) ]        
        # h           = [1.0:rand(num(ratio):64, 1)[1]]
        # interpolation = num(ratio)
        # decimation    = den(ratio)
        # hLen          = 2*interpolation
        # tapsPerφ      = 4
        
        x    = [ 1.0:rand(103:192, 1)[1] ]
        # x      = [ 1.0 : 20*decimation/interpolation ]

        # hLen = interpolation*2
        # h      = zeros( 2, interpolation )
        # h[1,:] = 1
        # h      = vec(reshape(h', interpolation*2, 1))
        
        # Multirate.filt( Multirate.FIRFilter( h, ratio ), x[1:min(100, length(x))])
        # xLen = length(x)
        # x = x[1:xLen - mod( xLen, den( ratio ) )]
        @test test_rational( h, x, ratio )
    end    
end

function debug_test()
    for i in 1:100
        # for j in 1:10; println(); end
        # println( "Test number $i")
        # interpolation = rand(2:32, 1)[1]
        # decimation    = rand(2:32, 1)[1]
        interpolation   = 3
        decimation      = 4
        # h             = rand(Float32, rand(1:128,1)[1] )
        # x             = rand(Float32, int(1e3)+rand( 1:100, 1 )[1] )
        # x               = [ 1.0:40 ]
        # h = [1.0:40]
        ratio           = interpolation//decimation
        # while ( num( ratio ) == 1 || den( ratio ) == 1 )
            # ratio = rand(2:32, 1)[1]//rand(2:32, 1)[1]
        # end
        
        
        h           = [ 1.0, 4*interpolation-1 ]        
        # h           = [1.0:rand(num(ratio):64, 1)[1]]
        # interpolation = num(ratio)
        # decimation    = den(ratio)
        # hLen          = 2*interpolation
        # tapsPerφ      = 4
        
        # x    = [ 1.0:rand(103:192, 1)[1] ]
        x      = [ 1.0 : 20*decimation/interpolation ]

        # hLen = interpolation*2
        # h      = zeros( 2, interpolation )
        # h[1,:] = 1
        # h      = vec(reshape(h', interpolation*2, 1))
        
        # Multirate.filt( Multirate.FIRFilter( h, ratio ), x[1:min(100, length(x))])
        # xLen = length(x)
        # x = x[1:xLen - mod( xLen, den( ratio ) )]
        @test test_rational( h, x, ratio )
    end    
end

function test_phasenext()
    
    for interpolation in 1:8
        for decimation in 1:8
            # println()
            # println()
            # println()
            # println( "ratio = $(interpolation//decimation)")
            ratio           = interpolation//decimation
            interpolation   = num(ratio)
            decimation      = den(ratio)
            x               = repmat( [1:interpolation], decimation )
            reference = [ x[n] for n = 1:decimation:length( x ) ]
            result = [ 1 ]
            for i in 2:interpolation
                append!( result, [ Multirate.nextphase( result[end], ratio ) ] )
            end
            @test areApprox( reference, result )
            # display( [ reference result ] )
        end
    end
end


# debug_test()
run_tests()
# test_phasenext()
