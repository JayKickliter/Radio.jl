reload(expanduser("~/.julia/v0.3/Radio/src/Filter/FIRFilter.jl"))

using Base.Test
import Multirate
import IPPDSP
import DSP


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
    
    println( "____ _ _  _ ____ _    ____    ____ ____ ___ ____" )
    println( "[__  | |\\ | | __ |    |___    |__/ |__|  |  |___" )
    println( "___] | | \\| |__] |___ |___    |  \\ |  |  |  |___" )
    println()
    println( "Testing single-rate fitering. xLen = $xLen, hLen = $hLen")

    @printf( "\n\tBase.filt\n\t\t")
    @time baseResult = Base.filt( h, 1.0, x )

    @printf( "\n\tDSP.firfilt\n\t\t")
    @time dspResult = DSP.firfilt( h, x )
    
    @printf( "\n\tMultirate.filt. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt filt. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    @time begin
        for i in 1:length(x1)
            y1[i] = Multirate.filt( self, x1[i:i] )[1]
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]


    areApprox( statefulResult, baseResult ) && areApprox( piecewiseResult, baseResult )
end




#==============================================================================#
#                      ___  ____ ____ _ _  _ ____ ___ ____                     #
#                      |  \ |___ |    | |\/| |__|  |  |___                     #
#                      |__/ |___ |___ | |  | |  |  |  |___                     #
#==============================================================================#

function test_decimation( h, x, decimation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    println(); println()
    println( "___  ____ ____ _ _  _ ____ ___ _ ____ _  _ " )
    println( "|  \\ |___ |    | |\\/| |__|  |  | |  | |\\ | " )
    println( "|__/ |___ |___ | |  | |  |  |  | |__| | \\| " )
    println()
    println( "Testing decimation. xLen = $xLen, hLen = $hLen, decimation = $decimation")

    @printf( "\n\tNaive decimation with Base.filt\n\t\t")
    @time begin
        baseResult   = Base.filt( h, one(eltype(h)), x )
        baseResult   = baseResult[1:decimation:end]
    end
    
    @printf( "\n\tNaive decimation with DSP.firfilt\n\t\t")
    hDSP = h
    eltype( hDSP ) == eltype( x ) || (hDSP = eltype( x )[ he for he in h ])
    @time begin
        dspResult = DSP.firfilt( hDSP, x )
        dspResult = dspResult[1:decimation:end]
    end

    @printf( "\n\tMultirate.filt decimation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, 1//decimation )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt decimation. Piecewise for first %d inputs.\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    if areApprox( baseResult, statefulResult ) && areApprox( baseResult, piecewiseResult )
        return true
    end
    
    display( [ baseResult statefulResult piecewiseResult ] )
    return false
    
end




#==============================================================================#
#               _ _  _ ___ ____ ____ ___  _    ____ ____ ___ ____              #
#               | |\ |  |  |___ |__/ |__] |    |  | |__|  |  |___              #
#               | | \|  |  |___ |  \ |    |___ |__| |  |  |  |___              #
#==============================================================================#

function test_interpolation( h, x, interpolation )
    xLen       = length( x )
    hLen       = length( h )
    pivotPoint = min( rand(50:150, 1)[1], ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]

    println(); println()
    println( "_ _  _ ___ ____ ____ ___  _    ____ ____ ___ _ ____ _  _ " )
    println( "| |\\ |  |  |___ |__/ |__] |    |  | |__|  |  | |  | |\\ | " )
    println( "| | \\|  |  |___ |  \\ |    |___ |__| |  |  |  | |__| | \\| " )
    println()                                                             
    println( "Testing interpolation. xLen = $xLen, hLen = $hLen, interpolation = $interpolation")

    @printf( "\n\tNaive interpolation with Base.filt\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        baseResult = Base.filt( h, one(eltype(h)), xZeroStuffed )
    end
    
    @printf( "\n\tNaive interpolation with DSP.firfilt\n\t\t")
    hDSP = h
    eltype( hDSP ) == eltype( x ) || (hDSP = eltype( x )[ he for he in h ])
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        dspResult = DSP.firfilt( hDSP, xZeroStuffed )
    end

    @printf( "\n\tMultirate.filt interpolation. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, interpolation//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt interpolation. Piecewise for first %d inputs\n\t\t", length( x1 ) )
    Multirate.reset( self )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    areApprox( baseResult, statefulResult ) && areApprox( piecewiseResult, baseResult )
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

    println(); println()
    println( "     ____ ____ ___ _ ____ _  _ ____ _    " )
    println( "     |__/ |__|  |  | |  | |\\ | |__| |    " )
    println( "     |  \\ |  |  |  | |__| | \\| |  | |___ " )
    println( "                                         " )
    println( "____ ____ ____ ____ _  _ ___  _    ____ ____ " )
    println( "|__/ |___ [__  |__| |\\/| |__] |    |___ |__/ " )
    println( "|  \\ |___ ___] |  | |  | |    |___ |___ |  \\ " )
    println()                                         
    println( "Testing rational resampling. xLen = $xLen, hLen = $hLen, ratio = $ratio")

    @printf( "\n\tNaive rational resampling with Base.filt\n\t\t")
    @time begin
        xStuffed   = zeros( resultType, length(x) * upfactor )
        baseResult = Array( resultType, int( ceil( length(x) * ratio )))

        for n = 0:length(x)-1;
            xStuffed[ n*upfactor+1 ] = x[ n+1 ]
        end 

        baseResult = Base.filt( h, one(eltype(h)), xStuffed );
        baseResult = [ baseResult[n] for n = 1:downfactor:length( baseResult ) ]
    end    

    @printf( "\n\tNaive rational resampling DSP.firfilt\n\t\t")
    hDSP = h
    eltype( hDSP ) == eltype( x ) || (hDSP = eltype( x )[ he for he in h ])
    @time begin
        xStuffed  = zeros( resultType, length(x) * upfactor )
        dspResult = Array( resultType, int( ceil( length(x) * ratio )))

        for n = 0:length(x)-1;
            xStuffed[ n*upfactor+1 ] = x[ n+1 ]
        end 

        dspResult = DSP.firfilt( hDSP, xStuffed );
        dspResult = [ dspResult[n] for n = 1:downfactor:length( dspResult ) ]
    end    

    @printf( "\n\tIPPDSP.filt\n\t\t")
    self = IPPDSP.FIRFilter( eltype(x), h, upfactor, downfactor )
    @time ippResult = IPPDSP.filt( self, x )

    @printf( "\n\tMultirate.filt rational resampling. Single pass over all %d elements.\n\t\t", xLen )
    self = Multirate.FIRFilter( h, ratio )
    @time singlepassResult = Multirate.filt( self, x )
    
    @printf( "\n\tMultirate.filt rational resampling. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, ratio )
    @time begin
        s1 = Multirate.filt( self, x1 )
        s2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ s1, s2 ]
                                                 
    @printf( "\n\tMultirate.filt rational. Piecewise for all %d inputs\n\t\t", length( x ) )
    self = Multirate.FIRFilter( h, ratio )
    # println( "pfb = $(self.kernel.pfb)")
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x)
            append!( y1, Multirate.filt( self, x[i:i] ) )
        end
    end
    piecewiseResult = y1
    
    if areApprox( baseResult, singlepassResult ) && areApprox( baseResult, statefulResult ) && areApprox( baseResult, piecewiseResult )
        return true
    end

    display( [  baseResult singlepassResult statefulResult ])

    return false
end




function test_all()    
    for interpolation in 1:32, decimation in 1:32
        h    = rand(Float32, rand(16:128,1)[1] )
        xLen = int(rand( 1000:2000, 1 )[1])
        xLen = xLen-mod( xLen, decimation )
        x    = rand( Float32, xLen )
    
        @test test_singlerate( h, x )
        @test test_decimation( h, x, decimation )
        @test test_interpolation( h, x, interpolation )
        
        ratio = interpolation//decimation
        
        if num(ratio) != 1 && den(ratio) != 1
            @test test_rational( h, x, ratio )           
        end
    end
end




function test_speed()
    ratio     = 29//27
    xLen      = int(1e6)
    xLen      = xLen - mod( xLen, den(ratio) )
    x         = rand( Float32, xLen )
    h         = rand( Float32, 56 )
    @test test_rational( h, x, ratio )
end




function test_phasenext()    
    for interpolation in 1:8
        for decimation in 1:8
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
        end
    end
end


# debug_test()
# run_tests()
# test_phasenext()
# test_speed()
test_all()
