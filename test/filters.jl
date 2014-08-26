reload("/Users/jkickliter/.julia/v0.3/Radio/src/Filter/FIRFilter.jl")

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

    @printf( "\n\tBase.filt\n\t\t")
    @time naiveResult = Base.filt( h, 1.0, x )

    @printf( "\n\tMultirate.filt filt. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
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

    @printf( "\n\tNaive decimation with Base.filt\n\t\t")
    @time begin
        naiveResult   = Base.filt( h, one(eltype(h)), x )
        naiveResult   = naiveResult[1:decimation:end]
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

    @printf( "\n\tNaive interpolation with Base.filt\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        naiveResult = Base.filt( h, one(eltype(h)), xZeroStuffed )
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

    println(); println()
    println( "Testing rational resampling. xLen = $xLen, hLen = $hLen, ratio = $ratio")

    @printf( "\n\tNaive rational resampling with Base.filt\n\t\t")
    @time begin
        xx          = zeros( resultType, length(x) * upfactor )
        naiveResult = Array( resultType, int( ceil( length(x) * ratio )))

        for n = 0:length(x)-1;
            xx[ n*upfactor+1 ] = x[ n+1 ]
        end

        naiveResultInterpolated = Base.filt( h, one(eltype(h)), xx );
        naiveResult = [ naiveResultInterpolated[n] for n = 1:downfactor:length( naiveResultInterpolated ) ]
    end

    @printf( "\n\tMultirate.filt rational resampling. length( x1 ) = %d, length( x2 ) = %d\n\t\t", length( x1 ), length( x2 ) )
    self = Multirate.FIRFilter( h, ratio )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]

    @printf( "\n\tMultirate.filt rational. Piecewise for first %d inputs, then x2\n\t\t", length( x1 ) )
    self = Multirate.FIRFilter( h, ratio )
    y1 = similar( x, 0 )
    @time begin
        for i in 1:length(x1)
            append!( y1, Multirate.filt( self, x1[i:i] ) )
        end
        y2 = Multirate.filt( self, x2 )
    end
    piecewiseResult = [ y1, y2 ]

    areApprox( naiveResult, statefulResult ) && areApprox( naiveResult, piecewiseResult )
end

h  = rand( rand(27:128,1)[1] );
x  = rand( int(1e6)+rand( 1:100, 1 )[1] );

interpolation = rand(1:32, 1)[1]
decimation    = rand(1:32, 1)[1]
ratio         = interpolation//decimation

Multirate.filt( Multirate.FIRFilter( h ), x[1:min(100, length(x))]);
@test test_singlerate( h, x );

Multirate.filt( Multirate.FIRFilter( h, 1//decimation ), x[1:min(100, length(x))]);
@test test_decimate( h, x, decimation );

Multirate.filt( Multirate.FIRFilter( h, interpolation//1 ), x[1:min(100, length(x))]);
@test test_interpolate( h, x, interpolation )

Multirate.filt( Multirate.FIRFilter( h, ratio ), x[1:min(100, length(x))]);
@test test_rational( h, x, ratio )
