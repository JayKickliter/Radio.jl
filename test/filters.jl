reload("/Users/jkickliter/.julia/v0.4/Radio/src/Filter/FIRFilter.jl")

using Base.Test
import Multirate





#==============================================================================#
#               ____ _ _  _ ____ _    ____    ____ ____ ___ ____               #
#               [__  | |\ | | __ |    |___    |__/ |__|  |  |___               #
#               ___] | | \| |__] |___ |___    |  \ |  |  |  |___               #
#==============================================================================#

function test_singlerate( h, x )
    @printf( "\nTesting single rate")
    xLen = length( x )
    pivotPoint = min( 100, ifloor( xLen/4 ))
    x1 = x[ 1 : pivotPoint ]
    x2 = x[ pivotPoint+1 : end ]    

    @printf( "\n\tBase's filt\n\t\t")
    @time naiveResult = Base.filt( h, 1.0, x )
    
    @printf( "\n\tMultirate's stateful filt\n\t\t")
    self = Multirate.FIRFilter( h, 1//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]
    
    @printf( "\n\tMultirate's stateful filt. Piecewise for first %d inpouts\n\t\t", length( x1 ) )
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
    @printf( "\nTesting decimation")
    xLen = length( x )
    pivotPoint = min( 100, ifloor( xLen/4 ))
    x1 = x[ 1 : pivotPoint ]
    x2 = x[ pivotPoint+1 : end ]    
    
    @printf( "\n\tNaive decimation\n\t\t")
    @time begin
        naiveResult   = Base.filt( h, one(eltype(h)), x )
        naiveResult   = naiveResult[1:decimation:end]
    end

    @printf( "\n\tMultirate's stateful decimation\n\t\t")
    self = Multirate.FIRFilter( h, 1//decimation )    
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end    
    statefulResult = [ y1, y2 ]
    
    @printf( "\n\tMultirate's stateful decimation. Piecewise for first %d inpouts\n\t\t", length( x1 ) )
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
    xLen = length( x )
    pivotPoint = min( 100, ifloor( xLen/4 ))
    x1 = x[ 1 : pivotPoint ]
    x2 = x[ pivotPoint+1 : end ]    
    
    @printf( "\nTesting interpolation")

    @printf( "\n\tNaive interpolation\n\t\t")
    @time begin
        xZeroStuffed = zeros( eltype(x), xLen * interpolation )
        for n = 0:xLen-1;
            xZeroStuffed[ n*interpolation+1 ] = x[ n+1 ]
        end
        naiveResult = Base.filt( h, one(eltype(h)), xZeroStuffed )        
    end
    
    @printf( "\n\tMultirate's stateful interpolation\n\t\t")
    self = Multirate.FIRFilter( h, interpolation//1 )
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )        
    end
    statefulResult = [ y1, y2 ]
    
    @printf( "\n\tMultirate's stateful interpolation. Piecewise for first %d inpouts\n\t\t", length( x1 ) )
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
    pivotPoint = min( 100, ifloor( xLen/4 ))
    x1         = x[ 1 : pivotPoint ]
    x2         = x[ pivotPoint+1 : end ]
    upfactor   = num( ratio )
    downfactor = den( ratio )
    
    @printf( "Multirate's stateful rational resampling\n\t")
    self = Multirate.FIRFilter( h, ratio );
    @time begin
        y1 = Multirate.filt( self, x1 )
        y2 = Multirate.filt( self, x2 )
    end
    statefulResult = [ y1, y2 ]
    
    @printf( "Naive resampling\n\t")
    @time begin
        xx         = zeros( length(x) * upfactor );
        naiveResult = similar( x, int( ceil( length(x) * ratio )))
        
        for n = 0:length(x)-1;
            xx[ n*upfactor+1 ] = x[ n+1 ];
        end
                
        naiveResultInterpolated = Base.filt( h, 1.0, xx );
        naiveResult = [ naiveResultInterpolated[n] for n = 1:downfactor:length( naiveResultInterpolated ) ]
    end
    
    display( [ naiveResult statefulResult ] )

    areApprox( naiveResult, statefulResult )
end

h = rand( 25 );
x = rand( int(100e3) );
# h = Float64[10:-1:1];
# x = Float64[1:100];
Th = eltype( h )
Tx = eltype( x )
interpolation = 3
decimation = 4

# Multirate.filt( Multirate.FIRFilter( h ), x[1:min(100, length(x))]);
# @test test_singlerate( h, x );

# Multirate.filt( Multirate.FIRFilter( h, 1//decimation ), x[1:min(100, length(x))]);
# @test test_decimate( h, x, decimation );

# Multirate.filt( Multirate.FIRFilter( h, interpolation//1 ), x[1:min(100, length(x))]);
# @test test_interpolate( h, x, interpolation )

ratio = interpolation//decimation
Multirate.filt( Multirate.FIRFilter( h, ratio ), x[1:min(100, length(x))]);
@test test_rational( h, x, ratio )
