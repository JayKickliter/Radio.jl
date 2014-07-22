import IPPDSP

function decimate!{T}( buffer::Vector{T}, h::Vector{T}, x::Vector{T}, M::Integer )
    xLen   = length( x )    
    hLen   = length( h )
    bufLen = length( buffer )
    yLen   = int(xLen / M)    
    
    xLen   % M == 0     || error( "signal length % decimation must be 0" )
    bufLen * M >= xLen  || error( "buffer lenght must be >= signal length * decimation" )        

    criticalYidx = int(ceil(hLen / M))
        
    xIdx = 1
    
    for yIdx = 1:criticalYidx
            
        acc  = zero(T)
        
        kMax = xIdx < hLen ? xIdx : hLen
                
        for k = 1:kMax
            @inbounds acc += h[ k ] * x[ xIdx+1-k ]
        end
        
        @inbounds buffer[yIdx] = acc
        xIdx += M
    end
    
    h = flipud(h)
    
    xIdx -= hLen
    
    for yIdx = criticalYidx+1:yLen
            
        acc  = zero(T)
                
        for k = 1:hLen
            @inbounds acc += h[ k ] * x[ xIdx+k ]
        end
        
        @inbounds buffer[yIdx] = acc
        xIdx += M
    end
        
    return buffer
end

decimate{T}( h::Vector{T}, x::Vector{T}, M::Integer ) = decimate!( similar(x, int(length(x)/M)), h::Vector{T}, x::Vector{T}, M::Integer )


function testdecimate( Th, Tx, hLen, xLen, factor )
    h = rand( Th, hLen )
    x = rand( Tx, xLen )
    
    bufLen = int( xLen/factor )
    
    nativeBuffer =  similar( x, bufLen )
    # ippBuffer =  similar( x, bufLen )
    
    @time decimate!( nativeBuffer, h, x, factor )    
    # @time IPPDSP.filt!( ippBuffer, h, x, 1, factor )
    
    # areApprox( nativeBuffer, ippBuffer )
end

factor = 7
xLen   = int(1e6*factor)
hLen   = 32
@profile testdecimate( Float32, Float32, 32, xLen, factor )
