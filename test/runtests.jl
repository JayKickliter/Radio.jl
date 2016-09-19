using Radio
using Base.Test

qam4 = QAM(4)
@test qam4.M == 4
@test qam4.m == 2
@test qam4.α ≈ √2
@test qam4.β ≈ -0.5
@test qam4.constellation ≈ [-1-im,-1+im,1-im,1+im] / √2f0
@test qam4.grayConstellation ≈ [-1-im,-1+im,1-im,1+im] / √2f0
@test qam4.bitsMap ≈ [0,1,2,3]

qam16 = QAM(16)
@test qam16.M == 16
@test qam16.m == 4
@test qam16.α ≈ 2/√10
@test qam16.β ≈ -1.5
@test qam16.constellation ≈ [complex(i,q)/√10f0 for q in -3:2:3, i in -3:2:3][:]
@test qam16.grayConstellation ≈ ([complex(i,q)/√10f0 for q in [-3,-1,3,1], i in [-3,-1,3,1]][:])
@test qam16.bitsMap == [i*4+q for q in [0,1,3,2], i in [0,1,3,2]][:]

qam64 = QAM(64)
@test qam64.M == 64
@test qam64.m == 6
@test qam64.α ≈ 2/√42
@test qam64.β ≈ -3.5
@test qam64.constellation ≈ [complex(i,q)/√42f0 for q in -7:2:7, i in -7:2:7][:]
@test qam64.grayConstellation ≈ ([complex(i,q)/√42f0 for q in [-7,-5,-1,-3,7,5,1,3], i in [-7,-5,-1,-3,7,5,1,3]][:])
@test qam64.bitsMap == [i*8+q for q in [0,1,3,2,6,7,5,4], i in [0,1,3,2,6,7,5,4]][:]
