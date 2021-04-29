using Acoustics,Test

single=Acoustics.Acoustic([1.0,2.0,3.0],convert(Float32,48000),"Bridgetower",convert(UInt16,1),0x04,3)

# write your own tests here
@test typeof(single.samples)==Array{Float64,1}
@test typeof(single.samplerate)==Float32
@test typeof(single.name)==String
@test typeof(single.channels)==UInt16
@test typeof(single.format)==UInt8
@test typeof(single.l_samples)==Int64
@test single.format==0x04


