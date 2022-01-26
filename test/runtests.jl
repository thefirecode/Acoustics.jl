using Acoustics,Test,DSP

single=Acoustics.Acoustic([1.0,2.0,3.0],convert(Float32,48000),"Bridgetower",convert(UInt16,1),0x04,3)

# write your own tests here
@test typeof(single.samples)==Array{Float64,1}
@test typeof(single.samplerate)==Float32
@test typeof(single.name)==String
@test typeof(single.channels)==UInt16
@test typeof(single.format)==UInt8
@test typeof(single.l_samples)==Int64
@test single.format==0x04

#Octave Band Filter Verification
@test filter_verify(1,44100,1)
@test filter_verify(1,48000,1)
@test filter_verify(1,88200,1)
@test filter_verify(1,96000,1)
#@test filter_verify(1,176400,1)
#@test filter_verify(1,192000,1)
@test filter_verify(2,44100,1)
@test filter_verify(2,48000,1)
@test filter_verify(2,88200,1)
@test filter_verify(2,96000,1)
#@test filter_verify(2,176400,1)
#@test filter_verify(2,192000,1)
@test filter_verify(3,44100,1)
@test filter_verify(3,48000,1)
@test filter_verify(3,88200,1)
@test filter_verify(3,96000,1)
#@test filter_verify(3,176400,1)
#@test filter_verify(3,192000,1)
@test filter_verify(4,44100,1)
@test filter_verify(4,48000,1)
@test filter_verify(4,88200,1)
@test filter_verify(4,96000,1)
#@test filter_verify(4,176400,1)
#@test filter_verify(4,192000,1)
@test filter_verify(5,44100,1)
@test filter_verify(5,48000,1)
@test filter_verify(5,88200,1)
@test filter_verify(5,96000,1)
#@test filter_verify(5,176400,1)
#@test filter_verify(5,192000,1)
@test filter_verify(6,44100,1)
@test filter_verify(6,48000,1)
@test filter_verify(6,88200,1)
@test filter_verify(6,96000,1)
#@test filter_verify(6,176400,1)
#@test filter_verify(6,192000,1)
