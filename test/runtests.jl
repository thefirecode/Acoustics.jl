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

#Octave Band Filter Verification
@test filter_verify(1,44100,1)
@test filter_verify(1,48000,1)
@test filter_verify(1,88200,1)
@test filter_verify(1,96000,1)
@test filter_verify(1,176400,1)
@test filter_verify(1,192000,1)

@test filter_verify(2,44100,1)
@test filter_verify(2,48000,1)
@test filter_verify(2,88200,1)
@test filter_verify(2,96000,1)
@test filter_verify(2,176400,1)
@test filter_verify(2,192000,1)

@test filter_verify(3,44100,1)
@test filter_verify(3,48000,1)
@test filter_verify(3,88200,1)
@test filter_verify(3,96000,1)
@test filter_verify(3,176400,1)
@test filter_verify(3,192000,1)

@test filter_verify(4,44100,1)
@test filter_verify(4,48000,1)
@test filter_verify(4,88200,1)
@test filter_verify(4,96000,1)
@test filter_verify(4,176400,1)
@test filter_verify(4,192000,1)

@test filter_verify(5,44100,1)
@test filter_verify(5,48000,1)
@test filter_verify(5,88200,1)
@test filter_verify(5,96000,1)
@test filter_verify(5,176400,1)
@test filter_verify(5,192000,1)

@test filter_verify(6,44100,1)
@test filter_verify(6,48000,1)
@test filter_verify(6,88200,1)
@test filter_verify(6,96000,1)
@test filter_verify(6,176400,1)
@test filter_verify(6,192000,1)

@test filter_verify(7,44100,1)
@test filter_verify(7,48000,1)
@test filter_verify(7,88200,1)
@test filter_verify(7,96000,1)
@test filter_verify(7,176400,1)
@test filter_verify(7,192000,1)

@test filter_verify(8,44100,1)
@test filter_verify(8,48000,1)
@test filter_verify(8,88200,1)
@test filter_verify(8,96000,1)
@test filter_verify(8,176400,1)
@test filter_verify(8,192000,1)

@test filter_verify(9,44100,1)
@test filter_verify(9,48000,1)
@test filter_verify(9,88200,1)
@test filter_verify(9,96000,1)
@test filter_verify(9,176400,1)
@test filter_verify(9,192000,1)

@test filter_verify(10,44100,1)
@test filter_verify(10,48000,1)
@test filter_verify(10,88200,1)
@test filter_verify(10,96000,1)
@test filter_verify(10,176400,1)
@test filter_verify(10,192000,1)

@test filter_verify(11,44100,1)
@test filter_verify(11,48000,1)
@test filter_verify(11,88200,1)
@test filter_verify(11,96000,1)
@test filter_verify(11,176400,1)
@test filter_verify(11,192000,1)

@test filter_verify(12,44100,1)
@test filter_verify(12,48000,1)
@test filter_verify(12,88200,1)
@test filter_verify(12,96000,1)
@test filter_verify(12,176400,1)
@test filter_verify(12,192000,1)

@test filter_verify(13,44100,1)
@test filter_verify(13,48000,1)
@test filter_verify(13,88200,1)
@test filter_verify(13,96000,1)
@test filter_verify(13,176400,1)
@test filter_verify(13,192000,1)

@test filter_verify(14,44100,1)
@test filter_verify(14,48000,1)
@test filter_verify(14,88200,1)
@test filter_verify(14,96000,1)
@test filter_verify(14,176400,1)
@test filter_verify(14,192000,1)

@test filter_verify(15,44100,1)
@test filter_verify(15,48000,1)
@test filter_verify(15,88200,1)
@test filter_verify(15,96000,1)
@test filter_verify(15,176400,1)
@test filter_verify(15,192000,1)

@test filter_verify(16,44100,1)
@test filter_verify(16,48000,1)
@test filter_verify(16,88200,1)
@test filter_verify(16,96000,1)
@test filter_verify(16,176400,1)
@test filter_verify(16,192000,1)
