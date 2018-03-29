# Acoustics

[![Build Status](https://travis-ci.org/thefirecode/Acoustics.jl.svg?branch=master)](https://travis-ci.org/thefirecode/Acoustics.jl)

[![Coverage Status](https://coveralls.io/repos/thefirecode/Acoustics.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/thefirecode/Acoustics.jl?branch=master)

[![codecov.io](http://codecov.io/github/thefirecode/Acoustics.jl/coverage.svg?branch=master)](http://codecov.io/github/thefirecode/Acoustics.jl?branch=master)
#Acoustics.jl

##Examples

**Usage example**
```
using Acoustics,LibSndFile
pwd() #this will show where julia is looking for files
cd("path you would like to go to")
a=load("test.wav") #it can be any file supported by
RT(a,60) #will take some time but, will return a broadband reverberation time, not the midband

example
using Acoustics,LibSndFile
cd("C:\\Users\\phenix")
a=load("Sports Hall UYork 441.aiff")[:,1]
RT(a,60)


```

**Generic Options in all function**
```
every function independent of function specific parameters

"Function name" (source,....(arbitary argument),weighting,band)
weighting - the weighting of the input audio files. the following weightings are defined z,a,b,c,d

example: RT(a,60,"z")

bands-this allows for multiband processing to be supported currently on 1/3 octave bands are supported

example: example: RT(a,60,"z","1/3")

```

##Functions
