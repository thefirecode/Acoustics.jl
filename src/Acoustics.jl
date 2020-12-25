module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,Reexport

export Leq,C,RT,D,Ts,sweep,deconvolve,EDT,acoustic_load,ST_late,ST_early,IACC,G,sweep_target

#this contian how to generate third octaves
include("bands.jl");

using Reexport
@reexport using .Bands


#=
0x00-omni-"omni"
0x01-omni(channel 1) figure 8 (channel 2)-"omni8"
0x02-G recording-"g"
0x03-Binaural-"bin"
0x04-Multichannel-"m" - assume multichannel omni
=#

struct Acoustic
	samples::Array{Float64}
	samplerate::Float32
	name::String
	channels::UInt16
	format::UInt8
	l_samples::Int64
end

function format_parser(x::String,chan::UInt16=0x0001)
(x=="m")&&return 0x04
chan>0x0002&&return 0x04
(chan==0x0001)&&return 0x00
(chan==0x0002)&(x=="omni8")&&return 0x01
(chan==0x0002)&(x=="g")&&return 0x02
(chan==0x0002)&(x=="bin")&&return 0x03
#errors this is just a section for error handling
(chan>0x0001)&(x=="omni")&&return error("Undefined format. Try another format")
((chan==0x0002)&!(x=="bin")&!(x=="g")&!(x=="omni8"))&&return error("Undefined format. Try another format")
!(chan==0x0002)&(x=="omni8")&&return error("Incorrect channel count")
!(chan==0x0002)&(x=="g")&&return error("Incorrect channel count")
!(chan==0x0002)&(x=="bin")&&return error("Incorrect channel count")
!(chan==0x0001)&(x=="")&&error("Unspecified Format")
end



"""
# Acoustic Load

acoustic_load(path,format) - Loads file for processing by other functions in Acoustics.jl

**path** - The location of the supported audio file. Must either be a fullpath or relative to the current working directory. to find the current working directory type pwd(). Look into julia shell for more information
**format** - This used to flag the channel order in a multichannel input for certian functions.

## Example

`julia>` a=acoustic_load("Center Hi Sweep-5 impulse_short.wav")

The wav file has been loaded into the variable a

`julia>` a.samples

this is an array of individual samples from the loaded files

`julia>` a.samplerate

this is a float contianing the samplerate


`julia>` a.name

this is a string of the file name of imported file

`julia>` a.channels

this is a unsigned integer

"""
function acoustic_load(path::String,format::String="")

	i=length(path)

	p_end=length(path)

	while !(path[i]=='.')

		p_end=i

		i-=1

	end

	p_end=p_end-2

	i=length(path)


	p_beg=1

	if((path[1]=='/'))

		while (!(path[i]=='/'))

			p_beg=i

			i-=1

		end
		temp=wavread(path)

	else

		path_h=pwd()
		head_path=string(path_h,'/',path)
		temp=wavread(head_path)


	end


	if ("wav"==path[(length(path)-2):end])||("WAV"==path[(length(path)-2):end])
	samples=temp[1]
	samplerate=temp[2]
	else


	end

	sizez=size(samples)
	l_samples=sizez[1]
	chan=convert(UInt16,sizez[2])
	format=format_parser(format,chan)
	name=path[p_beg:p_end]

	return Acoustic(samples,samplerate,name,chan,format,l_samples)


end
#=
function general(source,weighting="z",band="b" ;s=1)

	samplerate=source.samplerate

	l=source.l_samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


	if (band=="b")||(band=="B")

		return source

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]



	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]



	results=[filt(digitalfilter(bands[1],Butterworth(4)),source)]


	i=2

	while length(bands)>=i


		results=vcat(results,[filt(digitalfilter(bands[i],Butterworth(4)),source)])

		i+=1

	end



	return hcat(center,results)

	else

	end

end
=#

"""
# Leq -
`Leq(source,weighting,bands=0)`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
C is known as Clarity it is the balance between early and late eneregy in an impulse expressed in Decibels (dB). Rooms with a positive C value will have greater percieved definition or clarity in reverberance.The Just Noticeable Diffrence (JND) for clarity metrics is 1 dB.

### Recomendations
* Use 50ms for rooms that will be used for music
* Use 80ms for rooms that will be used for speech

Leq

See ISO-3382 for more information
"""
function Leq(source;weighting::String="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time]))/sum(abs2.(x[time:end])))

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"

			writecsv2("Leq_"*source.name*".csv",vcat(["Frequency (Hz)" "L"*weighting*"eq(dB)"],hcat(center,results)))

		else

		end

	end

end

"""
# C - Clarity
`C(source,time,weighting,bands=0)`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
C is known as Clarity it is the balance between early and late eneregy in an impulse expressed in Decibels (dB). Rooms with a positive C value will have greater percieved definition or clarity in reverberance.The Just Noticeable Diffrence (JND) for clarity metrics is 1 dB.

### Recomendations
* Use 50ms for rooms that will be used for music
* Use 80ms for rooms that will be used for speech



See ISO-3382 for more information
"""
function C(source,time::Number;weighting::String="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	real_time=time
	time=Int((time/1000.0)*source.samplerate)
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end



#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time]))/sum(abs2.(x[time:end])))

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)


		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "C"*string(real_time)*"(Weight="*weighting*") (dB)"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("C"*string(real_time)*"_"*name*".csv",file_output)

		else

		end

	end

end
"""
# D - Definition

`D(source,time,weighting,bands=0)` -> ratio (unitless)

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]


### Explation
D is known as Definition it is the balance between early and late eneregy in an impulse expressed as ratio. Rooms with a D ratio greater than one will have greater percieved definition or clarity in reverberance.The Just Noticeable Diffrence (JND) for definition metrics is 0,05.

### Recomendations
* Use 50ms for rooms that will be used for music
* Use 80ms for rooms that will be used for speech



See ISO-3382 for more information
"""
function D(source,time;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	real_time=time
	time=Int((time/1000.0)*source.samplerate)
	name=source.name
	source=source.samples


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
	f(x)=sum(abs2.(x[1:time]))/sum(abs2.(x[time:end]))

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "D"*string(real_time)*"(Weight="*weighting*")"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("D"*string(real_time)*"_"*name*".csv",file_output)

		else

		end

	end

end

"""
# RT - Reverberation Time

`RT(source,decay,weighting,bands=0)` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Decay** - The amount of level decay from steady state (dB)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]


### Explation
RT is known as Reverberation time it is the measure of decay from steady state to some level.

### Recomendations
* normally measured in T20 & T30


See ISO-3382 for more information

"""
function RT(source,decay;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	l=source.l_samples
	sampl_amount=Int(ceil(samplerate/1000))
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end



	function f(x)
		x=abs2.(x[:,1])

		max=sum(x[:,1])

		#takes only the first channel
		x=reverse((/).(x[:,1],max))


		target=[(10^(-5/10)),(10^(-((decay+5)/10)))]


#hi and low refer to level
		hi_range=1
		lo_range=1
		sampled_y=[1.0]
		sampled_x=[0.0]
		i=l-sampl_amount
		total=2

	if x[l]>target[2]

		return Inf

	else


		while 0<i

			sampled_y=vcat(sampled_y,sum(x[1:i]))

			i-=sampl_amount
		end

			sampled_y=vcat(sampled_y,x[1])


#the -5dB decay point

		total=1
		while (total>=target[1])&&(hi_range<l)
			hi_range+=1
			total=sampled_y[hi_range]

		end
#Starts the search at -5dB point
		lo_range=hi_range

		while (total>=target[2])&&(lo_range<l)
			lo_range+=1
			total=sampled_y[lo_range]
		end

	sequence=LinRange(0,length(sampled_y[hi_range:lo_range])*(sampl_amount/samplerate),length(sampled_y[hi_range:lo_range]))
	std_x=std(sequence)
	std_y=std(10*log.(10,sampled_y[hi_range:lo_range]))
	r=cor(sequence,10*log.(10,sampled_y[hi_range:lo_range]))

		m=r*(std_y/std_x)

			return -60/m

		end

	end

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "T"*string(decay)*" (Weight="*weighting*") (s)"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("T"*string(decay)*"_"*name*".csv",file_output)

		else

		end

	end

end

"""
 EDT - Early Decay Time

`EDT(source,decay,weighting,bands=0)` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]


### Explation
EDT is known as Early Decay Time it is the measure of decay from peak to 10dB down.

### Recomendations
* EDT is important as it correlates percieved reverberance

See ISO-3382 for more information

"""
function EDT(source;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	l=source.l_samples
	sampl_amount=Int(ceil(samplerate/1000))
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end



	function f(x)
		x=abs2.(x[:,1])

		max=sum(x[:,1])

		#takes only the first channel
		x=reverse((/).(x[:,1],max))


		target=[1,10^(-1)]


#hi and low refer to level
		hi_range=1
		lo_range=1
		sampled_y=[1.0]
		sampled_x=[0.0]
		i=l-sampl_amount
		total=2

	if x[l]>target[2]

		return Inf

	else


		while 0<i

			sampled_y=vcat(sampled_y,sum(x[1:i]))

			i-=sampl_amount
		end

			sampled_y=vcat(sampled_y,x[1])


#the -5dB decay point

		total=1
		while (total>=target[1])&&(hi_range<l)
			hi_range+=1
			total=sampled_y[hi_range]

		end
#Starts the search at -5dB point
		lo_range=hi_range

		while (total>=target[2])&&(lo_range<l)
			lo_range+=1
			total=sampled_y[lo_range]
		end

	sequence=LinRange(0,length(sampled_y[hi_range:lo_range])*(sampl_amount/samplerate),length(sampled_y[hi_range:lo_range]))
	std_x=std(sequence)
	std_y=std(10*log.(10,sampled_y[hi_range:lo_range]))
	r=cor(sequence,10*log.(10,sampled_y[hi_range:lo_range]))

		m=r*(std_y/std_x)

			return -60/m

		end

	end

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "Early Decay Time"*" (Weight="*weighting*") (s)"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("EDT_"*name*".csv",file_output)

		else

		end

	end

end

"""
# Ts- Centre Time
`Ts(source,weighting,band=0)`-> milliseconds (ms)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
Ts is the time centre which is centre of gravity of the squared impulse repose. It finds to find the center of energy of an impulse response. It is reported in milliseconds. Just Noticeable Diffrence (JND) for the Centre Time (Ts) metrics is 10ms.

### Recomendations
* Avoids the discrete division with C & D
* Good for finding strongest reflection location



See ISO-3382 for more information
"""
function Ts(source;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	l=source.l_samples
	t=LinRange(0,l/samplerate,l)
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end


	function f(x)
		x=abs2.(x)
		top=(*).(x,t)

		return (sum(top)/sum(x))*1000

	end


	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "Time Centre "*" (Weight="*weighting*") (ms)"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("Ts_"*name*".csv",file_output)

		else

		end

	end

end

"""
# ST_early- Early Support
`ST_early(source,weighting="z",bands=0)`-> ratio (dB)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
ST_early is the ratio of the reflectioned energy relative to the direct energy in the first 200th of a second to the first 10th of a second.

### Recomendations
* Is a stage derived measurement is describes of sound of players in an orchestra
* Describes the ease of hearing other members of the orchestra
* The start time is based on the arrival of the direct sound



See ISO-3382 for more information
"""
function ST_early(source;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end


time_1=Int(ceil(0.01*samplerate))
time_2=Int(ceil(0.02*samplerate))
time_3=Int(ceil(0.1*samplerate))

#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time_1]))/sum(abs2.(x[time_2:time_3])))

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "Early Support (dB"*weighting*")"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("ST_early_"*name*".csv",file_output)

		else

		end

	end

end

"""
# ST_late- Late Support
`ST_early(source,weighting="z",bands=0)`-> ratio (dB)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
ST_early is the ratio of the reflectioned energy relative to the direct energy in the first 10th of a second and 1 second.

### Recomendations
* Is a stage derived measurement is describes of sound of players in an orchestra
* Describes the percieved reverberance of a room
* The start time is based on the arrival of the direct sound



See ISO-3382 for more information
"""
function ST_late(source;weighting="z",bands::Int64=0,output="shell")

	samplerate=source.samplerate
	name=source.name
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")
		ifil=a_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="c")||(weighting=="C")
		ifil=c_wtd(samplerate)
		source=filt(ifil,source)
	elseif (weighting=="ccir")||(weighting=="CCIR")
		ifil=ccir(samplerate)
		source=filt(ifil,source)
	else
		return print("Weighting undefined")

	end


time_1=Int(ceil(0.01*samplerate))
time_2=Int(ceil(0.1*samplerate))
time_3=Int(ceil(samplerate))

#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time_1]))/sum(abs2.(x[time_2:time_3])))

	if (bands==0)||(bands==0)

		return f(source)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source)),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"
			header=["Frequency (Hz)" "Late Support (dB"*weighting*")"]
			file_output=vcat(header,string.(hcat(center,results)))

			writecsv2("ST_late_"*name*".csv",file_output)

		else

		end

	end

end

"""
# IACC- Inter-Aural Cross Correlation Coefficients
`IACC(source,weighting="z",bands=0 ;s=1)`-> time (ms)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
IACC is the point of maximum cross correlation between the left and right ears.

### Recomendations
* Correlates well with percieved spaciousness of a concert hall
* Measures the diffrence of sound arrive at each ear



See ISO-3382 for more information
"""
function IACC(source;weighting="z",bands::Int64=0)

	samplerate=source.samplerate
	name=source.name
	source=source.samples
	left=source[:,1]
	right=source[:,2]

	#x-left and y-right
	function f(x,y)
		x=fft(x)
		y=fft(y)
		top=(*).(x,conj.(y))
		bottom=sqrt(sum((^).(x,2))*sum((^).(y,2)))

		return maximum(real.((/).(top,bottom)))
	end

	return f(left,right)
end


function G(source,weighting="z",bands::Int64=0 ;s=1)

	samplerate=source.samplerate
	l=source.l_samples
	source=source.samples
	l=source[:,1]
	l_10=source[:,2]

	#x is direct and y is 10m

	function f(x,y)
		x=sum(abs2.(x))
		y=sum(abs2.(y))

		return log10(x/y)
	end

	return f(l,l_10)

end
"""
# J_LF- Early Lateral Fraction
`J_LF(source,weighting,band)`-> ratio

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
J_LF is the ratio between a figure-8 microphone microphone null pointed at the sound source and omnidirectional microphone at a measurement position in the first 80ms. It relates to the percieved width of a sound source.

### Recomendations
* It can be obtained from numerous methods
* It relates to the percieved width of a rooom

See ISO-3382 for more information
"""
function J_LF(source;weighting="z",bands::Int64=0)

	samplerate=source.samplerate
	name=source.name
	source=source.samples
	l_o=source[:,1]
	l_8=source[:,2]
	time_1=Int(ceil(samplerate*0.005))
	time_2=Int(ceil(samplerate*0.08))

#x is the omni & y is the figure 8
	function f(x,y)
		x=sum(abs2.(x[1:time_2]))
		y=sum(abs2.(y[time_1:time_2]))

		return (y/x)
	end

	return f(l_o,l_g)

end

"""
# L_j- Late Lateral Fraction
`L_j(source,weighting,band)`-> ratio (dB)

* **Source** - The audio file loaded by acoustic_load
* **Weighting** - The frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]

### Explation
L_j is the ratio between a figure-8 microphone microphone null pointed at the sound source and omnidirectional microphone at a measurement position for the length of the whole impulse. It relates to the percieved width of a sound source.

### Recomendations
* It can be obtained from numerous methods
* It relates to the percieved spaciousness of the room

See ISO-3382 for more information
"""
function L_j(source;weighting="z",bands::Int64=0)

	samplerate=source.samplerate
	name=source.name
	source=source.samples
	l_o=source[:,1]
	l_8=source[:,2]
	time_1=Int(ceil(samplerate*0.08))

#x is the omni & y is the figure 8
	function f(x,y)
		x=sum(abs2.(x))
		y=sum(abs2.(y[time_1:end]))

		return (y/x)
	end

	return f(l_o,l_g)

end

function swup(time,duration,f_1,f_2)
#time seconds
#duraction seconds
#f_1<f_2
	K=(duration*2*pi*f_1)/log((f_2)/(f_1))
	L=duration/log((f_2)/(f_1))
	return sin(K*(exp(time/L)-1))
end

#inverse generation
function iswup(time,duration,f_1,f_2)
#time seconds
#duraction seconds
#f_1<f_2
	L=duration/log((f_2)/(f_1))
	m=exp(-(time)/L)
	return m
end


"""
# sweep- Logarithmic Sine Sweep
`sweep(duration,silence_duration,f_1,f_2,samplerate,alpha=0.01)`-> sweep & inverse sweep in current working directory

* **Duration** - The duration of the sine sweep in seconds
* **Silence Duration** - the duration of the silence following the sweep
* **F_1** - The start frequency (Hz) of the sweep
* **F_2** - The end frequency (Hz) of the sweep
* **Samplerate** - The samplerate (Hz) of the sine sweep
* **α** -  The mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window. This parameter controls the tukey window.

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. If alpha is zero a click will be heard.

### Recomendations
* The default value for α will provide a perfectly pop free experience but, cause a loss of high frequencies above 94.9125% of f_2
* Use this to capture multiple channel impulse values.
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur due the change in frequency


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep(duration,silence_duration,f_1,f_2,samplerate,α=0.01)
	#duration in seconds
	values=string.([duration,silence_duration,f_1,f_2,α,samplerate])
	sequence=LinRange(0.0,duration,Int(floor(samplerate*duration)))
	window=tukey(Int(floor(samplerate*duration)),α)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(Int(floor(samplerate*silence_duration)))
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep-Duration_"*values[1]*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse_"*"Sweep-Duration_"*values[1]*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
end


"""
# sweep_target- adaptive Logarithmic Sine Sweep
`sweep_target(duration,silence_duration,f_1,f_2,samplerate,alpha=0.0003)`-> sweep & inverse sweep in current working directory

* **Duration** - This is the targeted duration of the sine sweep in seconds
* **Silence Duration** - The duration of the silence following the sweep
* **F_1** - The start frequency (Hz) of the sweep
* **F_2** - The end frequency (Hz) of the sweep
* **Samplerate** - The samplerate (Hz) of the sine sweep
* **α** -  The mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window. This parameter controls the tukey window.

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. If alpha is zero a click will be heard. This function differs from sweep in that it tries to find an optimal duration that will cause the sweep function to end on zero. Allowing smaller α values which will reduce high frequency loss in the impulse. This function will always generate a duration smaller than the input duration.

### Recomendations
* The default value for α will provide a perfectly pop free experience but, cause a loss of high frequencies above 99.825% of f_2
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur due the change in frequency


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep_target(duration,silence_duration,f_1,f_2,samplerate,α=0.0003)
	#duration in seconds
	values=string.([duration,silence_duration,f_1,f_2,α,samplerate])
	r=f_1*((exp(-log(f_2/f_1))-1)/log(f_2/f_1))
	n_target=ceil(r*duration)
	duration_target=n_target/r
	sequence=LinRange(0.0,duration_target,Int(floor(samplerate*duration_target)))
	window=tukey(Int(floor(samplerate*duration_target)),α)
	sweep=(*).(swup.(sequence,duration_target,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration_target,f_1,f_2),iswup.(sequence,duration_target,f_1,f_2),window)
	silence=zeros(Int(floor(samplerate*silence_duration)))
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep-Duration_"*string(duration_target)*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse_"*"Sweep-Duration_"*string(duration_target)*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
end



"""
# Deconvolve - Generates impulse responses from sine sweeps
`deconvolve(inverse,measured;title="",output="file")`-> N-channel Wav file impulse

* **Inverse** - This file should be monophonic inverse sweep file.
* **Measured** - This should be the N channel capture sweep
* **Title** - If you want to name the file something besides the name of the measured sweep with impulse appended
* **Output** - The "file" argument saves the impulse to a wave file. The "acoustic_load" argument allows you to store the results to a variable. file is the default.

### Explation
Deconvolve converts a measured sweep into an impulse response using Logarithmic Sine Sweep Method.

### Recomendations
* Audio must be loaded with acoustic_load
* The first argument must be a mono file inverse sweep
* The measured file must be identical number of samples per channel as inverse file
* Ardour will make sample accurate edits
* Audacity will not make sample accurate edits
* Type pwd() to find the where impulses are saved


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function deconvolve(inverse,measured;title::String="",output="file")

	l=measured.l_samples
	title=String(title)
	samplerate=measured.samplerate
	colmn=measured.channels
	format=measured.format

	if length(title)==0

		title=measured.name
	else
		title=title

	end


	if measured.l_samples<inverse.l_samples

		error("Measured sweep has less samples than the generated sweep")

	else measured.l_samples==inverse.l_samples

	inverse=rfft(inverse.samples)
	measured=rfft(measured.samples)
	imp=(*).(measured,inverse)
	rimp=irfft(imp,l)

	if output=="file"
		return wavwrite(rimp,title*"-impulse.wav",Fs=samplerate)
	elseif output=="acoustic_load"
		return Acoustic(rimp,samplerate,title*"-impulse.wav",colmn,format,l)
	end

	end

end




end # module
