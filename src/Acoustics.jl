module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,Reexport,DataFrames

export Leq,C,RT,D,Ts,sweep,deconvolve,EDT,acoustic_load,ST_late,ST_early,IACC,G,sweep_target,peak_loc,seq_create,seq_deconvolve,parseval_crop,gain,Leqm

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
	samples::Array{<:AbstractFloat}
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

function acoustic_loader(path::String,format::String="")

	if isdir(path)
		loc=[]
		for file in readdir(path)
			#remove hidden files
			if ('.'==file[1])||(isdir(joinpath(path,file)))

			else

				loc=vcat(loc,[joinpath(path,file)])
			end
		end

		l_loc=length(loc)
		swpz=Vector{Acoustic}(undef,length(loc))

		for index in 1:l_loc
			t_path=loc[index]
			i=length(t_path)
			p_end=length(t_path)

			#Get the name Automted
			while !(t_path[i]=='.')
				p_end=i
				i-=1
			end

			p_end=p_end-2
			i=length(path)

			p_beg=1

			if((t_path[1]=='/'))

				while (!(t_path[i]=='/'))
					p_beg=i
					i-=1
				end



			else
				path_h=pwd()
				head_path=string(path_h,'/',t_path)
				temp=wavread(head_path)
			end
				temp=wavread(t_path,format="native")
				if typeof(temp[1])==Matrix{Int32}
					temp=wavread(t_path)
					if ("wav"==t_path[(length(path)-2):end])||("WAV"==t_path[(length(t_path)-2):end])
						samples=temp[1]
						samplerate=temp[2]
					else
					end

				else
					temp=wavread(t_path)
					if ("wav"==t_path[(length(path)-2):end])||("WAV"==t_path[(length(t_path)-2):end])
						samples=Float32.(temp[1])
						samplerate=temp[2]
					else
					end

				end

			if ("wav"==t_path[(length(path)-2):end])||("WAV"==t_path[(length(t_path)-2):end])
				samples=temp[1]
				samplerate=temp[2]
			else


			end

			sizez=size(samples)
			l_samples=sizez[1]
			chan=convert(UInt16,sizez[2])
			format=format_parser(format,chan)
			name=t_path[p_beg:p_end]

			swpz[i]=Acoustic(samples,samplerate,name,chan,format,l_samples)
		end

		return swpz

	else
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
			temp=wavread(path,format="native")
			if typeof(temp[1])==Matrix{Int32}
				temp=wavread(path)
				if ("wav"==path[(length(path)-2):end])||("WAV"==path[(length(path)-2):end])
					samples=temp[1]
					samplerate=temp[2]
				else
				end

			else
				temp=wavread(path)
				if ("wav"==path[(length(path)-2):end])||("WAV"==path[(length(path)-2):end])
					samples=Float32.(temp[1])
					samplerate=temp[2]
				else
				end

			end

		else
			path_h=pwd()
			head_path=string(path_h,'/',path)
			temp=wavread(head_path,format="native")
			if typeof(temp[1])==Matrix{Int32}
				temp=wavread(head_path)
				if ("wav"==head_path[(length(head_path)-2):end])||("WAV"==head_path[(length(head_path)-2):end])
					samples=temp[1]
					samplerate=temp[2]
				else
				end

			else
				temp=wavread(head_path)
				if ("wav"==head_path[(length(head_path)-2):end])||("WAV"==head_path[(length(head_path)-2):end])
					samples=Float32.(temp[1])
					samplerate=temp[2]
				else
				end

			end
		end


		sizez=size(samples)
		l_samples=sizez[1]
		chan=convert(UInt16,sizez[2])
		format=format_parser(format,chan)
		name=path[p_beg:p_end]

		return Acoustic(samples,samplerate,name,chan,format,l_samples)


	end

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

	if isdir(path)
		loc=[]
		for file in readdir(path)
			if ('.'==file[1])||(isdir(joinpath(path,file)))

			else
				loc=vcat(loc,[joinpath(path,file)])
			end
		end

		l_loc=length(loc)
		swpz=Vector{Acoustic}(undef,length(loc))

		for i in 1:l_loc
			swpz[i]=acoustic_loader(loc[i],format)
		end

		return swpz

	else
		return acoustic_loader(path,format)

	end

end



"""
# Leq - time-averaged sound level or equivalent continuous sound level
`Leq(source,weighting,bands=0,ref=1,output="shell")`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **ref** - This is the reference level this will scale your recording to the proper level
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

### Explation
This is the average sound level over the whole file

### Recomendations
* Use this to find the signal to noise ratio
* Use this to find longterm noise level

See ANSI/ASA S1.1-2013 for more information
"""
function Leq(source;weighting::String="z",bands::Int64=0,ref=1,output="shell")

	samplerate=source.samplerate
	l=source.l_samples
	source=source.samples
	#chan=source.channels

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
function f(x,n,p_0)

sig=abs2.(x)
power=sum(sig)
powerdb=10*log10(power)
norm=10*log10(n*p_0^2)

	return powerdb-norm

end

	if (bands==0)||(bands==0)

		return f(source,l,ref)

	else (bands>0)

	gen_bands=generateband(bands,samplerate)
	bands=gen_bands[1]
	center=gen_bands[2]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(3)),source),l,ref),bands)

		if output=="shell"

			return hcat(center,results)

		elseif output=="file"

			writecsv2("Lt_"*source.name*".csv",vcat(["Frequency (Hz)" "L"*weighting*"eq(dB)"],hcat(center,results)))

		else

		end

	end

end

"""
# Leqm - time-averaged sound level or equivalent continuous sound level multi channel
`Leq(source,weighting,bands=0,ref=1,output="shell",combination)`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **ref** - This is the reference level this will scale your recording to the proper level
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.
* **Combination** - when "individual" it returns a text file with the calculation for each channel and "averaged" is the average of all the channels


### Explation
This is the average sound level over the whole file

### Recomendations
* Use this to find the signal to noise ratio
* Use this to find longterm noise level

See ANSI/ASA S1.1-2013 for more information
"""
function Leqm end
function Leqm(source;weighting::String="z",bands::Int=0,ref=1,output="shell",combination::String="individual")

	if (weighting=="z")||(weighting=="Z")
	elseif (weighting=="a")||(weighting=="A")
	elseif (weighting=="c")||(weighting=="C")
	elseif (weighting=="ccir")||(weighting=="CCIR")
	else
		error("Weighting undefined")
	end

	samplerate=source.samplerate
	l=source.l_samples
	samples=source.samples
	chan=source.channels

	#f(x) is the defined function
	function f(x,n,p_0)

	sig=abs2.(x)
	power=sum(sig)
	powerdb=10*log10(power)
	norm=10*log10(n*p_0^2)

		return powerdb-norm

	end

	results=Dict()

	for channels in 1:chan
		if (weighting=="z")||(weighting=="Z")
			filtered=samples[:,channels]
		elseif (weighting=="a")||(weighting=="A")
			ifil=a_wtd(samplerate)
			filtered=filt(ifil,samples[:,channels])
		elseif (weighting=="c")||(weighting=="C")
			ifil=c_wtd(samplerate)
			filtered=filt(ifil,samples[:,channels])
		elseif (weighting=="ccir")||(weighting=="CCIR")
			ifil=ccir(samplerate)
			filtered=filt(ifil,samples[:,channels])
		else
			error("Weighting undefined")
		end
		value_name="Channel "*string(channels)
		temp=f(filtered,l,ref)
		merge!(results,Dict([(value_name,temp)]))

	end
	#temporary dataframe array
	t_results=DataFrame(results)
	if combination=="individual"
		results=DataFrame()
		results.Bands=["Broadband"]
		results=t_results
	elseif combination=="averaged"
		indi=t_results[:,1]
		for i in 2:chan
				indi=hcat(indi,t_results[:,i])
		end
		indi=db2amp.(indi)
		avg=mean(indi)
		variance=var(indi)
		results=DataFrame()
		results.Bands=["Broadband"]
		results."Averaged (dB)"=[amp2db(avg)]
		results."Standard Deviation (dB)"=[powerdb(variance)]

	else
	end



	return results
end

function Leqm(source;weighting::String="z",bands::Int=0,ref=1,output="shell",combination::String="individual")

	if (weighting=="z")||(weighting=="Z")
	elseif (weighting=="a")||(weighting=="A")
	elseif (weighting=="c")||(weighting=="C")
	elseif (weighting=="ccir")||(weighting=="CCIR")
	else
		error("Weighting undefined")
	end

	samplerate=source.samplerate
	l=source.l_samples
	samples=source.samples
	chan=source.channels

	#f(x) is the defined function
	function f(x,n,p_0)

	sig=abs2.(x)
	power=sum(sig)
	powerdb=10*log10(power)
	norm=10*log10(n*p_0^2)

		return powerdb-norm

	end

	results=Dict()

	for channels in 1:chan
		if (weighting=="z")||(weighting=="Z")
			filtered=samples[:,channels]
		elseif (weighting=="a")||(weighting=="A")
			ifil=a_wtd(samplerate)
			filtered=filt(ifil,samples[:,channels])
		elseif (weighting=="c")||(weighting=="C")
			ifil=c_wtd(samplerate)
			filtered=filt(ifil,samples[:,channels])
		elseif (weighting=="ccir")||(weighting=="CCIR")
			ifil=ccir(samplerate)
			filtered=filt(ifil,samples[:,channels])
		else
			error("Weighting undefined")
		end
		value_name="Channel "*string(channels)
		temp=f(filtered,l,ref)
		merge!(results,Dict([(value_name,temp)]))

	end
	#temporary dataframe array
	t_results=DataFrame(results)
	if combination=="individual"
		results=DataFrame()
		results.Bands=["Broadband"]
		results=t_results
	elseif combination=="averaged"
		indi=t_results[:,1]
		for i in 2:chan
				indi=hcat(indi,t_results[:,i])
		end
		indi=db2amp.(indi)
		avg=mean(indi)
		variance=var(indi)
		results=DataFrame()
		results.Bands=["Broadband"]
		results."Averaged (dB)"=[amp2db(avg)]
		results."Standard Deviation (dB)"=[powerdb(variance)]

	else
	end



	return results
end

"""
# C - Clarity
`C(source,time,weighting,bands=0,output="shell")`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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

`D(source,time,weighting,bands=0,output="shell")` -> ratio (unitless)

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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

`RT(source,decay,weighting,bands=0,output="shell")` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Decay** - The amount of level decay from steady state (dB)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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

`EDT(source,decay,weighting,bands=0,output="shell")` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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
# Ts - Centre Time
`Ts(source,weighting,band=0,output="shell")`-> milliseconds (ms)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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
# ST_early - Early Support
`ST_early(source,weighting="z",bands=0,output="shell")`-> ratio (dB)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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
# ST_late - Late Support
`ST_early(source,weighting="z",bands=0)`-> ratio (dB)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **Output** - can be either "shell" which returns the results to the terminal or "file" which will put the results in a file.

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
# IACC - Inter-Aural Cross Correlation Coefficients
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
# J_LF - Early Lateral Fraction
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
# L_j - Late Lateral Fraction
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
	K=(duration*2*f_1)/log((f_2)/(f_1))
	L=duration/log((f_2)/(f_1))
	return sinpi(K*(exp(time/L)-1))
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
# sweep - Logarithmic Sine Sweep
`sweep(duration,silence_duration,f_1,f_2,samplerate,α=0.01;inv_norm="p")`-> sweep & inverse sweep in current working directory

* **Duration** - The duration of the sine sweep in seconds
* **Silence Duration** - the duration of the silence following the sweep
* **F_1** - The start frequency (Hz) of the sweep
* **F_2** - The end frequency (Hz) of the sweep
* **Samplerate** - The samplerate (Hz) of the sine sweep
* **α** -  The mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window. This parameter controls the tukey window.
* **inv_norm** - This is the normalization gain applied to the inverse sweep. The default value of "p" will normalize to the peak passband amplitude. "a" will normalize to the average passband amplitude. "e" will normalize to the total energy. All other values will cause the signal to be unnormalized.

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. If alpha is zero a click will be heard.

### Recomendations
* The default value for α will provide a perfectly pop free experience but, cause a loss of high frequencies above 94.9125% of f_2
* Use this to capture multiple channel impulse values.
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur due the change in frequency
* Leave the inverse sweep normalization at either "p" or "a" the unnormalized sweep applies a pretty large gain in the passband and which will cause impulses to clip.
* If measuring the absolute energy is important then use "e" other wise it cause your signal to be extremely quiet for no benefit.


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep(duration,silence_duration,f_1,f_2,samplerate,α=0.01;inv_norm="p")
	#duration in seconds
	values=string.([duration,silence_duration,f_1,f_2,α,samplerate])
	sequence=range(0.0,duration,length=Int(floor(samplerate*duration)))
	window=tukey(Int(floor(samplerate*duration)),α)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(Int(floor(samplerate*silence_duration)))
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)

	#Amplitude Normalization
	sweep_length=length(sweep)
	padding_length=nextfastfft(2*sweep_length)-sweep_length
	padding=zeros(padding_length)
	sweep_fft=vcat(sweep,padding)
	isweep_fft=vcat(isweep,padding)
	sweep_fft=fft(sweep_fft)
	isweep_fft=fft(isweep_fft)
	conv=(*).(sweep_fft,isweep_fft)
	length_conv=padding_length+sweep_length
	bin=LinRange(0,(length_conv-1),length_conv)
	norm_bin=(/).(bin,length_conv)
	freq_bin=(*).(samplerate,norm_bin)
	low_index=0
	hi_index=0

	if inv_norm=="p"
		for n in 1:1:Int(floor(0.5*length_conv)+1)

			if freq_bin[n]==f_1
				low_index=n
			else
				if 1<n
					if freq_bin[n-1]<f_1<freq_bin[n+1]
						low_index=n
					end
				end

			end

			if f_2==(samplerate*0.5)
				hi_index=Int(floor(0.5*length_conv)+1)
			else
				if freq_bin[n]==f_2
					hi_index=n
				else
					if 1<n
						if freq_bin[n-1]<f_2<freq_bin[n+1]
							hi_index=n
						end
					end
				end
			end

		end

		mag_conv=abs.(conv[low_index:hi_index])
		mag_conv=maximum(mag_conv)
		avg_amp=mag_conv
		println("Low End: ",freq_bin[low_index],", High End: ",freq_bin[hi_index]," Peak Passband Amplitude: ",avg_amp)
		isweep=(/).(isweep,avg_amp)
	elseif inv_norm=="a"
		for n in 1:1:Int(floor(0.5*length_conv)+1)

			if freq_bin[n]==f_1
				low_index=n
			else
				if 1<n
					if freq_bin[n-1]<f_1<freq_bin[n+1]
						low_index=n
					end
				end

			end

			if f_2==(samplerate*0.5)
				hi_index=Int(floor(0.5*length_conv)+1)
			else
				if freq_bin[n]==f_2
					hi_index=n
				else
					if 1<n
						if freq_bin[n-1]<f_2<freq_bin[n+1]
							hi_index=n
						end
					end
				end
			end

		end

		mag_conv=abs.(conv[low_index:hi_index])
		mag_conv=sum(mag_conv)
		avg_amp=mag_conv/(bin[hi_index]-bin[low_index])
		println("Low End: ",freq_bin[low_index],", High End: ",freq_bin[hi_index]," Average Passband Amplitude: ",avg_amp)
		isweep=(/).(isweep,avg_amp)
	elseif inv_norm=="e"
		energy=abs2.(conv)
		energy=sum(energy)/(2*pi)
		isweep=(/).(isweep,energy)
		println("Total Energy: ",energy)

	else
		println("Unnormalized")

	end

	#work on the output
	wavwrite(sweep,"Sweep-Duration_"*values[1]*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse_"*"Sweep-Duration_"*values[1]*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
end



"""
# sweep_target - adaptive Logarithmic Sine Sweep
`sweep_target(duration,silence_duration,f_1,f_2,samplerate,α=0.0003;inv_norm="p")`-> sweep & inverse sweep in current working directory

* **Duration** - This is the targeted duration of the sine sweep in seconds
* **Silence Duration** - The duration of the silence following the sweep
* **F_1** - The start frequency (Hz) of the sweep
* **F_2** - The end frequency (Hz) of the sweep
* **Samplerate** - The samplerate (Hz) of the sine sweep
* **α** -  The mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window. This parameter controls the tukey window.
* **inv_norm** - This is the normalization gain applied to the inverse sweep. The default value of "p" will normalize to the peak passband amplitude. "a" will normalize to the average passband amplitude. "e" will normalize to the total energy. All other values will cause the signal to be unnormalized.

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. If alpha is zero a click will be heard. This function differs from sweep in that it tries to find an optimal duration that will cause the sweep function to end on zero. Allowing smaller α values which will reduce high frequency loss in the impulse. This function will always generate a duration smaller than the input duration.

### Recomendations
* The default value for α will provide a perfectly pop free experience but, cause a loss of high frequencies above 99.825% of f_2
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur due the change in frequency
* Leave the inverse sweep normalization at either "p" or "a" the unnormalized sweep applies a pretty large gain in the passband and which will cause impulses to clip.
* If measuring the absolute energy is important then use "e" other wise it cause your signal to be extremely quiet for no benefit.


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep_target(duration,silence_duration,f_1,f_2,samplerate,α=0.01;inv_norm="p")
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

	#Amplitude Normalization
	sweep_length=length(sweep)
	padding_length=nextfastfft(2*sweep_length)-sweep_length
	padding=zeros(padding_length)
	sweep_fft=vcat(sweep,padding)
	isweep_fft=vcat(isweep,padding)
	sweep_fft=fft(sweep_fft)
	isweep_fft=fft(isweep_fft)
	conv=(*).(sweep_fft,isweep_fft)
	length_conv=padding_length+sweep_length
	bin=LinRange(0,(length_conv-1),length_conv)
	norm_bin=(/).(bin,length_conv)
	freq_bin=(*).(samplerate,norm_bin)
	low_index=0
	hi_index=0

	if inv_norm=="p"
		for n in 1:1:Int(floor(0.5*length_conv)+1)

			if freq_bin[n]==f_1
				low_index=n
			else
				if 1<n
					if freq_bin[n-1]<f_1<freq_bin[n+1]
						low_index=n
					end
				end

			end

			if f_2==(samplerate*0.5)
				hi_index=Int(floor(0.5*length_conv)+1)
			else
				if freq_bin[n]==f_2
					hi_index=n
				else
					if 1<n
						if freq_bin[n-1]<f_2<freq_bin[n+1]
							hi_index=n
						end
					end
				end
			end

		end

		mag_conv=abs.(conv[low_index:hi_index])
		mag_conv=maximum(mag_conv)
		avg_amp=mag_conv
		println("Low End: ",freq_bin[low_index],", High End: ",freq_bin[hi_index]," Peak Passband Amplitude: ",avg_amp)
		isweep=(/).(isweep,avg_amp)
	elseif inv_norm=="a"
		for n in 1:1:Int(floor(0.5*length_conv)+1)

			if freq_bin[n]==f_1
				low_index=n
			else
				if 1<n
					if freq_bin[n-1]<f_1<freq_bin[n+1]
						low_index=n
					end
				end

			end

			if f_2==(samplerate*0.5)
				hi_index=Int(floor(0.5*length_conv)+1)
			else
				if freq_bin[n]==f_2
					hi_index=n
				else
					if 1<n
						if freq_bin[n-1]<f_2<freq_bin[n+1]
							hi_index=n
						end
					end
				end
			end

		end

		mag_conv=abs.(conv[low_index:hi_index])
		mag_conv=sum(mag_conv)
		avg_amp=mag_conv/(bin[hi_index]-bin[low_index])
		println("Low End: ",freq_bin[low_index],", High End: ",freq_bin[hi_index]," Average Passband Amplitude: ",avg_amp)
		isweep=(/).(isweep,avg_amp)
	elseif inv_norm=="e"
		energy=abs2.(conv)
		energy=sum(energy)/(2*pi)
		isweep=(/).(isweep,energy)
		println("Total Energy: ",energy)

	else
		println("Unnormalized")

	end

	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep-Duration_"*string(duration_target)*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse_"*"Sweep-Duration_"*string(duration_target)*"_Silence_"*values[2]*"_Low-Frequency_"*values[3]*"_High-Frequency_"*values[4]*"_Alpha_"*values[5]*"_Fs_"*values[6]*".wav",Fs=samplerate)
end



"""
# Deconvolve - Generates impulse responses from sine sweeps
`deconvolve(inverse,measured;title::String="",norm::String="u",norm_o=1,lp="h",output="file")`-> N-channel Wav file impulse

* **Inverse** - This file should be monophonic inverse sweep file.
* **Measured** - This should be the N channel captured sweep
* **Title** - (optional) If you want to name the file something besides the name of the measured sweep with impulse appended
* **Output** - (optional)The "file" argument saves the impulse to a wave file. The "acoustic_load" argument allows you to store the results to a variable. file is the default.
* **norm** - (optional) Controls the final gain applied to the impulse with strings.l = number of samples, n = peak amplitude, u = unnormalized, o = user defined normalized with with norm_0 setting the inverse amplitude
* **norm_o** - (optional) This allow is the custom gain input. This is the level the signal is divided by norm_o so (final impulse)/norm_o. norm_o is in amplitude not decbels.
* **lp** - 1 -full doubled sided impulse & lp>1- cuts from that sample to the end
### Explation
Deconvolve converts a measured sweep into an impulse response using Logarithmic Sine Sweep Method.

### Recomendations
* Audio must be loaded with acoustic_load
* The first argument must be a mono file inverse sweep
* Ardour will make sample accurate edits
* Audacity will not make sample accurate edits try to make measure sweep slighly longer
* Type pwd() to find the where impulses are saved
* Do not use unnormalized as the signal will clip
* Using other normalization you like peak amplitude will destroy the amplitude relationship between impulses.

See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function  deconvolve(inverse,measured;title::String="",norm::String="u",norm_o=1,lp::Int=1,output="file")

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


	if (((measured.l_samples-inverse.l_samples)/measured.l_samples)<-0.01)

		error("The Measure sweep has a sample length difference greater than 1%. Try re-editing your measured sweep")

	else measured.l_samples==inverse.l_samples

	#find the next small product
	padnum=nextprod([2,3,5,7,11,13,17],2*l)-l

	inver_zeros=zeros(typeof(inverse.samples[1]),(nextprod([2,3,5,7,11,13,17],2*l)-inverse.l_samples),1)
	inver_pad=vcat(inverse.samples,inver_zeros)
	inverse=rfft(inver_pad)

	mea_zeros=zeros(typeof(measured.samples[1]),padnum,colmn)
	mea_pad=vcat(measured.samples,mea_zeros)
	measured=rfft(mea_pad)
	imp=(*).(measured,inverse)
	rimp=irfft(imp,padnum+l)

	if lp>0

		rimp=rimp[lp:end,:]
	else
		error("you input a negative value it cannot not be cut from here")
	end


	if norm=="l"
	#normalized by number of samples
	normalizer=l
	rimp=(/).(rimp,normalizer)
	elseif norm=="n"
	#Peak amplitude normalized
	normalizer=maximum(abs.(rimp))
	rimp=(/).(rimp,normalizer)
	elseif norm=="u"
	#unnormalized
	else norm=="o"
	#other
	normalizer=norm_o
	rimp=(/).(rimp,normalizer)

	end


	if output=="file"
		return wavwrite(rimp,title*"-impulse.wav",Fs=samplerate)
	elseif output=="acoustic_load"
		return Acoustic(rimp,samplerate,title*"-impulse.wav",colmn,format,l)
	end

	end

end

"""
# Peak_loc - gets index of maximum sample value
`peak_loc(imp)`-> (samples,seconds)

* **Imp** - This file should be monophonic full length impulse response.

### Explation
Find the index of the maximum samples level.

### Recomendations
* This function can be used to calibrate the start time of your impulse responses and remove the beging time reversed section
* you should subtract this file by 1 or two samples so you do not cut off the begining of impulse responses.

"""
function peak_loc(imp)

	l=imp.l_samples
	samplerate=imp.samplerate
	colmn=imp.channels
	format=imp.format
	samples=imp.samples
	max_samp=maximum(imp.samples)
	index=1
	while !(samples[index]==max_samp)

		index=index+1
	end

	return (index,Float64(index/samplerate))

end

"""
# seq_create - create multichannel sweep sequence
`seq_create(source,channels::Int64)`-> sweep & inverse sweep in current working directory

* **source** - this is a generated mono sweep loaded by acoustic_load
* **channels** - this is the number of channels to sequence over

### Explation
This function will create a sweep played on each channel independently. This sequence of sweeps will allow easy capturing of channel crossfeed or what is better know as "true stereo" impuleses.

### Recomendations
* make sure the silence following the generated sweep is longer than the expected as that is the only buffer between channels
* Set the level at playback this function will never adjust levels
* When capturing sweep make sure you have atleast -6dB below peak as reverberation can add to the level.

"""
function seq_create(source,channels::Int64)
	samplerate=source.samplerate
	l=source.l_samples
	name=source.name
	source=source.samples
	final=zeros(channels*l,channels)


	for index=1:1:channels
		first_i=(l*index)-l+1
		last_i=index*l
		final[first_i:last_i,index]=source
	end

	wavwrite(final,name*"_"*string(channels)*"_sequence.wav",Fs=samplerate)

end

"""
# Seq_Deconvolve - Generates a set of impulse responses from sine sweeps
`seq_deconvolve(inverse,measured;title::String="",norm::String="u",norm_o=1,lp="h",output="file")`-> N-channel Wav file impulse

* **Inverse** - This file should be monophonic inverse sweep file.
* **Measured** - This should be the N channel captured sweep seq
* **Title** - If you want to name the file something besides the name of the measured sweep with impulse appended
* **Output** - (optional)The "file" argument saves the impulse to a wave file. The "acoustic_load" argument allows you to store the results to a variable. file is the default.
* **norm** - (optional) Controls the final gain applied to the impulse with strings.l = number of samples, n = peak amplitude, u = unnormalized, o = user defined normalized with with norm_0 setting the inverse amplitude
* **norm_o** - (optional) This allow is the custom gain input. This is the level the signal is divided by norm_o so (final impulse)/norm_o. norm_o is in amplitude not decbels.
* **lp** - 1 -full doubled sided impulse & lp>1- cuts from that sample to the end
### Explation
Deconvolve converts a measured sequence sweep into an impulse response using Logarithmic Sine Sweep Method.

### Recomendations
* Audio must be loaded with acoustic_load
* The first argument must be a mono file inverse sweep
* Ardour will make sample accurate edits
* Audacity will not make sample accurate edits try to make measure sweep slighly longer
* Type pwd() to find the where impulses are saved
* Do not use unnormalized as the signal will clip
* Using other normalization you like peak amplitude will destroy the amplitude relationship between impulses. Do not change normalization for a sequence.

See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function seq_deconvolve(inverse,measured;title::String="",norm::String="u",norm_o=1,lp::Int=1,output="file")
	msamplerate=measured.samplerate
	mchan=measured.channels

	isamplerate=inverse.samplerate
	ichan=inverse.channels
	l=inverse.l_samples

	if (msamplerate==isamplerate)&&(mchan>ichan)

		if (output=="acoustic_load")&&(title=="")
			impz=Vector{Acoustic}(undef,mchan)
			for index=1:1:mchan
				first_i=(l*index)-l+1
				last_i=index*l

				temp=Acoustic(measured.samples[first_i:last_i,:],measured.samplerate,(measured.name)*"_chan_"*string(index)*".wav",mchan,measured.format,l)
				impz[index]=deconvolve(inverse,temp,norm=norm,norm_o=norm_o,lp=lp,output="acoustic_load")
			end

			return impz

		elseif (output=="acoustic_load")&&!(title=="")
			impz=Vector{Acoustic}(undef,mchan)
			for index=1:1:mchan
				first_i=(l*index)-l+1
				last_i=index*l

				temp=Acoustic(measured.samples[first_i:last_i,:],measured.samplerate,(measured.name),mchan,measured.format,l)
				impz[index]=deconvolve(inverse,temp,title*"_chan_"*string(index)*".wav",norm=norm,norm_o=norm_o,lp=lp,output="acoustic_load")
			end

			return impz

		end
		#this is if it returns as a file
		if (output=="file")&&(title=="")
			for index=1:1:mchan
				first_i=(l*index)-l+1
				last_i=index*l

				temp=Acoustic(measured.samples[first_i:last_i,:],measured.samplerate,(measured.name)*"_chan_"*string(index)*".wav",mchan,measured.format,l)
				deconvolve(inverse,temp,norm=norm,norm_o=norm_o,lp=lp,output="file")
			end

		elseif (output=="file")&&!(title=="")
			for index=1:1:mchan
				first_i=(l*index)-l+1
				last_i=index*l

				temp=Acoustic(measured.samples[first_i:last_i,:],measured.samplerate,(measured.name),mchan,measured.format,l)
				deconvolve(inverse,temp,title=title*"_chan_"*string(index)*".wav",norm=norm,norm_o=norm_o,lp=lp,output="file")
			end

		end

	else
		error("Samplerates do not match or the sweep does not have more channels than the inverse")

	end

end

function parseval_crop(imp;output="file",precision=7,cut_threshold=120)
	l=imp.l_samples
	samplerate=imp.samplerate
	colmn=imp.channels
	format=imp.format
	samples=imp.samples
	chan=imp.channels
	stop=10.0^(-precision)

	#precalculates the first step & the crop vector for each all channels
	crop=[]
	step_int=1+(l-1)/2
	step_int=floor(step_int)
	step_i=Int(step_int)

	for i in 1:chan
		#Finds the total energy for the channel and precalculate squared amplitude
		samp_thres=threshold.(samples[:,i],cutoff=cut_threshold)
		samp_sqr=abs2.(samp_thres)
		tot_e=sum(samp_sqr)
		tot_e=round(tot_e,digits=precision)

		n=2
		step=step_i
		samp_e=sum(samp_sqr[1:step])
		max_l=l
		#This is designed to take massive steps to find a crop point search range
		while (step>1)&&(tot_e<=round(samp_e,digits=precision))
			step=1+(max_l-1)/factorial(big(n))
			step=floor(big(step))
			step=Int(step)
			samp_e=sum(samp_sqr[1:step])
			n=n+1
		end

		step_lo=step
		step_hi=1+(max_l-1)/factorial(big(n-2))
		step_hi=floor(big(step_hi))
		step_hi=Int(step_hi)

		samp_e=sum(samp_sqr[1:(step_hi-1)])
		#checks if you have already found the value in the first iteration
		if (tot_e-samp_e)>=stop
			crop=vcat(crop,[step_hi])
		else
			#Now Trying to shrink crop search range
			n=1
			max_l=step_hi-1
			step=max_l

			while (step>step_lo)&&(tot_e<=round(samp_e,digits=precision))
				step=1+(max_l-1)/factorial(big(n))
				step=floor(big(step))
				step=Int(step)
				samp_e=sum(samp_sqr[1:step])
				n=n+1
			end

			step_lo=step
			step_hi=1+(max_l-1)/factorial(big(n-2))
			step_hi=floor(big(step_hi))
			step_hi=Int(step_hi)
			step=step_hi


			#brute force search for crop point
			while (step>step_lo)&&!((tot_e-samp_e)>=stop)
				samp_e=sum(samp_sqr[1:step])
				step=step-1
			end
			crop=vcat(crop,[step+1])

		end
	end
	final_crop=maximum(crop)
	cropped=samples[1:final_crop,:]

	if output=="file"
		return wavwrite(cropped,imp.name*".wav",Fs=samplerate)
	elseif output=="acoustic_load"
		return Acoustic(cropped,samplerate,imp.name*".wav",colmn,format,l)
	end

end

"""
scale -"amp" is just linear amplitude and "dB" is decibel level

"""
function gain(source,factor=1.0;scale="amp",output="file")
	l=source.l_samples
	samplerate=source.samplerate
	format=source.format
	samples=source.samples
	chan=source.channels
	if scale=="amp"
		g_samples=factor*samples

	elseif (scale=="db")||(scale=="dB")||(scale=="DB")
		raw=db2amp(factor)
		g_samples=raw*samples

	else
	end



	if output=="file"
		return wavwrite(g_samples,source.name*".wav",Fs=samplerate)
	elseif output=="acoustic_load"
		return Acoustic(g_samples,samplerate,source.name*".wav",chan,format,l)
	end

end

function threshold(x;cutoff=120.0)
	tp=typeof(x)
	thres=abs(cutoff)
	thres=db2amp(-thres)

	if abs(cutoff)>144.49439791871097
		return x
	elseif abs(x)<thres
		return convert(tp,0)
	else
		return x
	end

end

end # module
