module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,Reexport,DataFrames,DataStructures

export L,acoustic_load,filter_verify,acoustic_save,sweep,sweep_target
export WAVE_FORMAT_PCM, WAVE_FORMAT_IEEE_FLOAT, WAVE_FORMAT_ALAW, WAVE_FORMAT_MULAW #exporting WAVE constants
#this contian how to generate third octaves
include("bands.jl");

using Reexport
@reexport using .Bands

#=
Measurement – outputs non-audio values when given audios
Generation – output a with no input
Utility – inputs audio and outputs audio


=#

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

struct Measure
	name::String
	measurement::String
	format::UInt8
	parameter::Dict
	weighting::String
	bands::UInt8
	h_frequency::Float64
	l_frequency::Float64
	samplerate::Float32
	channels::UInt16
	results::DataFrame
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



#weighting function
function frequency_weighting(x::Array{<:AbstractFloat,1},weight::String,Fs::Float32)

	samples=x
	if lowercase(weight)=="z"
		return samples
	elseif lowercase(weight)=="a"
		w_fil=a_wtd(Fs)
		samples=filt(w_fil,samples)
		return samples
	elseif lowercase(weight)=="c"
		w_fil=c_wtd(Fs)
		samples=filt(w_fil,samples)
		return samples
	elseif lowercase(weight)=="ccir"
		w_fil=ccir(Fs)
		samples=filt(w_fil,samples)
		return samples
	else
		error("Unknown Weighting")
	end

end

"""
# Acoustic Load

`acoustic_load(path,format)` -> Loads file(s) in a collection for processing

* **path** - The location of the supported audio file. Must either be a fullpath or relative to the current working directory. to find the current working directory type pwd(). Look into julia shell for more information
* **format** - This used to flag the channel order in a multichannel input for certian functions.

## Example

`julia>` a=acoustic_load("Center Hi Sweep-5 impulse_short.wav")

The wav file(s) has been loaded into the variable a as a linked list. Call collect(a) to convert it to a vector this is not needed most of the time

`julia>` a.samples

this is an array of individual samples from the loaded files

`julia>` a.samplerate

this is a float contianing the samplerate

`julia>` a.name

this is a string of the file name of imported file

`julia>` a.channels

this is a unsigned integer

"""
function acoustic_load end
function acoustic_load(path::String,format::String="")
	if isdir(path)
		loc=MutableLinkedList{String}()
		for file in readdir(path)
			#remove hidden files
			fullp=joinpath(path,file) # See stores the full path of index to reduce recalculation
			if ('.'==file[1])||(isdir(fullp))
				#left blank as we do not want hidden files and directories to be added to the list of files
			else
				fext=splitext(file)[2] #gets the file extension
				fext=fext=lowercase(fext)
				#=
				Adds only wave files to the array location
				=#
				if (".wav"==fext)
					push!(loc,fullp)
				end
			end
		end
		
		#=
		It either says no WAVE files are found or it creates an empty linked list to store the 
		=#
		if length(loc)==0
			error("Directory contians no WAVE files")
		else
			sigz=MutableLinkedList{Acoustic}()
				#do the creation of the acoustic load list
			for file in loc
				name=splitext(basename(file))[1] #get the file name without extension
				temp=wavread(file,format="native")
				samplerate=temp[2] #converts samplerate to proper format
				
				#Deterimine if it will gain anything from being stored as a 64bit float if not will save meory space
				if typeof(temp[1])==Matrix{Int32}
					samples=temp[1]
				else
					samples=Float32.(temp[1])
				end
				#=
				Handles the inputing of acoustic load data
				=#
				sizez=size(samples)
				l_samples=sizez[1]
				chan=convert(UInt16,sizez[2])
				fformat=format_parser(format,chan) #this needs to be fformat (files format) so it can reuse the format the user input
				push!(sigz,Acoustic(samples,samplerate,name,chan,fformat,l_samples))
			end
			
			return sigz
		end
	
	else
		#this is the case of a single file
		file=basename(path) #gets the file name including extension
		if ('.'==file[1])
			error("Is hidden file")
		else
			fext=splitext(file)[2] #gets the file extension
			fext=lowercase(fext) #makes the extension all lowercase
			if (".wav"==fext)
				sigz=MutableLinkedList{Acoustic}() #created an empty linked list
				name=splitext(file)[1] #get the file name without extension
				temp=wavread(path,format="native")
				samplerate=temp[2] #converts samplerate to proper format
				
				#Deterimine if it will gain anything from being stored as a 64bit float if not will save meory space
				if typeof(temp[1])==Matrix{Int32}
					samples=temp[1]
				else
					samples=Float32.(temp[1])
				end
				
				sizez=size(samples)
				l_samples=sizez[1]
				chan=convert(UInt16,sizez[2])
				format=format_parser(format,chan)
				push!(sigz,Acoustic(samples,samplerate,name,chan,format,l_samples))
				return sigz
			else
				error("Not a supported file")
			end
		
		end
		
		
	
	end
	
	

end


"""
# Acoustic Save

`acoustic_save(signal,path="";bits=0,sformat=0)` - exports acoustics collection to audio files

* **signal** - an Acoustic collection
* **path** - the file path to save the files. it defaults to the current workinf directory
* **bits** - specifies the number of bits to be used to encode each sample; the default (0) is an automatic choice based on the values of sform.
* **sformat** - controls the type of encoding used in the file. The options are `WAVE_FORMAT_PCM`, `WAVE_FORMAT_IEEE_FLOAT` ,`WAVE_FORMAT_ALAW` ,` WAVE_FORMAT_MULAW`


## Example

`julia>` a=acoustic_write(a)

### Recomendations
This is merely a pretty warper for wavwrite for the WAV library. You will want to process them if your saving them to something besides float as toyu will get an inexact error.

"""
function acoustic_save end
function acoustic_save(signal,path="";bits=0,sformat=0)

	#if a path is not specified it uses the current working directory
	
	if (typeof(signal)==DataStructures.MutableLinkedList{Acoustics.Acoustic})||(typeof(signal)==Vector{Acoustics.Acoustic})
		
		if path==""
			path=pwd()
		else
			
		end
		
		for file in signal
			
			fullp=joinpath(path,file.name*".wav")
			println("Saving : "*fullp)
			wavwrite(file.samples,fullp,Fs=file.samplerate, nbits=bits, compression=sformat)
		end
		
	else
		error("Unsupported datatype")
	end


end

"""
# L - time-averaged sound level or equivalent continuous sound level
`L(source,channels,tweight,weighting,bands=0,ref=1,h_frequency,l_frequency)`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Channels** - "all" or vector of channels to select
* **Tweight** - selects the averaging time interval (eq, equivalent continuous & time-averaged)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - The number of integer octave band subdivisions aka 1/bands octave band. Where 0 is broadband. [Default 0]
* **ref** - This is the reference level this will scale your recording to the proper level
* **h_frequency** - This is the highest represented frequency
* **l_frequency** - This is the lowest represented frequency

### Explation
This is the average sound level over a specific interval over selected channels and included the overall average of these channels

### Recomendations
* Use this to find the signal to noise ratio
* Use this to find longterm noise level

See ANSI/ASA S1.1-2013 for more information
"""
function Leq(p::Array{<:AbstractFloat},l::Int64,reference_level::AbstractFloat=1)
	sig=abs2.(p)
	sig=sum(sig)
	sig=pow2db(sig)
	scaler=abs2(l)*reference_level
	scaler=pow2db(scaler)
	output=sig-scaler
	output=oftype(p[1],output)
	return output
end

#=
implement weighting
and return measure format
=#
function L end
function L(source::Acoustic,channels::String,tweight::String,weighting::String,bands::Integer=1,ref::AbstractFloat=1.0,h_frequency::Real=20000.0,l_frequency::Real=20.0)
	samps=source.samples
	chan=source.channels
	samp_l=source.l_samples
	samplerate=source.samplerate

	if (h_frequency==20000.0)&&(l_frequency==20.0)
		raw_gen=generateband(bands,samplerate)
	elseif !(h_frequency==20000.0)
		raw_gen=generateband(bands,samplerate,h_frequency)
	else
		raw_gen=generateband(bands,samplerate,h_frequency,l_frequency)
	end

	band_filters=raw_gen[1]
	center_freq=raw_gen[2]
	bands_num=length(center_freq)
	output=Array{Float64, 2}(undef,bands_num,chan)

	if lowercase(weighting)=="z"
		if (lowercase(channels)=="all")
			if (lowercase(tweight)=="eq")||(lowercase(tweight)=="continuous")
				for i in 1:chan
					samp_chan=samps[:,i]
					band_index=1
					for fil in band_filters
						samp_filter=filt(fil,samp_chan)
						output[band_index,i]=Leq(samp_filter,samp_l,ref)
						band_index+=1
					end
				end
			else
			end
		else
			error("Unknown argument for channels")
		end

		return output
	else
		if (lowercase(channels)=="all")
			if (lowercase(tweight)=="eq")||(lowercase(tweight)=="continuous")
				for i in 1:chan
					samp_chan=samps[:,i]
					samp_chan=frequency_weighting(samp_chan,weighting,samplerate)
					band_index=1
					for fil in band_filters
						samp_filter=filt(fil,samp_chan)
						output[band_index,i]=Leq(samp_filter,samp_l,ref)
						band_index+=1
					end
				end
			else
			end
		else
			error("Unknown argument for channels")
		end

		return output
	end

end
function L(source::Acoustic,channels::Vector{Int64},tweight::String,weighting::String,bands::Integer=1,ref::AbstractFloat=1) end
function L(source::Acoustic,channels::String,tweight::String,weighting::String,bands::String,ref::AbstractFloat=1) end
function L(source::Acoustic,channels::Vector{Int64},tweight::String,weighting::String,bands::String,ref::AbstractFloat=1) end

#allows acoustics vector
function L(source::Vector{Acoustic},channels::String,tweight::String,weighting::String,bands::Integer=0,ref::AbstractFloat=1) end
function L(source::Vector{Acoustic},channels::Vector{Int64},tweight::String,weighting::String,bands::Integer=0,ref::AbstractFloat=1) end
function L(source::Vector{Acoustic},channels::String,tweight::String,weighting::String,bands::String,ref::AbstractFloat=1) end
function L(source::Vector{Acoustic},channels::Vector{Int64},tweight::String,weighting::String,bands::String,ref::AbstractFloat=1) end


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
function C end

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
function D end

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
function RT end

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
function EDT end


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
function Ts end

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
function ST_early end

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
function ST_late end

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
function IACC end


function G end

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
function J_LF end


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
function L_j end


function swup(time,K,L)

	return sinpi(K*(exp(time/L)-1))
end

#inverse generation
function iswup(time,L)
	m=exp(-(time)/L)
	return m
end



"""
# sweep - Logarithmic Sine Sweep
`sweep(duration,f_1,f_2,samplerate)`-> sweep & inverse sweep in an acoustic load

* **Duration** - The duration of the sine sweep in seconds
* **F_1** - The start frequency (Hz) of the sweep
* **F_2** - The end frequency (Hz) of the sweep
* **Samplerate** - The samplerate (Hz) of the sine sweep

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. This merely generates the sweep and the inverse it useful for 

### Recomendations
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur due the change in frequency


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep end
function sweep(duration,f_1,f_2,samplerate)
	K=(duration*2*f_1)/log((f_2)/(f_1))
	L=duration/log((f_2)/(f_1))
	sigz=MutableLinkedList{Acoustic}()
	values=string.([duration,f_1,f_2,samplerate])
	n_samp=Int(floor(samplerate*duration))
	sequence=range(0.0,duration,length=n_samp)
	sweep=swup.(sequence,K,L)
	isweep=(*).(sweep[end:-1:1],iswup.(sequence,L))
	sname="Sweep-Duration_"*values[1]*"_f1_"*values[2]*"_f2_"*values[3]*"_Fs_"*values[4]
	isname="Inverse-Sweep-Duration_"*values[1]*"_f1_"*values[2]*"_f2_"*values[3]*"_Fs_"*values[4]
	swpl=Acoustic(sweep,Float32(samplerate),sname,UInt16(1),0x00,n_samp)
	iswpl=Acoustic(isweep,Float32(samplerate),isname,UInt16(1),0x00,n_samp)
	push!(sigz,swpl,iswpl)
	return sigz

end
"""  
# sweep_target - adaptive Logarithmic Sine Sweep


`sweep_target(duration,t60e,predelay,ichan,f_L,f_H,sformat::DataType,samplerate,α=0.0003;inv_norm="p")`-> sweep & inverse sweep

* **Duration** - This is the targeted duration of the sine sweep in seconds
* **t60e** - a best guess of what the maximum reverberation time in seconds
* **predelay** - the maximum travel time in seconds between transmission and reception of the sweep.
* **ichan** - The number of input channels on playback speakers or hardware devices. If you have two playback sources then ichan=2. If capturing a stereo device ichan=2.
* **f_L** - The lowest frequency (Hz) you want to measure in the flat output region 
* **f_H** - The highest frequency you want to measure in the flat output region
* **pbits** - what is the significant bit depth(CD would be 15,24bit audio is 23 and floating point audio would be 24 as that is the size of the mantissa of the floating point number)
* **Samplerate** - The samplerate (Hz) of the sine sweep
* **α** -  The mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window. This parameter controls the tukey window.
* **inv_norm** - This is the normalization gain applied to the inverse sweep. The default value of "p" will normalize to the peak passband amplitude. "a" will normalize to the average passband amplitude. "e" will normalize to the total energy. All other values will cause the signal to be unnormalized.

### Explation
The Logarithmic Sine Sweep is method for generating impulse responses. The longer the sweep more ambient noise suppresion. If alpha is zero a click will be heard. This function differs from sweep in that it tries to find an optimal sweep frequencies edges. Sweep Target also finds the optimal silence duration to stop

### Recomendations
* type pwd() to find the current working directory
* Avoid going all the way up to the nyquist frequency aliasing can occur
* Leave the inverse sweep normalization at either "p" or "a" the unnormalized sweep applies a pretty large gain in the passband and which will cause impulses to clip.
* If measuring the absolute energy is important then use "e" other wise it cause your signal to be extremely quiet for no benefit.


See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
See "SURROUND SOUND IMPULSE RESPONSE Measurement with the Exponential Sine Sweep; Application in Convolution Reverb" by Madeline Carson,Hudson Giesbrecht & Tim Perry for more information (ω_1 needs to be switched with ω_2)
"""
function sweep_target end
function sweep_target(duration,t60e,predelay,ichan,f_L,f_H,pbits::Int,samplerate,α=0.0003;inv_norm="p")
	values=string.([duration,t60e,predelay,ichan,f_L,f_H,pbits,α])
	#calculate the optimal values for the sweep
	edgedom=2*α-2
	N_lh=log(f_L*f_H)*α
	N_L=log(f_L)*2
	N_H=log(f_H)*2
	N_2=(N_lh-N_H)/edgedom
	
	
	if N_2>log(samplerate)
		error("Try a smaller α parameter or reduce high frequency parameter. Current Setting will cause aliasing f_2 : ",exp(N_2))
	end
	
	N_1=(N_lh-N_L)/edgedom
		
	f_1=exp(N_1)
	f_2=exp(N_2)
	f2div1=N_2-N_1 #f2div1=log((f_2)/(f_1))
	#found the optimal
	K=(duration*2*f_1)/f2div1 
	L=duration/f2div1
		
	#calculate optimum silence length and get it to a 10th of a second
	optil=predelay
	optil=optil+log(10,cbrt(2.0^pbits))*t60e
	optil=ceil(optil,digits=1)
	silence=zeros(Int(floor(samplerate*optil)))
	
	#calculate things need for sweep
	sequence=range(0.0,duration,length=Int(floor(samplerate*duration))) #generating the time seteps
	window=tukey(Int(floor(samplerate*duration)),α) #calculates the window
	sweep=(*).(swup.(sequence,K,L),window)
	isweep=(*).(sweep[end:-1:1],iswup.(sequence,L),window)
	sweep=vcat(sweep,silence)
	isweep=vcat(isweep,silence)
	
		
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
	
	#prep the arrays for output
	sigz=MutableLinkedList{Acoustic}()
	values=string.([duration,t60e,predelay,ichan,f_L,f_H,pbits,α])
	values=string.([duration,f_L,f_H,t60e,pbits,α])
	sname="Sweep-time_"*values[1]*"_fl_"*values[2]*"_fH_"*values[3]*"_RTe_"*values[4]*"_bit_"*values[5]*"_α_"*values[6]
	isname="Inverse-Sweep-time_"*values[1]*"_fl_"*values[2]*"_fH_"*values[3]*"_RTe_"*values[4]*"_bit_"*values[5]*"_α_"*values[6]
	
	if ichan>1
		pad=zeros((ichan-1)*sweep_length)
		swp=vcat(sweep,pad)
		temp=swp
		for dex in 1:(ichan-1)
			temp=hcat(temp,shiftsignal(swp,sweep_length*dex))
		end
		swpl=Acoustic(temp,Float32(samplerate),sname,UInt16(ichan),0x04,sweep_length*ichan)
	else
		swpl=Acoustic(sweep,Float32(samplerate),sname,UInt16(1),0x04,sweep_length)
	end
	
	iswpl=Acoustic(isweep,Float32(samplerate),isname,UInt16(1),0x00,sweep_length)

	push!(sigz,swpl,iswpl)	

	return sigz
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
function  deconvolve end

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
function peak_loc end


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
function seq_deconvolve end

function parseval_crop end


"""
scale -"amp" is just linear amplitude and "dB" is decibel level

"""
function gain end


end # module
