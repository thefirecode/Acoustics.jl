module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,Reexport,DataFrames

export L,acoustic_load,filter_verify

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
function acoustic_load end
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
function sweep end

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
function sweep_target end

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
function seq_create end

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
