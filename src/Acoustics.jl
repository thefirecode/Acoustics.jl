module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,FFTW

export general,C,RT,D,Ts,sweep,deconvolve_complex,deconvolve,EDT,acoustic_load,ST_late,ST_early,IACC,G


"""
# Acoustic Load 

acoustic_load(path) - Loads file for processing by other functions in Acoustics.jl

**path** - The location of the supported audio file. Must either be a fullpath or relative to the current working directory. to find the current working directory type pwd(). Look into julia shell for more information


## Example

`julia>` a=acoustic_load("Center Hi Sweep-5 impulse_short.wav") 

The wav file has been loaded into the variable a

`julia>` a.samples

this is an array of individual samples from the loaded files

`julia>` a.samplerate

this is a float contianing the samplerate


`julia>` a.name

this is a string of the file name of imported file




"""
function acoustic_load(path)

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
	samplerate=1.0*Int(temp[2])
	
	else

		
	end

	return (samples=samples,samplerate=samplerate,name=path[p_beg:p_end])


end

function general(source,weighting="z",band="b" ;s=1)

	samplerate=source.samplerate

	l=length(source.samples)

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

"""
# C - Clarity
`C(source,time,weighting,bands)`-> dB

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - "b" (Broadband),"1/3" (third octave bands) [Default b]

### Explation
C is known as Clarity it is the balance between early and late eneregy in an impulse expressed in Decibels (dB). Rooms with a positive C value will have greater percieved definition or clarity in reverberance.The Just Noticeable Diffrence (JND) for clarity metrics is 1 dB.

### Recomendations
* Use 50ms for rooms that will be used for music 
* Use 80ms for rooms that will be used for speech 



See ISO-3382 for more information
"""
function C(source,time,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	time=Int((time/1000.0)*Int(source.samplerate))
	source=source.samples


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time]))/sum(abs2.(x[time:l])))

	if (bands=="b")||(bands=="B")


		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]


	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)


	else
		#edge condition
	end

end
"""
# D - Definition

`D(source,time,weighting,bands)` -> ratio (unitless)

* **Source** - the audio file loaded by acoustic_load
* **Time** - (milliseconds)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - "b" (Broadband),"1/3" (third octave bands) [Default b]


### Explation 
D is known as Definition it is the balance between early and late eneregy in an impulse expressed as ratio. Rooms with a D ratio greater than one will have greater percieved definition or clarity in reverberance.The Just Noticeable Diffrence (JND) for definition metrics is 0,05.

### Recomendations
* Use 50ms for rooms that will be used for music 
* Use 80ms for rooms that will be used for speech 



See ISO-3382 for more information
"""
function D(source,time,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	time=Int((time/1000.0)*Int(source.samplerate))
	source=source.samples


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
	f(x)=sum(abs2.(x[1:time]))/sum(abs2.(x[time:l]))

	if (bands=="b")||(bands=="B")

		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]


	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)



	return hcat(center,results)

	else

	end

end

"""
# RT - Reverberation Time

`RT(source,decay,weighting,bands)` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Decay** - The amount of level decay from steady state (dB)
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - "b" (Broadband),"1/3" (third octave bands) [Default b]


### Explation 
RT is known as Reverberation time it is the measure of decay from steady state to some level. 

### Recomendations
* normally measured in T20 & T30 


See ISO-3382 for more information

"""
function RT(source,decay,weighting="z",band="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	sampl_amount=Int(ceil(samplerate/1000))
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

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

	if (band=="b")||(band=="B")

		return f(source)

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]



	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)

	else

	end

end

"""
 EDT - Early Decay Time

`EDT(source,decay,weighting,bands)` -> time (seconds)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - "b" (Broadband),"1/3" (third octave bands) [Default b]


### Explation 
EDT is known as Early Decay Time it is the measure of decay from peak to 10dB down. 

### Recomendations
* EDT is important as it correlates percieved reverberance

See ISO-3382 for more information

"""
function EDT(source,weighting="z",band="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	sampl_amount=Int(ceil(samplerate/1000))
	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

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

	if (band=="b")||(band=="B")

		return f(source)

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]



	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)

	else

	end

end

"""
# Ts- Centre Time
`Ts(source,weighting,band)`-> milliseconds (ms)

* **Source** - the audio file loaded by acoustic_load
* **Weighting** - the frquency band weightings (Z,A,C,CCIR) [Default Z]
* **Bands** - "b" (Broadband),"1/3" (third octave bands) [Default b]

### Explation
Ts is the time centre which is centre of gravity of the squared impulse repose. It finds to find the center of energy of an impulse response. It is reported in milliseconds. Just Noticeable Diffrence (JND) for the Centre Time (Ts) metrics is 10ms.

### Recomendations
* Avoids the discrete division with C & D 
* Good for finding strongest reflection location



See ISO-3382 for more information
"""
function Ts(source,weighting="z",band="b" ;s=1)
	
	samplerate=source.samplerate

	l=length(source.samples)

	t=LinRange(0,l/samplerate,l)

	source=source.samples

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end

	function f(x)
		x=abs2.(x)
		top=(*).(x,t)

		return (sum(top)/sum(x))*0.001

	end


	if (band=="b")||(band=="B")

		return f(source)

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]


	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]





	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)

	else

	end

end


function ST_early(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	source=source.samples


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end

time_1=Int(ceil(0.01*samplerate))
time_2=Int(ceil(0.02*samplerate))
time_3=Int(ceil(0.1*samplerate))

#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time_1]))/sum(abs2.(x[time_2:time_3])))

	if (bands=="b")||(bands=="B")


		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]


	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)


	else
		#edge condition
	end

end


function ST_late(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
	source=source.samples


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end

time_1=Int(ceil(0.01*samplerate))
time_2=Int(ceil(0.1*samplerate))
time_3=Int(ceil(samplerate))

#f(x) is the defined function
f(x)=10*log(10,sum(abs2.(x[1:time_1]))/sum(abs2.(x[time_2:time_3])))

	if (bands=="b")||(bands=="B")


		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2,14.1;fs=samplerate),Bandpass(14.1,17.8;fs=samplerate),Bandpass(17.8,22.4;fs=samplerate),Bandpass(22.4,28.2;fs=samplerate),Bandpass(28.2,35.5;fs=samplerate),Bandpass(35.5,44.7;fs=samplerate),Bandpass(44.7,56.2;fs=samplerate),Bandpass(56.2,70.8;fs=samplerate),Bandpass(70.8,89.1;fs=samplerate),Bandpass(89.1,112;fs=samplerate),Bandpass(112,141;fs=samplerate),Bandpass(141,178;fs=samplerate),Bandpass(178,224;fs=samplerate),Bandpass(224,282;fs=samplerate),Bandpass(282,355;fs=samplerate),Bandpass(355,447;fs=samplerate),Bandpass(447,562;fs=samplerate),Bandpass(562,708;fs=samplerate),Bandpass(708,891;fs=samplerate),Bandpass(891,1122;fs=samplerate),Bandpass(1122,1413;fs=samplerate),Bandpass(1413,1778;fs=samplerate),Bandpass(1778,2239;fs=samplerate),Bandpass(2239,2818;fs=samplerate),Bandpass(2818,3548;fs=samplerate),Bandpass(3548,4467;fs=samplerate),Bandpass(4467,5623;fs=samplerate),Bandpass(5623,7079;fs=samplerate),Bandpass(7079,8913;fs=samplerate),Bandpass(8913,11220;fs=samplerate),Bandpass(11220,14130;fs=samplerate),Bandpass(14130,17780;fs=samplerate),Bandpass(17780, samplerate*0.5*0.99;fs=samplerate)]


	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(filt(digitalfilter(x,Butterworth(4)),source)),bands)

	return hcat(center,results)


	else
		#edge condition
	end

end


function IACC(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
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


function G(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
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

function J_LF(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
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


function L_j(source,weighting="z",bands="b" ;s=1)

	samplerate=source.samplerate
	l=length(source.samples)
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
	K=(duration*2*pi*f_1)/log((f_1)/(f_2))
	L=duration/log((f_1)/(f_2))

	return sin(K*(exp(-time/L)-1))

end

#inverse generation
function iswup(time,duration,f_1,f_2)
#time seconds
#duraction seconds
#f_1<f_2
	L=duration/log10((f_2)/(f_1))
	m=exp((-time)/L)

	return m

end


"""
# sweep- Logarithmic Sine Sweep
`sweep(duration,silence_duration,f_1,f_2,samplerate,alpha=0.01)`-> sweep & inverse sweep in current working directory

* **Duration** - the duration of the sine sweep in seconds
* **Silence Duration** - the duration of the silence following the sweep
* **F_1** - the start frequency (Hz) of the sweep
* **F_2** - the end frequency (Hz) of the sweep
* **Samplerate** - the samplerate (Hz) of the sine sweep
* **α** -  the mix between a boxcar window and a Hann Window. An α=0 is a boxcar and an α=1 is a Hann window 

### Explation
Creates an expontential sine and inverse sweep in the current working directory. type pwd() to find the current working directory

### Recomendations
* Avoid going all the way up to the nyquist frequency aliasing can occure due the change in frequency  
* make sure your sweep long enough to avoid aliasing



See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
"""
function sweep(duration,silence_duration,f_1,f_2,samplerate,alpha=0.01)
	#duration in seconds
	values=string.([duration,silence_duration,f_1,f_2,alpha,samplerate])
	sequence=LinRange(0.0,duration,Int(floor(samplerate*duration)))
	window=tukey(Int(floor(samplerate*duration)),alpha)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(Int(floor(samplerate*silence_duration)))
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep_Duration_"*values[1]*"_Silence_"*values[2]*"_Low Frequency_"*values[3]*"_High Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse_"*"Sweep_Duration_"*values[1]*"_Silence_"*values[2]*"_Low_Frequency_"*values[3]*"_High_Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate)
end


"""
# deconvolve- generates an impulsem response from an inverse sweep
`deconvolve(inverse,measured,name="")`-> sweepname_impuse.wav

* **Inverse** - the inverse sweep generated by sweep() 
* **Measured** - the captured sweep generated by sweep() it must have the same number of samples per channel as inverse
* **Name** - An optional name

### Explation
deconvolve is turns a sine sweep into an impulse response. A captured sweep must be turned into an impulse to be used inside this library. Using this function will create a WAV file but, will not load the file into memory. inverse should always be the generated file and only needs to be monophonic while the measured sweep can be multichannel. from a signal processing standpoint the channels can be swapped but, the functions automatic name comes from the measured channel.

### Recomendations
* Ardour will make sample accurate edits
* Audacity will not make sample accurate edits
* Type pwd() to find the where impulses are saved



See "Simultaneous measurement of impulse response and distortion with a swept-sine technique" by Angelo Farina for more information
"""
function deconvolve(inverse,measured,name="")

	l=length(measured.samples[:,1])
	name=String(name)
	samplerate=measured.samplerate
	colmn=size(measured.samples)[2]

	if length(name)==0

		name=measured.name
	else
		name=name

	end


	if size(measured.samples)[2]<size(inverse.samples)[2]

	print("fail")


	else size(measured.samples)[2]==size(inverse.samples)[2]
	
	inverse=fft(inverse.samples)
	measured=fft(measured.samples)
	imp=(*).(measured,inverse)
	rimp=ifft(imp)
	#nomalization
	rimp=(/).(rimp,l)


	return wavwrite(rimp,name*"_impulse.wav",Fs=samplerate)

	end

end




end # module
