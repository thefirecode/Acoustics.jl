#sample variable.samplerate
#clarity ,Lq, RT time,strength
module Acoustics
using DSP,Distributions,WAV,Dierckx

export general,C,L,RT,D,Ts,sweep,sweep_windowed,deconvolve_complex,deconvolve,RT_cal

function general(source,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)

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

function C(source,time,weighting="z",bands="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)
	l=length(source)
	time=Int((time/1000.0)*Int(source.samplerate))
	source=abs2.(source)


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

function D(source,time,weighting="z",bands="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)
	l=length(source)
	time=Int((time/1000.0)*Int(source.samplerate))


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
	f(x) =sum(abs2.(x[1:time]))/sum(abs2.(x[time:l]))

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

function L(source,percent,weighting="z",band="b" ;s=1)
# percentage in 0 to 100
#source is the audio file
#central limit approximation the lower and
#the levels that won't be exceeded percentage of the time

	samplerate=1.0*Int(source.samplerate)

	l=length(source)
	percent=percent/100.0

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end

	f(x)=quantile(Normal(mean(abs.(x)),std(abs.(x))),percent)

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

function RT(source,decay,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)
	sampl_amount=Int(ceil(samplerate/1000))
	sequence=linspace(0,l/samplerate,l)

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


	function f(x)

		x=abs2.(1.0*x[:,1])

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
		time=sampl_amount

	if x[l]>target[2]

		return Inf

	else

#samples the points to use for regression
		while 0<i

			sampled_y=vcat(sampled_y,sum(x[1:i]))
			sampled_x=vcat(sampled_x,sequence[time])

			i-=sampl_amount
			time+=sampl_amount

		end

			sampled_y=vcat(sampled_y,x[1])
			sampled_x=vcat(sampled_x,sequence[l])
			#finishes adding all the points for regression
			model=Spline1D(sampled_x,sampled_y)

#generates the regression for the shroeder plot
			schroeder=abs.(evaluate(model,sequence))


#the -5dB decay point

		total=1
		while (total>=target[1])&&(hi_range<l)
			hi_range+=1
			total=schroeder[hi_range]

		end

		lo_range=hi_range
		while (total>=target[2])&&(lo_range<l)
			lo_range+=1
			total=schroeder[lo_range]
		end

		c,m=linreg(sequence[hi_range:lo_range],10*log.(10,schroeder[hi_range:lo_range]))

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

function Ts(source,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)

	t=linspace(0,l/samplerate,l)

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

		return sum(top)/sum(x)

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
	K=(duration*2*pi*f_1)/log((f_1)/(f_2))
	L=duration/log((f_1)/(f_2))
	m=(2*pi*f_1)/((K/L)*exp((-time)/L))

	return m

end

function sweep(duration,silence_duration,f_1,f_2,samplerate)
	#duration in seconds
	sequence=linspace(0,duration,samplerate*duration)
	sweep=swup.(sequence,duration,f_1,f_2)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2))
	silence=zeros(samplerate*silence_duration)
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
	wavwrite(isweep,"inverse_sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function sweep_windowed(duration,silence_duration,f_1,f_2,alpha,samplerate)
	#duration in seconds
	values=String.([duration,silence_duration,f_1,f_2,alpha,samplerate])
	sequence=linspace(0,duration,samplerate*duration)
	window=tukey(samplerate*duration,alpha)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(samplerate*silence_duration)
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep Duration_"*values[1]*"_Silence_"*values[2]*"_Low Frequency_"*values[3]*"_High Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
	wavwrite(isweep,"Inverse "*"Sweep Duration_"*values[1]*"_Silence_"*values[2]*"_Low Frequency_"*values[3]*"_High Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function deconvolve_complex(sweep,measured,name="")

	l=length(sweep)
	name=String(name)

	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)
	iimp=imag(ifft(imp))
	iimp=(/).(iimp,l)

	return wavwrite(rimp[:,1],name*" impulse real.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM),wavwrite(iimp[:,1],name*" impulse imag.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)

end

function deconvolve(sweep,measured,name="")

	l=length(sweep)
	name=String(name)

	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)


	return wavwrite(rimp[:,1],name*" impulse.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function RT_cal(RT,length,samplerate)
	values=string.([RT,length])
	sequence=linspace(0,length,round(length*samplerate))
	rng=MersenneTwister(1234)
	noise=randn!(rng,zeros(round(length*samplerate)))
	max=maximum([abs(maximum(noise)),abs(minimum(noise))])
	noise=(/).(noise,max)
	sequence=10.^(-(3/RT)*sequence)
	out=(*).(sequence,noise)


	return wavwrite(out,"RT Calibration Length "*values[1]*" Reverberation Time "*values[2]*" .wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)

end

end
