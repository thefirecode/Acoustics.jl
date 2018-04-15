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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]



	results=[filt(digitalfilter(bands[1],Butterworth(2)),source)]


	i=2

	while length(bands)>=i


		results=vcat(results,[filt(digitalfilter(bands[i],Butterworth(2)),source)])

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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]




	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)



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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

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
		x=abs2.(x)
		max=sum(x)
		#takes only the first channel
		x=reverse((/).(x[:,1],max))
		target=[(10^(-5/10)),(10^(-((decay+5)/10)))]


#hi and low refer to level
		hi_range=1
		lo_range=1
		sampled_y=[1.0]
		sampled_x=[0.0]
		i=sampl_amount
		total=2

	if x[l]>target[2]

		return Inf

	else

#samples the points to use for regression
		while i < l

			sampled_y=vcat(sampled_y,sum(x[i:l]))
			sampled_x=vcat(sampled_x,sequence[i])

			i+=sampl_amount

		end

			sampled_y=vcat(sampled_y,x[l])
			sampled_x=vcat(sampled_x,sequence[l])
			#finishes adding all the points for regression
			model=Spline1D(sampled_x,sampled_y)

#generates the regression for the shroeder plot
			schroeder=evaluate(model,sequence)


#the -5dB decay point

		total=1
		while (total>=target[1])&&(hi_range<l)
			hi_range+=1
			total=schroeder[hi_range]

		end

		total=1
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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

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

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]



center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

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
