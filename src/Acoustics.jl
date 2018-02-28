#sample variable.samplerate
#clarity ,Lq, RT time,strength
module Acoustics
using DSP,MP3,LibSndFile,Distributions,WAV

export c_,l_,RT_,RT_multi_,c_multi_,sweep,deconvolve,swup,iswup,deconvolve_complex,l_distribution

"time in millisecond
source is the audio file preloaded using MP3 or LibSndFile
"
function c_(time,source)
	time=Int((time/1000.0)*Int(source.samplerate))
	source=abs2(source)

	return 10*log(sum(source[1:time])/sum(source[time:length(source)]))
end

function cr_(time,source,rate)
	time=Int((time/1000.0)*Int(rate))
	source=abs2(source)

	return 10*log(sum(source[1:time])/sum(source[time:length(source)]))
end

function c_multi_(time,source)
	samplerate=1.0*Int(source.samplerate)

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]


center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=[cr_(time,filt(digitalfilter(bands[1],Butterworth(2)),source),samplerate)]

	i=2

	while length(bands)>=i


		results=vcat(results,[cr_(time,filt(digitalfilter(bands[i],Butterworth(2)),source),samplerate)])

		i+=1

	end

	results=hcat(center,results)

	return results
end


"percentage is a percentage from 0-100
source is the audio file preloaded using MP3 or LibSndFile
"
function l_(percent,source)
# percentage in 0 to 100
#source is the audio file
#central limit approximation the lower and
#the levels that won't be exceeded percentage of the time
	source=abs(source)
	percent=percent/100.0
	a=mean(source)
	b=std(source)
	r=Normal(a,b)
	value=quantile(r,percent)
	return value
end

function l_distribution(source)
	sequence=round.(linspace(0,100,100))
	l_dist=l_(sequence,source)

	return l_dist
end

"Decay is the level decay in decibel
source is the audio file preloaded using MP3 or LibSndFile
"
function RT_(decay,source)
#this is so an RT60,RT30 over whatever can be taken
	samplerate=1.0*Int(source.samplerate)
	max=sum(abs2(source))
	target=max*(10^(-decay/(10)))
	total=max


	#Begining processing

	i=1


	if abs2(source[length(source)])> target
		i=Inf

	elseif 	sum(abs2(source[Int(ceil((1-10^(-5.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-5.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-4.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-4.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-3.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-3.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-2.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-2.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-1.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-1.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.5))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.5))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.4))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.4))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.3))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.3))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.2))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.2))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.1))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.1))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.05))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.05))*length(source)))

	else
		i=1

	end


	while total>=target

		total=sum(abs2(source[i:length(source)]))

		i+=1

	end


	return i/samplerate
end

function RTR_(decay,source,rate)
#this is so an RT60,RT30 over whatever can be taken
	samplerate=1.0*Int(rate)
	max=sum(abs2(source))
	target=max*(10^(-decay/(10)))
	total=max


	#Begining processing

	i=1


	if abs2(source[length(source)])> target
		i=Inf

	elseif 	sum(abs2(source[Int(ceil((1-10^(-5.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-5.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-4.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-4.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-3.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-3.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-2.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-2.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-1.0))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-1.0))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.5))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.5))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.4))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.4))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.3))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.3))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.2))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.2))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.1))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.1))*length(source)))

	elseif sum(abs2(source[Int(ceil((1-10^(-0.05))*length(source))):length(source)]))>target
		i=Int(ceil((1-10^(-0.05))*length(source)))

	else
		i=1

	end
	deconvolve

	while total>=target

		total=sum(abs2(source[i:length(source)]))

		i+=1

	end


	return i/samplerate
end

function RT_multi_(decay,source)

	samplerate=1.0*Int(source.samplerate)

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=[RTR_(decay,filt(digitalfilter(bands[1],Butterworth(2)),source),samplerate)]

	i=2

	while length(bands)>=i


		results=vcat(results,[RTR_(decay,filt(digitalfilter(bands[i],Butterworth(2)),source),samplerate)])

		i+=1

	end


	results=hcat(center,results)


	return results


end

function swup(time,duration,f_1,f_2)
#time seconds
#duraction seconds
	K=(duration*2*pi*f_1)/log((f_1)/(f_2))
	L=duration/log((f_1)/(f_2))

	return sin(K*(exp(-time/L)-1))

end

#inverse generation
function iswup(time,duration,f_1,f_2)
#time seconds
#duraction seconds
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

	return wavwrite(sweep,"sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM),wavwrite(isweep,"inverse_sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function deconvolve_complex(sweep,measured)

#	if length(sweep)>length(measured)
#		measured=vcat(measured,zeros(length(sweep)-length(measured)))

#	else
#		sweep=vcat(sweep,zeros(length(measured)-length(sweep)))
#	end


	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,maximum(rimp[:,1]))
	iimp=imag(ifft(imp))
	iimp=(/).(iimp,maximum(rimp[:,1]))

	return save("impulse real.wav",rimp[:,1]),save("impulse imag.wav",iimp[:,1])

end

function deconvolve(sweep,measured)

#	if length(sweep)>length(measured)
#		measured=vcat(measured,zeros(length(sweep)-length(measured)))

#	else
#		sweep=vcat(sweep,zeros(length(measured)-length(sweep)))
#	end

	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,maximum(rimp[:,1]))


	return save("impulse.wav",rimp[:,1])

end


end
