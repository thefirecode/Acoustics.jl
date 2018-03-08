#sample variable.samplerate
#clarity ,Lq, RT time,strength
module Acoustics
using DSP,Distributions,WAV

export general,c,l,RT,RT_parallel,RT_multi_,d,Ts,sweep,sweep_windowed,deconvolve_complex,deconvolve

function general(source,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

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

function c(source,time,weighting="z",bands="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)
	l=length(source)
	time=Int((time/1000.0)*Int(source.samplerate))


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
f(x)=10*log(10,sum(x[1:time])/sum(x[time:l]))

	if (bands=="b")||(bands=="B")

		source=abs2.(source)

		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]




	results=pmap(x->f(abs2.(filt(digitalfilter(x,Butterworth(2)),source))),bands)

	return hcat(center,results)




	else
		#edge condition
	end

end

function d(source,time,weighting="z",bands="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)
	l=length(source)
	time=Int((time/1000.0)*Int(source.samplerate))


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


#f(x) is the defined function
	f(x) =sum(x[1:time])/sum(x[time:l])

	if (bands=="b")||(bands=="B")

		source=abs2.(source)

		return f(source)


	elseif (bands=="1/3")||(bands=="1/3")

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]


	results=pmap(x->f(abs2.(filt(digitalfilter(x,Butterworth(2)),source))),bands)



	return hcat(center,results)

	else

	end

end

function l(source,percent,weighting="z",band="b" ;s=1)
# percentage in 0 to 100
#source is the audio file
#central limit approximation the lower and
#the levels that won't be exceeded percentage of the time

	samplerate=1.0*Int(source.samplerate)

	l=length(source)
	percent=percent/100.0

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

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

	results=pmap(x->f(abs.(filt(digitalfilter(x,Butterworth(2)),source))),bands)

	return hcat(center,results)

	else

	end

end

function RT_parallel(source,decay,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)
	l_lg=logspace(0,log(10,l),32)


	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end


	function f(x)
		max=sum(x)
		target=[max*10^(-5/20),max*10^(-decay/20)]
		search=round.(logspace(0,log(10,l),100))
		i=1

		while sum(abs2(x[1:search[i]]))>=target[1]
			i_5=search[i]
			i+=1
		end

		while sum(abs2(x[1:i_5]))>=target[1]
				i_5+=1
		end

			i=1
		while sum(abs2(x[1:search[i]]))>=target[2]
			i_target=search[i]
			i+=1
		end

		while sum(abs2(x[1:i_target]))>=target[2]
			i_target+=1
		end
		a=20*log(10,sum(abs2(x[1:i_5])))
		b=20*log(10,sum(abs2(x[1:i_target])))
		m=(a-b)/(i_target-i_5)
		b=a-(m*i_5)

		return (decay-b)/m
	end


	if (band=="b")||(band=="B")

		return f(source)

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]



	results=pmap(x->f(abs2.(filt(digitalfilter(x,Butterworth(2)),source))),bands)



	return hcat(center,results)

	else

	end

end

function RT(source,decay,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="a")||(weighting=="A")

	elseif (weighting=="b")||(weighting=="B")

	elseif (weighting=="c")||(weighting=="C")

	elseif (weighting=="d")||(weighting=="D")

	else
		return print("Weighting undefined")

	end



	function f(x)
		max=sum(abs2.(x))
		target=max*(10^(-decay/(20.0)))
		total=max

		i=1


		if abs2(x[l])> target
			i=Inf

		elseif 	sum(abs2.(x[Int(ceil((1-10^(-5.0))*l)):l]))>target
			i=Int(ceil((1-10^(-5.0))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-4.0))*l)):l]))>target
			i=Int(ceil((1-10^(-4.0))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-3.0))*l)):l]))>target
			i=Int(ceil((1-10^(-3.0))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-2.0))*l)):l]))>target
			i=Int(ceil((1-10^(-2.0))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-1.0))*l)):l]))>target
			i=Int(ceil((1-10^(-1.0))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.5))*l)):l]))>target
			i=Int(ceil((1-10^(-0.5))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.4))*l)):l]))>target
			i=Int(ceil((1-10^(-0.4))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.3))*l)):l]))>target
			i=Int(ceil((1-10^(-0.3))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.2))*l)):l]))>target
			i=Int(ceil((1-10^(-0.2))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.1))*l)):l]))>target
			i=Int(ceil((1-10^(-0.1))*l))

		elseif sum(abs2.(x[Int(ceil((1-10^(-0.05))*l)):l]))>target
			i=Int(ceil((1-10^(-0.05))*l))

		else
			i=1

		end


		while total>=target

			total=sum(abs2.(x[i:l]))

			i+=1

		end


		return i/samplerate
	end

	if (band=="b")||(band=="B")

		return f(source)

	elseif (band=="1/3")||(band=="1/3")

	bands=[Bandpass(11.2/Int(samplerate),14.1/Int(samplerate)),Bandpass(14.1/Int(samplerate),17.8/	Int(samplerate)),Bandpass(17.8/Int(samplerate),22.4/Int(samplerate)),Bandpass(22.4/Int(samplerate),28.2/Int(samplerate)),Bandpass(28.2/Int(samplerate),35.5/Int(samplerate)),Bandpass(35.5/Int(samplerate),44.7/Int(samplerate)),Bandpass(44.7/Int(samplerate),56.2/Int(samplerate)),Bandpass(56.2/Int(samplerate),70.8/Int(samplerate)),Bandpass(70.8/Int(samplerate),89.1/Int(samplerate)),Bandpass(89.1/Int(samplerate),112/Int(samplerate)),Bandpass(112/Int(samplerate),141/Int(samplerate)),Bandpass(141/Int(samplerate),178/Int(samplerate)),Bandpass(178/Int(samplerate),224/Int(samplerate)),Bandpass(224/Int(samplerate),282/Int(samplerate)),Bandpass(282/Int(samplerate),355/Int(samplerate)),Bandpass(355/Int(samplerate),447/Int(samplerate)),Bandpass(447/Int(samplerate),562/Int(samplerate)),Bandpass(562/Int(samplerate),708/Int(samplerate)),Bandpass(708/Int(samplerate),891/Int(samplerate)),Bandpass(891/Int(samplerate),1122/Int(samplerate)),Bandpass(1122/Int(samplerate),1413/Int(samplerate)),Bandpass(1413/Int(samplerate),1778/Int(samplerate)),Bandpass(1778/Int(samplerate),2239/Int(samplerate)),Bandpass(2239/Int(samplerate),2818/Int(samplerate)),Bandpass(2818/Int(samplerate),3548/Int(samplerate)),Bandpass(3548/Int(samplerate),4467/Int(samplerate)),Bandpass(4467/Int(samplerate),5623/Int(samplerate)),Bandpass(5623/Int(samplerate),7079/Int(samplerate)),Bandpass(7079/Int(samplerate),8913/Int(samplerate)),Bandpass(8913/Int(samplerate),11220/Int(samplerate)),Bandpass(11220/Int(samplerate),14130/Int(samplerate)),Bandpass(14130/Int(samplerate),17780/Int(samplerate)),Bandpass(17780/Int(samplerate),0.5)]

	center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]



	results=[filt(digitalfilter(bands[1],Butterworth(2)),source)]


	i=2

	while length(bands)>=i


		results=vcat(results,[filt(digitalfilter(bands[i],Butterworth(2)),source)])

		i+=1

	end



	return hcat(center,f.(results))

	else

	end

end


#stop gap functions to be removed in future release once perforamce problems solved
function RTR_(decay,source,rate)
#this is so an RT60,RT30 over whatever can be taken
	samplerate=1.0*Int(rate)
	max=sum(abs2(source))
	target=max*(10^(-decay/(20)))
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



function Ts(source,weighting="z",band="b" ;s=1)

	samplerate=1.0*Int(source.samplerate)

	l=length(source)

	t=linspace(0,l/samplerate,l)

	if (weighting=="z")||(weighting=="Z")

	elseif (weighting=="a")||(weighting=="A")

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



	results=[filt(digitalfilter(bands[1],Butterworth(2)),source)]


	i=2

	while length(bands)>=i


		results=vcat(results,[filt(digitalfilter(bands[i],Butterworth(2)),source)])

		i+=1

	end



	return hcat(center,f.(results))

	else

	end

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
	wavwrite(sweep,"sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
	wavwrite(isweep,"inverse_sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function sweep_windowed(duration,silence_duration,f_1,f_2,alpha,samplerate)
	#duration in seconds
	sequence=linspace(0,duration,samplerate*duration)
	window=tukey(samplerate*duration,alpha)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(samplerate*silence_duration)
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
	wavwrite(isweep,"inverse_sweep.wav",Fs=samplerate,nbits=32,compression=WAVE_FORMAT_PCM)
end

function deconvolve_complex(sweep,measured)

	l=length(sweep)

	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)
	iimp=imag(ifft(imp))
	iimp=(/).(iimp,l)

	return save("impulse real.wav",rimp[:,1]),save("impulse imag.wav",iimp[:,1])

end

function deconvolve(sweep,measured)

	l=length(sweep)

	sweep=fft(sweep)
	measured=fft(measured)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)


	return save("impulse.wav",rimp[:,1])
end



end
