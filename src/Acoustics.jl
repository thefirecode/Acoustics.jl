module Acoustics

using DSP,WAV,ReadWriteDlm2,FFTW,Statistics,Distributed,FFTW

export general,C,RT,D,Ts,sweep,deconvolve_complex,deconvolve,EDT,acoustic_load,ST_late,ST_early,IACC,G

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

	while (!(path[i]=='/'))

		p_beg=i

		i-=1
		
	end


	if ("wav"==path[(length(path)-2):end])||("WAV"==path[(length(path)-2):end])

	temp=wavread(path)
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

function sweep(duration,silence_duration,f_1,f_2,samplerate,alpha=0.01)
	#duration in seconds
	values=string.([duration,silence_duration,f_1,f_2,alpha,samplerate])
	sequence=LinRange(0,duration,samplerate*duration)
	window=tukey(samplerate*duration,alpha)
	sweep=(*).(swup.(sequence,duration,f_1,f_2),window)
	isweep=(*).(swup.(sequence[end:-1:1],duration,f_1,f_2),iswup.(sequence,duration,f_1,f_2),window)
	silence=zeros(samplerate*silence_duration)
	sweep=vcat(sweep,silence)
	isweep=vcat(silence,isweep)
	wavwrite(sweep,"Sweep Duration_"*values[1]*"_Silence_"*values[2]*"_Low Frequency_"*values[3]*"_High Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate)
	wavwrite(isweep,"Inverse "*"Sweep Duration_"*values[1]*"_Silence_"*values[2]*"_Low Frequency_"*values[3]*"_High Frequency_"*values[4]*"_Alpha_"*values[5]*"Fs"*values[6]*".wav",Fs=samplerate)
end

function deconvolve_complex(sweep,measured)

	l=length(sweep.samples)
	name=String(name)
	samplerate=measured.samplerate

	if length(name)==0

		name=measured.name
	else
		name=name

	end

	sweep=fft(sweep.samples)
	measured=fft(measured.samples)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)
	iimp=imag(ifft(imp))
	iimp=(/).(iimp,l)

	return wavwrite(rimp[:,1],name*" impulse real.wav",Fs=samplerate),wavwrite(iimp[:,1],name*" impulse imag.wav",Fs=samplerate)

end

function deconvolve(sweep,measured,name="")

	l=length(sweep.samples)
	name=String(name)
	samplerate=measured.samplerate

	if length(name)==0

		name=measured.name
	else
		name=name

	end

	sweep=fft(sweep.samples)
	measured=fft(measured.samples)
	imp=(*).(measured,sweep)
	rimp=real(ifft(imp))
	#nomalization
	rimp=(/).(rimp,l)


	return wavwrite(rimp[:,1],name*" impulse.wav",Fs=samplerate)
end




end # module
