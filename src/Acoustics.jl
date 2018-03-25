#sample variable.samplerate
#clarity ,Lq, RT time,strength
module Acoustics
using DSP,Distributions,WAV

export general,c,l,RT,RT_parallel,RT_multi_,d,Ts,sweep,sweep_windowed,deconvolve_complex,deconvolve,sinc

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
	source=abs2.(source)


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

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

	return hcat(center,results)

	else

	end

end

function RT_parallel(source,decay,weighting="z",band="b" ;s=1)

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

#numerical derivative method taken from "Numerical Evaluation of Derivatives of Functions"(2014) by S. B. Sahoo, M. Acharya and B. P. Acharya
	function f(x)
		samplerate=1.0*Int(x.samplerate)
		x=abs2.(x)
		max=sum(x)
		x=reverse((/).(x,max))
		target=[(10^(-5/20)),(10^(-((decay+5)/20)))]
		z=[complex(0,-1),1,complex(0,1),-1]
		h=[1,complex(0,1),-1,complex(0,-1)]

		taylor=[1,x[1]]
#hi and low refer to level
		hi_range=1
		lo_range=1
		numerator=0
		top_real=0.0
		top_imag=0.0
		denominator=0
		final=0
		intermediate=0
		total=20
		s=0.0000000000002

		#20-point coefficient from DLMF nist
		nodes=[0.076526521133497333755,0.227785851141645078080,0.373706088715419560673,0.510867001950827098004,0.636053680726515025453,0.746331906460150792614,0.839116971822218823395,0.912234428251325905868,0.963971927277913791268,0.993128599185094924786]
		weights=[0.152753387130725850698,0.149172986472603746788,0.142096109318382051329,0.131688638449176626898,0.118194531961518417312,0.101930119817240435037,0.083276741576704748725,0.062672048334109063570,0.040601429800386941331,0.017614007139152118312]

		sequence=linspace(0,l/samplerate,l)

	if x[1]>target[2]

		return Inf

	else

			for m=2:12 #this selects the derivative
				for j=1:4 #this selects the range

					for i=1:length(nodes) #this evaluates the function on one of the ranges

						#intermediate+=h[j]*sum((*).(weights[i],(*).(sinc.((*).((-).(z[j]+h[j]*nodes[i],sequence),samplerate)),(/).(1,(^).((-).(z[j]+h[j]*nodes[i],sequence),m+1)))))
						#intermediate+=h[j]*sum((*).(weights[i],(*).(sinc.((*).((-).(z[j]+h[j]*(-nodes[i]),sequence),samplerate)),(/).(1,(^).((-).(z[j]+h[j]*(-nodes[i]),sequence),m+1)))))
						numerator=(*).((-).(z[j]+s*h[j]*nodes[i],sequence),samplerate)
						top_real=real.(numerator)%(2*pi)
						top_imag=complex(0,imag.(numerator)%(2*pi))
						numerator=((real,imaginary)->sin(real)*cos(imaginary)+cos(real)*sin(imaginary)).(top_real,top_imag)


						print("Numerator Positive:")
						print("\n")
						print(sum(numerator))
						print("\n")

						denominator=(^).((-).(z[j]+s*h[j]*nodes[i],sequence),m+2)
						denominator=(/).(1,denominator)
						intermediate=(*).(numerator,denominator)
						intermediate=(*).(weights[i],intermediate)
						intermediate=sum(intermediate)
						final+=h[j]*intermediate

						#negative nodes
						numerator=(*).((-).(z[j]-s*h[j]*nodes[i],sequence),samplerate)
						top_real=real.(numerator)%(2*pi)
						top_imag=complex(0,imag.(numerator)%(2*pi))
						numerator=((real,imaginary)->sin(real)*cos(imaginary)+cos(real)*sin(imaginary)).(top_real,top_imag)





						denominator=(^).((-).(z[j]-s*h[j]*nodes[i],sequence),m+2)
						denominator=(/).(1,denominator)

						intermediate=(*).(numerator,denominator,x)
						intermediate=(*).(weights[i],intermediate)
						intermediate=sum(intermediate)
						final+=h[j]*intermediate

					end

				end
				final=imag((factorial(m)/complex(0,2*pi))*final)
				taylor=vcat(taylor,final)

				final=0
			end
		end
		#print(taylor)
		#print("\n")

		schroeder=(x->taylor[1]+taylor[2]*x+(taylor[3]/2)*x^2+(taylor[4]/6)*x^3+(taylor[5]/24)*x^4+(taylor[6]/120)*x^5+(taylor[7]/720)*x^6+(taylor[8]/5040)*x^7+(taylor[9]/362880)*x^8+(taylor[10]/362880)*x^9+(taylor[11]/3628800)*x^10+(taylor[12]/39916800)*x^11+(taylor[13]/479001600)*x^12).(sequence)


#the -5dB decay point
		if 	schroeder[Int(ceil((1-10^(-5.0))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-5.0))*l))

		elseif schroeder[Int(ceil((1-10^(-4.0))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-4.0))*l))

		elseif schroeder[Int(ceil((1-10^(-3.0))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-3.0))*l))

		elseif schroeder[Int(ceil((1-10^(-2.0))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-2.0))*l))

		elseif schroeder[Int(ceil((1-10^(-1.0))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-1.0))*l))

		elseif schroeder[Int(ceil((1-10^(-0.5))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.5))*l))

		elseif schroeder[Int(ceil((1-10^(-0.4))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.4))*l))

		elseif schroeder[Int(ceil((1-10^(-0.3))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.3))*l))

		elseif schroeder[Int(ceil((1-10^(-0.2))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.2))*l))

		elseif schroeder[Int(ceil((1-10^(-0.1))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.1))*l))

		elseif schroeder[Int(ceil((1-10^(-0.05))*l))]>target[1]
			hi_range=Int(ceil((1-10^(-0.05))*l))

		else
			hi_range=1
		end


		while total>=target[1]
			total=schroeder[hi_range]
			hi_range+=1
		end


		#decay level
		if schroeder[Int(ceil((1-10^(-0.05))*l))]>target[2]
					lo_range=Int(ceil((1-10^(-0.05))*l))

		elseif schroeder[Int(ceil((1-10^(-0.1))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-0.1))*l))

		elseif schroeder[Int(ceil((1-10^(-0.2))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-0.2))*l))

		elseif schroeder[Int(ceil((1-10^(-0.3))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-0.3))*l))

		elseif schroeder[Int(ceil((1-10^(-0.4))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-0.4))*l))

		elseif schroeder[Int(ceil((1-10^(-0.5))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-0.5))*l))

		elseif schroeder[Int(ceil((1-10^(-1.0))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-1.0))*l))

		elseif schroeder[Int(ceil((1-10^(-2.0))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-2.0))*l))

		elseif schroeder[Int(ceil((1-10^(-3.0))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-3.0))*l))

		elseif schroeder[Int(ceil((1-10^(-4.0))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-4.0))*l))

		elseif 	schroeder[Int(ceil((1-10^(-5.0))*l))]>target[2]
			lo_range=Int(ceil((1-10^(-5.0))*l))

		else
			lo_range=1
		end


		while total>=target[2]
			total=schroeder[lo_range]
			lo_range+=1
		end

		c,m=linreg(schroeder[hi_range:lo_range],sequence[hi_range:lo_range])

		return (10^(-decay/20)-c)/m

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



	return hcat(center,results)


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



center=[12.5,16,20,25,31.5,40,50,63,80,100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,2500,3150,4000,5000,6300,8000,10000,12500,16000,20000]

	results=pmap(x->f(filt(digitalfilter(x,Butterworth(2)),source)),bands)

	return hcat(center,results)

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
