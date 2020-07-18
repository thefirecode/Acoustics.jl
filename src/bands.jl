module Bands
using DSP
export generateband,ccir,edge
const G=10^(3/10)
const f_r=1000

#class 1 frequency band
const attenuation_1=[-0.4 0.4;-0.4 0.5;-0.4 0.7;-0.4 1.4;16.6 Inf;40.5 Inf;60 Inf;70 Inf];
const attenuation_2=[-0.6 0.6;-0.6 0.7;-0.6 0.9;-0.4 1.7;15.6 Inf;39.5 Inf;54 Inf;60 Inf];
const octave_band=[1,G^(1/8),G^(1/4),G^(3/8),G,G^2,G^3,G^4];
#=
Class 1
Elliptic(3,0.4,70)

Class 2
Elliptic(2,0.6,60)
=#

#=
Bandpass(10^1.5,10^3.9,48000)

=#

#=
Chebyshev1(2,0.2)

 cwgt=digitalfilter(Bandpass(31.5,10^3.9,fs=48000),Butterworth(2))
 cwgt=digitalfilter(Bandpass(20.466,12300,fs=48000),Butterworth(2))

cwgt=digitalfilter(Bandpass(lower(0,1),upper(0,1),fs=48000),Butterworth(2))
=#

#find edge and nyquist checker
#band find
#round(center(13,3)/10^floor(log10(center(13,3))),digits=1)*10^floor(log10(center(13,3)))

#band search
#=
passed class 1 & 2  for octave bands
work 1,3 but, not at 8
Butterworth(3)

=#
b_441k=[0.74183,-0.59185,-0.10701,0.22549,-0.14035,0.094936,-0.24852,0.10223,0.00128,0.001277,-0.039937,-0.13344,0.16029,-0.021732,-0.14884,0.034306,0.028027,0.036722];
a_441k=[1,-2.5531,3.8457,-3.9494,3.2271,-2.1032,1.059,-0.23062,-0.20659,0.3679,-0.33573,0.12494,0.22498,-0.50994,0.48634,-0.31209,0.12765,-0.029342];

ccir_441k=PolynomialRatio(b_441k,a_441k);

b_48k=[0.59277,-0.64278,0.1597,0.010847,-0.0093632,0.045265,-0.14758,0.012525,0.034439,-0.025168,0.021699,-0.043273,-0.069614,0.13181,-0.038812,-0.076451,0.011089,0.019465,0.014415];
a_48k=[1,-3.1159,5.4642,-6.6845,6.5435,-5.4023,3.8997,-2.4889,1.4974,-0.93754,0.70747,-0.62996,0.48873,-0.15423,-0.21824,0.35243,-0.26741,0.11689,-0.026244];

ccir_48k=PolynomialRatio(b_48k,a_48k);

b_96k=[0.058848,-0.46442,1.8276,-4.8221,9.7329,-16.357,24.36,-33.507,43.348,-52.856,60.711,-66.007,68.446,-67.733,63.135,-53.99,40.526,-23.813,4.8912,15.607,-36.79,56.895,-73.89,86.666,-95.231,99.78,-100.04,95.77,-87.469,76.227,-62.832,47.438,-30.193,11.856,6.4667,-24.039,40.664,-56.14,69.828,-80.892,88.819,-93.485,94.739,-92.123,85.144,-73.797,58.772,-41.174,22.155,-2.8324,-15.591,31.86,-44.982,54.487,-60.349,62.69,-61.646,57.519,-50.95,42.841,-34.089,25.434,-17.501,10.857,-5.9149,2.7481,-1.052,0.31721,-0.070364,0.010145,-0.00070572];
a_96k=[1,-11.469,65.039,-243.45,678.85,-1513,2832.2,-4630.6,6827.3,-9296.4,11865,-14292,16293,-17614,18080,-17577,16009,-13316,9548.7,-4904.8,-318.5,5791.3,-11156,16001,-19908,22560,-23823,23722,-22357,19881,-16523,12593,-8424.1,4309.2,-487.57,-2830.2,5465.9,-7313.3,8366.4,-8700.6,8431.4,-7690.4,6625.5,-5392.6,4131.6,-2943.2,1888.6,-1003.8,309.4,190.58,-510.96,680.85,-734.61,705.79,-625.26,520.02,-410.77,310.27,-224.44,155,-101.68,63.036,-36.716,19.897,-9.862,4.3541,-1.6505,0.51104,-0.12018,0.018969,-0.0015016];

ccir_96k=PolynomialRatio(b_96k,a_96k);

function center(x::Int64,b::Int64)
    if iseven(b)==true
        return f_r*(G^((2*x+1)/(2*b)))
    else
        return f_r*(G^(x/b))
    end
end

function lower(x::Int64,b::Int64)
    if iseven(b)==true
        return (f_r*(G^((2*x-2)/(2*b))))
    else
        return f_r*(G^((x-0.5)/b))
    end
end

function upper(x::Int64,b::Int64)
    if iseven(b)==true
        return (f_r*(G^((2*x+2)/(2*b))))
    else
        return f_r*(G^((x+0.5)/b))
    end
end

function edge(b)
i=1
d=-1
x=0
c=center(0,b)
	while c<=20000.0
		i+=1
		c=center(i,b)	
	end

x=center(0,b)

	while 20<=c
		d-=1
		c=center(d,b)
	end

	i=i-1
	d=d+1
	#temp=Int64.((*).(floor.(sign.(LinRange(d,i,(i-d)+1))),ceil.(abs.(LinRange(d,i,(i-d)+1)))))
	temp=(*).(sign.(LinRange(d,i,(i-d)+1)),abs.(LinRange(d,i,(i-d)+1)))
	temp=round.(temp)
	temp=Int64.(temp)
		
	return temp
end

function generateband(bands::Int64,samplerate::Real)
	edgez=edge(bands)
	centerz=center.(edgez,bands)
	centerz_pow=log.(10,centerz)
	centerz_pow=floor.(centerz_pow)
	centerz_pow=(^).(10,centerz_pow)
	centerz_right=(/).(centerz,centerz_pow)
	center_right=round.(centerz_right,digits=2)
	centerz=(*).(center_right,centerz_pow)
	
	return (vcat(Bandpass.(lower.(edgez[begin:(end-1)],bands),upper.(edgez[begin:(end-1)],bands),fs=samplerate),[Highpass(lower(edgez[end],bands),fs=samplerate)]),centerz)

end

function ccir(samplerate::Real)

	if samplerate==44100
		return ccir_441k
	elseif samplerate==48000
		return ccir_48k
	elseif samplerate==96000
		return ccir_96k
	else
		return error("CCIR is not defined at this samplerate")
	end

end

end # end module definition
