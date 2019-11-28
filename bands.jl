const G=10^(3/10)
const f_r=1000

#class 1 frequency band
const attenuation=[-0.4 0.4;-0.4 0.5;-0.4 0.7;-0.4 1.4;16.6 Inf;40.5 Inf;60 Inf;70 Inf];
const octave_band=[1,G^(1/8),G^(1/4),G^(3/8),G,G^2,G^3,G^4];


function fractional_octave(octave,b)

    return 1+(((G^(1/(2*b))-1)/(G^(1/2)-1))*(octave-1))

end

function center(x,b)

    if iseven(b)==true

        return f_r*(G^((2*x+1)/(2*b)))
    else

        return f_r*(G^(x/b))
    end

end

function lower(x,b)
    if iseven(b)==true

        return (f_r*(G^((2*x)/(2*b))))
    else

        return f_r*(G^((x-0.5)/b))
    end

end

function upper(x,b)
    if iseven(b)==true

        return (f_r*(G^((2*x+2)/(2*b))))
    else

        return f_r*(G^((x+0.5)/b))
    end

end
