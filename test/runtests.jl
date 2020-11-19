using Acoustics
using Base.Test


RT_indef(x)=-((2.0^(-(6.0*x)-1))*(5.0^(-6.0*x)))/(3.0*log(10.0));

single=(samples=(^).(10,((*).(-3,LinRange(0,10,480000)))),samplerate=48000.0,name="test"); #ten second RT60 1s

# write your own tests here
#@test 1 == 2
@test round(RT(single,20),digits=2)==1.0
@test round(C(single,80),digits=2)==round((10*log(10,(RT_indef(0)-RT_indef(0.08))/(RT_indef(0.08)-RT_indef(10)))),digits=2)
@test round(D(single,80),digits=2)==round(((RT_indef(0)-RT_indef(0.08))/(RT_indef(0.08)-RT_indef(10))),digits=2)
#@test round(EDT(single),digits=2)==(
#@test Ts
#@test ST_early
#@test ST_late
#@test IACC #stero
#@test G
#@test J_LF #stereo
#@test L_j

