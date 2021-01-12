#JULIA_NUM_THREADS=32 julia '/home/phenix/Documents/effectiveness_batch.jl'
using Pkg
Pkg.activate("/home/phenix/Acoustics.jl")
using Acoustics
inv=acoustic_load("/home/phenix/Documents/effectiveness_test/Inverse_Sweep-Duration_30.0_Silence_30.0_Low-Frequency_45.0_High-Frequency_6000.0_Alpha_0.001_Fs_44100.0.wav")
path="/home/phenix/Documents/effectiveness_test/"
#impulse and reverberation
cd("/home/phenix/Documents/effectiveness_test/processed/")
println(pwd())

for (root, dirs, files) in walkdir(path*"Helmholtz_ceiling/")
    println("Files in $root")
    for file in files

            println(joinpath(root, file)) # path to files
            f=joinpath(root, file)
            swp=acoustic_load(f)
	    
            imp=deconvolve(inv,swp,output="acoustic_load",norm="o",norm_o=1.508959209281072e6)
            #Measure Reverberation times
            RT(imp,20,bands=3,output="file")
            RT(imp,30,bands=3,output="file")
            #measure signal level
            #Leq(imp,bands=3,output="file")
            println("File Done!")
    end
end

prinln("Finally Finished")
