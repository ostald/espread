#include("espread.jl")
include("setup.jl")

using Distributed
procs = addprocs(10)

@everywhere include("espread.jl")

#f = @spawnat 2  main(4000, 2, alt0, 20, loc_gmag, loc_geod)
#ff = fetch(f)

i = 2
for E0 in e_energy
    for lim_pitch_deg in pitch_limits
        lim_pitch = lim_pitch_deg/180*pi
        f = @spawnat i main(Int(E0), N_electrons, alt0, lim_pitch, loc_gmag, loc_geod)
        i = i+1
    end
end
