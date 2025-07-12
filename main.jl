#include("espread.jl")
include("setup.jl")

using Distributed
procs = addprocs(10)

@everywhere include("espread.jl")

#f = @spawnat 2  main(4000, 2, alt0, 20, loc_gmag, loc_geod)
#ff = fetch(f)

"""
f2  = @spawnat 2  main( 500, N_electrons, alt0, 20, loc_gmag, loc_geod)
f3  = @spawnat 3  main( 500, N_electrons, alt0, 90, loc_gmag, loc_geod)
f4  = @spawnat 4  main(1000, N_electrons, alt0, 20, loc_gmag, loc_geod)
f5  = @spawnat 5  main(1000, N_electrons, alt0, 90, loc_gmag, loc_geod)
f6  = @spawnat 6  main(2000, N_electrons, alt0, 20, loc_gmag, loc_geod)
f7  = @spawnat 7  main(2000, N_electrons, alt0, 90, loc_gmag, loc_geod)
f8  = @spawnat 8  main(4000, N_electrons, alt0, 20, loc_gmag, loc_geod)
f9  = @spawnat 9  main(4000, N_electrons, alt0, 90, loc_gmag, loc_geod)
f10 = @spawnat 10 main(8000, N_electrons, alt0, 20, loc_gmag, loc_geod)
f11 = @spawnat 11 main(8000, N_electrons, alt0, 90, loc_gmag, loc_geod)
"""

i_proc = 2
for E0 in e_energy
    for pl in pitch_limits_deg
        @spawnat i_proc  main(E0, N_electrons, alt0, pl, loc_gmag, loc_geod, c)
        i_proc = i_proc+1
        break
    end
end
