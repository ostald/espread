include("setup.jl")

using Distributed
nprocesses = 10 
prcs = addprocs(nprocesses)

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

#wrkrs = Vector{Any}(undef, nprocesses)

#include("espread.jl")
#main(E0, 10, alt0, pitch_lim, loc_gmag, loc_geod, c)
#println(workers())

global i_proc = 1
for E0 in e_energy
    for pitch_lim in pitch_limits_deg
        @spawnat workers()[i_proc] main(E0, N_electrons, alt0, pitch_lim, loc_gmag, loc_geod, c)
        global i_proc = i_proc+1
    end
end


# kill all workers:
#rmprocs(Distributed.workers(), waitfor = 1)
"""
for i in workers()
    w = Distributed.worker_from_id(i)
    kill(w.config.process, Base.SIGKILL)
end
"""
# check on silent failure:
# Distributed.worker_from_id(2)