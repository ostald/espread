include("setup.jl")

using Distributed
nprocesses = 10 
prcs = addprocs(nprocesses)

@everywhere include("espread.jl")

global i_proc = 1
for E0 in e_energy
    for pitch_lim in pitch_limits_deg
        @time @spawnat workers()[i_proc] main(E0, N_electrons, alt0, pitch_lim, loc_gmag, loc_geod, c)
        global i_proc = i_proc+1
    end
end

#include("setup.jl")
#include("espread.jl")
#main(E0, 10, alt0, pitch_lim, loc_gmag, loc_geod, c)


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