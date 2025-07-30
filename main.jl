include("setup.jl")

using Distributed
prcs = addprocs(nprocesses)

@everywhere include("espread.jl")

#global i_proc = 1
for batch in 1:100
    for E0 in e_energy
        @distributed for pitch_lim in pitch_limits_deg
                @time main(E0, N_electrons, alt0, pitch_lim, loc_gmag, loc_geod, c, res_dir; batch=batch)
                print("E0 = ", E0, "\npitch_lim = ", pitch_lim, "\nbatch = ", batch, "\n\n")
            #global i_proc = i_proc+1
        end
    end
end


include("setup.jl")
include("espread.jl")
@profview main(E0, 10, alt0, lim_pitch_deg, loc_gmag, loc_geod, c, res_dir)


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