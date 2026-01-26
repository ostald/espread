include("util.jl")
include("setup.jl")
mkdir(res_dir)
save_commit_hash(res_dir)


using Distributed
prcs = addprocs(nprocesses, env=["JULIA_WORKER_TIMEOUT" => "500"])

for p in workers()
    remotecall(include, p, "espread.jl" )
    sleep(10)
end

sleep(60)
#@everywhere include("espread.jl")



#global i_proc = 1
for batch in 2:100
    for E0 in e_energy
        @distributed for pitch_lim_deg in pitch_limits_deg
                @time main(E0, N_electrons, alt0, pitch_lim_deg, loc_gmag, loc_geod, c, res_dir, b_model, nPerGyro; batch=batch)
                print("E0 = ", E0, "\npitch_lim_deg = ", pitch_lim_deg, "\nbatch = ", batch, "\n\n")
            #global i_proc = i_proc+1
        end
    end
end


"""
include("setup.jl")
mkdir(res_dir)
include("espread.jl")
batch = 0
@profview main(E0, 10, alt0, lim_pitch_deg, loc_gmag, loc_geod, c, res_dir, b_model, nPerGyro)
"""

# kill all workers:
#rmprocs(Distributed.workers(), waitfor = 1)

"""
#for i in workers()
#    w = Distributed.worker_from_id(i)
#    kill(w.config.process, Base.SIGKILL)
#end
"""

# check on silent failure:
# Distributed.worker_from_id(2)