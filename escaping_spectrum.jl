#escaping Energy

using WGLMakie
using Serialization
#using DataFrames
include("analysis_util.jl")
include("constants.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d[1:end-8] for d in dir_con_raw])

get_escaping_df = true


E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]
df_escaping = nothing 

if get_escaping_df
    for r in runs
        println("Processing run: ", r)
        filter_crit = r

        files = filter(x-> contains(x, filter_crit), dir_con_raw)
        #files = files[1:2]
        
        df_escaping = DataFrame()
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

        @time for (id, file) in enumerate(files)
            E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(file)

            #df.E0 = E_ev.(norm.(df.v0))
            #df.E_end = E_ev.(norm.(df.v))
            #sanity_checks(df)

            ## primary electrons:
            df.alt0 = altitude.(df.r0)
            df.alt_end = altitude.(df.r)
            #df_alt0 = filter(:alt0 => x -> x > 599e3, df)

            #choose ending altitude for primary electrons, starting altitude for secondaries:
            #df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)

            #choose ending position for primary electrons, starting position for secondaries:
            df.pos_earth_centered = ifelse.(df.generation .== 1, df.r, df.r0)

            if contains(dir, "conic")
                p0 = [0, 0, c.re]
                #df.pos = [p - p0 for p in df.pos_earth_centered]

                for particle in eachrow(filter(:alt_end => x-> x > 599e3, df)) #escaping
                    push!(df_escaping, (pos = (particle.r - p0) .* [1, 1, 1e3], generation = particle.generation, v = particle.v,)) 
                end

            elseif contains(dir, "dipole")
                error("not implemented yet")
            else
                error("unknown field model")
            end
        end
    end
    if !isdir(joinpath(dir, "escaping"))
        mkdir(joinpath(dir, "escaping"))
    end
    open(joinpath(dir, "escaping", "escaping_" * r * ".bin"), "w") do io
        serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df_escaping])
    end
end

"""

dir = "results/r4_conicB_2025-09-05T14:19:27.566/escaping"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d for d in dir_con_raw])

#for r in runs
    println("Processing run: ", r)
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(r)


end
"""