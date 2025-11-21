#escaping particles Energy spectrum analysis

using WGLMakie
using Serialization
using StatsBase
#using DataFrames
include("analysis_util.jl")
include("constants.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d[1:end-8] for d in dir_con_raw])

get_escaping_df = true


E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]
df_escaping = DataFrame() 

if get_escaping_df
    for r in runs
        println("Processing run: ", r)
        filter_crit = r

        files = filter(x-> contains(x, filter_crit), dir_con_raw)
        #files = files[1:2]
        
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

        @time for (id, file) in enumerate(files)
            println("  Processing file: ", file)
            E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(joinpath(dir, file))

            df.alt_end = altitude.(df.r)

            for r in eachrow(filter(:alt_end => x-> x > 599e3, df))
                push!(df_escaping, r)
            end
        end
    end
    if !isdir(joinpath(dir, "escaping"))
        mkdir(joinpath(dir, "escaping"))
    end
    open(joinpath(dir, "escaping", "escaping_" * r * ".bin"), "w") do io
        serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals])
        serialize(io, df_escaping)
    end

end



dir_con = readdir(joinpath(dir, "escaping"))
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)
runs_ = unique([d for d in dir_con_raw])

for r_ in runs_
    println("Processing run: ", r_)

    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df_escaping = [0, 0, 0, 0, 0 , 0, 0]
    io = open(joinpath(dir, "escaping", r_), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
    df = deserialize(io)
    close(io)
    df.E = E_ev.(norm.(df.v))
    his_E = fit(Histogram, df.E, collect(0:100:E0))
    E_middle = his_E.edges[1][1:end-1] .+ diff(his_E.edges[1])

    hist(df.E, bins = 0:20:E0)

end

