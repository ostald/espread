using DataFrames
using Serialization

function load_result(file)
    df = DataFrame(generation = Int[],
            idx_scatter = Int[], 
            r0 = Vector{Float64}[], 
            v0 = Vector{Float64}[], 
            status = Int[], 
            r = Vector{Float64}[], 
            v = Vector{Float64}[]
            )

    io = open(joinpath(dir,  file), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
    data = [deserialize(io) for _ in 1:Int(1e6) if !eof(io)]
    close(io)
    for d in data
        push!(df, d)
    end
    return E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df
end



function sanity_checks(df)
    if size(filter(row -> :status .== 0, df))[1] > 0
        error("Failure of Boris mover!")
    end

    if any(nonunique(df))
        error("Dublicates found.")
    end

    #unfinished propagation:
    unfinished = filter(row -> row.E0 < 12.072 && row.status != -1, df)
    if size(unfinished)[1] > 0
        error("Unfinished propagation found.")
    end
end


function phase_angle(p)
    return atan(p[1], p[2])
end
