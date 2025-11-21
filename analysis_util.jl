using DataFrames
using Serialization
using StatsBase

function load_result(file)
    df = DataFrame(generation = Int[],
            idx_scatter = Int[], 
            r0 = Vector{Float64}[], 
            v0 = Vector{Float64}[], 
            status = Int[], 
            r = Vector{Float64}[], 
            v = Vector{Float64}[]
            )

    io = open(file, "r")
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



function rebin(h, new_edges)
    #1. make sure we have bare counts as weightshis
    @assert(h.isdensity == false)

    #2. create new histogram, assign new edges, weights set to zero 
    h_new = zero(h)
    h_new.edges = new_edges
    h_new.weights = zeros(Float64, map(length, new_edges) .- 1...)    
    
    #3. check that new edges are subset of old edges
    for (ne, oe) in zip(h_new.edges, h.edges)
        for e in ne
            if !any(oe .== e)
                error("New edges must be a subset of old edges")
            end
        end
    end
    #check that new edges are sorted
    for ne in h_new.edges
        issorted(ne) || error("New edges must be sorted")
    end
    #check that new edges have the same endpoints as old edges
    for (d, (ne, oe)) in enumerate(zip(h_new.edges, h.edges))
        if ne[1] != oe[1] || ne[end] != oe[end]
            error("New edges must have the same endpoints as old edges for dimension $d")
        end
    end    

    #4. iterate over new histogram bins, find corresponding old bins, add old weight to new weight
    inds = CartesianIndices(size(h_new.weights))
    for I in inds
        i_lower_edge = I.I
        i_upper_edge = I.I .+ 1
        old_i_lower_edge = map((ne, oe, i) -> findfirst(x -> x == ne[i], oe), h_new.edges, h.edges, i_lower_edge)
        old_i_upper_edge = map((ne, oe, i) -> findfirst(x -> x == ne[i], oe), h_new.edges, h.edges, i_upper_edge)
        #sum over old bins
        #old_inds = CartesianIndices(Tuple(map((ol, ou) -> ol:ou-1, old_i_lower_edge, old_i_upper_edge)))
        old_inds = map((ol, ou) -> ol:ou-1, old_i_lower_edge, old_i_upper_edge)
        h_new.weights[I] = sum(h.weights[old_inds...])
        #for old_I in old_inds
        #    h_new.weights[I] += h.weights[old_I]    
        #end
    end
    return h_new
end


"""
#test rebinning
r = "results/r4_conicB_2025-09-05T14:19:27.566/hist_summed/h_xyz_8000.0eV_90.0deg_summed.hist"
println(r)
io = open(r, "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
close(io)
h = his_xyz
new_edges = (h.edges[1], h.edges[2], filter(x -> mod(x, 1000) == 0, h.edges[3]))
@time his_xyz_rebinned = rebin(his_xyz, new_edges);
"""


function normalize_histogram_density(h::Histogram, volumes)
    @assert(h.isdensity == false)
    # normalise by density is the same as dividing by bin volume
    # check:
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true
    return Histogram(deepcopy(h.edges), h.weights ./ volumes, h.closed, true)
end