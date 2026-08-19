using LinearAlgebra
using StatsBase
include("analysis_util.jl")
include("constants.jl")


dir = "results/r18_pitchAngle_2026-08-16T02:40:58.062/"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

if !isdir(joinpath(dir, "hist_deltaPitch"))
    mkdir(joinpath(dir, "hist_deltaPitch"))
end

for file in dir_con_raw
    #file = dir_con_raw[2]
    println("Processing file: ", file)
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(joinpath(dir, file))
    n_electrons = size(filter(x -> x.generation == 1, df), 1)
    if n_electrons != 1000
        println("number of electrons in this batch: ", n_electrons)
    end
    #@assert n_electrons == 1000

    df.E0 = E_ev.(norm.(df.v0))
    df.E_end = E_ev.(norm.(df.v))
    sanity_checks(df)

    ## primary electrons:
    df.alt0 = altitude.(df.r0)
    df.alt_end = altitude.(df.r)
    df_alt0 = filter(:alt0 => x -> x > 599e3, df)

    #choose ending altitude for primary electrons, starting altitude for secondaries:
    df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)

    #choose ending position for primary electrons, starting position for secondaries:
    df.pos_earth_centered = ifelse.(df.alt0 .> 599e3, df.r, df.r0)
    #if contains(dir, "conic")
        p0 = [0, 0, c.re]
        df.pos = [p - p0 for p in df.pos_earth_centered]
    #elseif contains(dir, "dipole")
    #    error("not implemented yet")
    #else
    #    error("unknown field model")
    #end


    if contains(dir, "dipole")
        error("not implemented yet. phase, rdist calculation need to take gyrocenter into account.")
    end
    df.phase = phase_angle.(df.pos)
    df.rdist = [norm(p[1:2]) for p in df.pos]

    df.pitch = pitch_angle.(df.v0)


    #define histogram bins
    h_edges = collect(80e3:1e3:600e3)
    r_edges = collect(0:0.1:40)
    p_edges = collect(-pi:pi/36:pi)

    h_middle = h_edges[1:end-1] + diff(h_edges)/2
    r_middle = r_edges[1:end-1] + diff(r_edges)/2
    p_middle = p_edges[1:end-1] + diff(p_edges)/2

    his_hrp = fit(Histogram, (df.alt, df.rdist, df.phase), (h_edges, r_edges, p_edges))
    open(joinpath(dir, "hist_deltaPitch", "h_hrp_$(file[5:end-4]).hist"), "w") do io
        serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp])
    end


    x_edges = collect(-20.5:0.1:20.5)
    y_edges = copy(x_edges)
    z_edges = collect(80e3:1e3:600e3) 
    #intervalls of 1km bc of resolution of atmosphere? => can maybe be updated now. or maybe due to size of histogram?

    x_middle = x_edges[1:end-1] + diff(x_edges)/2
    y_middle = y_edges[1:end-1] + diff(y_edges)/2
    z_middle = z_edges[1:end-1] + diff(z_edges)/2

    his_xyz = fit(Histogram, ([p[1] for p in df.pos], [p[2] for p in df.pos], [p[3] for p in df.pos]), (x_edges, y_edges, z_edges))
    open(joinpath(dir, "hist_deltaPitch", "h_xyz_$(file[5:end-4]).hist"), "w") do io
        serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz])
    end

end


if !isdir(joinpath(dir, "hist_deltaPitch_summed"))
    mkdir(joinpath(dir, "hist_deltaPitch_summed"))
end
dir_con = readdir(joinpath(dir, "hist_deltaPitch"))
runs = unique([d[1:end-10] for d in dir_con])
#radial_runs = filter(x-> contains(x, "hrp"), runs)
#cartesian_runs = filter(x-> contains(x, "xyz"), runs)

for r in runs
    list_h_run = filter(x-> contains(x, r), dir_con)
    println("Summing histograms for run: ", r)
    #sum all histograms
    his_summed = nothing
    ne_pitch_summed = nothing
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = 0, 0, 0, 0, 0, 0
    for file in list_h_run
        println("Processing histogram file: ", file)
        io = open(joinpath(dir, "hist_deltaPitch", file), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_i = deserialize(io)
        close(io)

        if his_summed === nothing
            his_summed = his_i
            #ne_pitch_summed = ne_pitch.weights
        else
            his_summed.edges == his_i.edges || error("Histogram edges do not match!")
            his_summed.weights .+= his_i.weights
            #ne_pitch_summed .+= ne_pitch.weights
        end
    end

    open(joinpath(dir, "hist_deltaPitch_summed", r * "_summed.hist"), "w") do io
        serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_summed])
    end
end
