using WGLMakie
#using DataFrames
include("analysis_util.jl")
include("constants.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d[1:end-8] for d in dir_con_raw])

#for r in runs[9:end]
    #println("Processing run: ", r)
    #filter_crit = r
    filter_crit = runs[7]

    files = filter(x-> contains(x, filter_crit), dir_con_raw)
    files = files[1:2]

    df_comb = DataFrame()
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

    @time for (id, file) in enumerate(files)
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(file)
        
        #df.E0 = E_ev.(norm.(df.v0))
        #df.E_end = E_ev.(norm.(df.v))
        #sanity_checks(df)

        ## primary electrons:
        df.alt0 = altitude.(df.r0)
        #df.alt_end = altitude.(df.r)
        #df_alt0 = filter(:alt0 => x -> x > 599e3, df)

        #choose ending altitude for primary electrons, starting altitude for secondaries:
        #df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)

        #choose ending position for primary electrons, starting position for secondaries:
        df.pos_earth_centered = ifelse.(df.alt0 .> 599e3, df.r, df.r0)
        if contains(dir, "conic")
            p0 = [0, 0, c.re]
            #df.pos = [p - p0 for p in df.pos_earth_centered]
            for p in df.pos_earth_centered
                push!(df_comb, (pos = p - p0,))
            end
        elseif contains(dir, "dipole")
            error("not implemented yet")
        else
            error("unknown field model")
        end
    end


    df = df_comb

    
    # 3d point plot of endpoints (earth centered coordinates)
    fig = Figure()
    #sleep(2)
    ax = Axis3(fig[1, 1], aspect=(1, 1, 1), viewmode=:fit)
    #sleep(2)
    sc = scatter!(ax,Point3.(df.pos), markersize = 1, alpha = 0.25, transparency = true)#,
        #axis=(limits=(nothing, nothing, nothing),),)
    sleep(2)
    zlims!(ax, 60e3, 6.4e5)
    ylims!(ax, -30, 30)
    xlims!(ax, -30, 30)
    sleep(2)
    """
    fig, ax, sc = scatter(Point3.(df.pos), markersize = 1, 
        axis3d=(limits=((-30, 30), (-30, 30), (60e3, 6.4e5)),
            title = "Electron Position\n$E0 eV, $lim_pitch_deg deg",
            xlabel = "x [m]",
            ylabel = "y [m]",
            zlabel = "z [km]"),
    )
    """
    #display(fig)
    save(joinpath(dir, "plots", "position_$(E0)_$(lim_pitch_deg)_v2.png"), fig)

    nframes = 360
    framerate = 30
    az_iterator = range(1.275*pi, (1.275+2)*pi, length=nframes)
    el_iterator = ones(size(az_iterator))*pi/8
    el_iterator[90:180] .= range(pi/8, 0, length=length(el_iterator[90:180]))
    el_iterator[180:270] .= 0
    el_iterator[270:360] .= range(0, pi/8, length=length(el_iterator[270:360]))

    record(fig, "point_cloud_$(E0)eV_$(lim_pitch_deg)deg_v2.mp4", zip(az_iterator, el_iterator);
            framerate = framerate) do (az, el)
        ax.azimuth = az
        ax.elevation = el
    end

end