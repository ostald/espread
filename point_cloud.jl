using WGLMakie
#using DataFrames
include("analysis_util.jl")
include("constants.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir = "results/r9_pitchAngle2026-01-26T11:09:00.654/"
dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d[1:end-8] for d in dir_con_raw])

#for r in runs[9:end]
    #println("Processing run: ", r)
    #filter_crit = r
    filter_crit = runs[10]

    files = filter(x-> contains(x, filter_crit), dir_con_raw)
    files = files[1:2]

    df = DataFrame()
    df_ion = DataFrame()
    df_primary_injection = DataFrame()
    df_escaping = DataFrame()
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

    #@time for (id, file) in enumerate(files)
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(joinpath(dir, file))
        
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
         
        #if contains(dir, "conic")
            p0 = [0, 0, c.re]
            #df.pos = [p - p0 for p in df.pos_earth_centered]
            for particle in eachrow(filter(:alt_end => x-> x < 599e3, df)) #non-escaping
                push!(df_ion, (pos = (particle.pos_earth_centered - p0) .* [1, 1, 1e3] , generation = particle.generation,))
            end

            for particle in eachrow(filter(:alt_end => x-> x > 599e3, df)) #escaping
                push!(df_escaping, (pos = (particle.r - p0) .* [1, 1, 1e3], generation = particle.generation, v = particle.v,)) 
            end

            for particle in eachrow(filter(:generation => x-> x == 1, df)) #primary
                push!(df_primary_injection, (pos = (particle.r0 - p0) .* [1, 1, 1e3],)) 
            end

        #elseif contains(dir, "dipole")
        #    error("not implemented yet")
        #else
        #    error("unknown field model")
        #end
    end

    fig, ax, h = hist(df_ion.pos)

    alpha = 0.25
    markersize = 1
    
    # 3d point plot of endpoints (earth centered coordinates)
    fig = Figure()
    sleep(2)
    ax = Axis3(fig[1, 1],
        aspect=(1, 1, 1), 
        viewmode=:fit, 
        title = "Electron Position\n$E0 eV, $lim_pitch_deg deg",
        xlabel = "x [m]",
        ylabel = "y [m]",
        zlabel = "z [km]",
        )
    #sleep(2)
    sc = scatter!(ax,Point3.(df_primary_injection.pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "Injected Primary" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:v => x-> E_ev(norm(x)) > E0-0.01, filter(:generation => x-> x == 1, df_escaping)).pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Primary\nno Energy loss" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:v => x-> E_ev(norm(x)) > E0-11.9 && E_ev(norm(x)) < E0-0.01, filter(:generation => x-> x == 1, df_escaping)).pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Primary\nonce scattered, non-ionizing" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:v => x-> E_ev(norm(x)) > E0-40 && E_ev(norm(x)) < E0-20, filter(:generation => x-> x == 1, df_escaping)).pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Primary\nonce ionizing" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:v => x-> E_ev(norm(x)) < E0-40 && E_ev(norm(x)) < E0-80, filter(:generation => x-> x == 1, df_escaping)).pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Primary\ntwice ionizing" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:v => x-> E_ev(norm(x)) > 0, filter(:generation => x-> x == 1, df_escaping)).pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Primary\nall" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:generation => x-> x > 1, df_escaping).pos), 
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping Secondaries" => (; markersize = 15, alpha = 1))#,
    
    sc = scatter!(ax,Point3.(filter(:generation => x-> x > 0, df_ion).pos), 
        markersize = markersize, alpha = alpha, transparency = true, 
        label = "End Primary" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:generation => x-> x == 2, df_ion).pos), 
        markersize = markersize, alpha = alpha, transparency = true, 
        label = "Ionizations Secondaries" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:generation => x-> x == 3, df_ion).pos), 
        markersize = markersize, alpha = alpha, transparency = true, 
        label = "Ionizations Tertiaries" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(filter(:generation => x-> x > 3, df_ion).pos), 
        markersize = markersize, alpha = alpha, transparency = true, 
        label = "Ionizations higher order" => (; markersize = 15, alpha = 1))#,
    #sc = scatter!(ax,Point3.(df.pos), markersize = 1, alpha = 0.25, transparency = true)#,
        #axis=(limits=(nothing, nothing, nothing),),)
    #sleep(2)
    #zlims!(ax, 599.9999, 600.0001)
    #zlims!(ax, 60, 6.4e2)
    ylims!(ax, -30, 30)
    xlims!(ax, -30, 30)
    sleep(2)
    axislegend() 
    save(joinpath(dir, "plots", "position_$(E0)_$(lim_pitch_deg)_pitchAngle.png"), fig)

#end


for r in runs[9:end]
    println("Processing run: ", r)
    filter_crit = r
    #filter_crit = runs[8]

    files = filter(x-> contains(x, filter_crit), dir_con_raw)
    #files = files[1:2]

    df = DataFrame()
    df_ion = DataFrame()
    df_primary_injection = DataFrame()
    df_escaping = DataFrame()
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

    @time for (id, file) in enumerate(files)
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(joinpath(dir, file))
        
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
            for particle in eachrow(filter(:alt_end => x-> x < 599e3, df)) #non-escaping
                push!(df_ion, (pos = (particle.pos_earth_centered - p0) .* [1, 1, 1e-3], generation = particle.generation,))
            end

            for particle in eachrow(filter(:alt_end => x-> x > 599e3, df)) #escaping
                push!(df_escaping, (pos = (particle.r - p0) .* [1, 1, 1e-3], generation = particle.generation, v = particle.v,)) 
            end

            for particle in eachrow(filter(:generation => x-> x == 1, df)) #primary
                push!(df_primary_injection, (pos = (particle.r0 - p0) .* [1, 1, 1e-3],)) 
            end

        elseif contains(dir, "dipole")
            error("not implemented yet")
        else
            error("unknown field model")
        end
    end


    alpha = 0.25
    markersize = 1
    
    # 3d point plot of endpoints (earth centered coordinates)
    fig = Figure()
    sleep(2)
    ax = Axis3(fig[1, 1],
        aspect=(1, 1, 1), 
        viewmode=:fit, 
        title = "Electron Position\n$E0 eV, $lim_pitch_deg deg",
        xlabel = "x [m]",
        ylabel = "y [m]",
        zlabel = "z [km]",
        )
    #sleep(2)
    sc = scatter!(ax,Point3.(df_primary_injection.pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "Injected Primary" => (; markersize = 15, alpha = 1))#,
    sc = scatter!(ax,Point3.(df_escaping.pos),
        markersize = markersize*2, alpha = alpha*2, transparency = true, 
        label = "escaping" => (; markersize = 15, alpha = 1))#,
    
    sc = scatter!(ax,Point3.(df_ion.pos), 
        markersize = markersize, alpha = alpha, transparency = true, 
        label = "Ionizations" => (; markersize = 15, alpha = 1))#,
    
    sleep(2)
    zlims!(ax, 60, 6.4e2)
    ylims!(ax, -30, 30)
    xlims!(ax, -30, 30)
    sleep(2)
    axislegend() 
   
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