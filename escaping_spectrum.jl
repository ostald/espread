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

if get_escaping_df
    for r in runs
        println("Processing run: ", r)
        filter_crit = r

        files = filter(x-> contains(x, filter_crit), dir_con_raw)
        #files = files[1:2]
        
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]
        df_escaping = DataFrame() 

        @time for (id, file) in enumerate(files)
            println("  Processing file: ", file)
            E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df = load_result(joinpath(dir, file))

            df.alt_end = altitude.(df.r)

            for r in eachrow(filter(:alt_end => x-> x > 599e3, df))
                push!(df_escaping, r)
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
end



dir_con = readdir(joinpath(dir, "escaping"))
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)
runs_ = unique([d for d in dir_con_raw])

E_bins = 10 .^collect(0:0.1:log10(8000*1.5))

##
WGLMakie.activate!()
f_esc_E_all = Figure()
sleep(1)
ax_esc = Axis(f_esc_E_all[1, 1], 
    xscale = log10,
    yscale = log10,
    limits = ((10, nothing), (1, 1e6)),
    xlabel = "Energy [eV]", 
    ylabel = "Counts [1]")

#for r_ in runs_
    println("Processing run: ", r_)

    io = open(joinpath(dir, "escaping", r_), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
    df = deserialize(io)
    close(io)
    df.E = E_ev.(norm.(df.v))
    his_E = fit(Histogram, df.E, collect(0:100:E0))
    E_middle = his_E.edges[1][1:end-1] .+ diff(his_E.edges[1])
    #open(joinpath(dir, "escaping", "escaping_" * r * ".bin"), "w") do io
    #    serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals])
    #    serialize(io, df_escaping)
    #end

    #stephist!(ax_esc, df.E, bins = 0:20:E0, label="$(E0)eV_$(lim_pitch_deg)_deg")
    #stephist!(ax_esc, df.E, bins = 10 .^collect(0:0.1:log10(8000*1.5)), label="$(E0)eV_$(lim_pitch_deg)_deg")
    s = stairs!(ax_esc, [his_E.edges[1].+1e-3; E0], [1e-30; his_E.weights; 1e-30], label="$(E0)eV_$(lim_pitch_deg)_deg")

end
ylims!(ax_esc, 1, 1e6)
#ax_esc.yscale[] = log10
ax_esc.xscale[] = log10
##
axislegend(ax_esc)
display(f_esc_E_all)
save(joinpath(dir, "plots", "escaping_hist_allE_pitch.png"), f_esc_E, px_per_unit = 3)





using CairoMakie
CairoMakie.activate!()



f_esc_E = Figure(size = (900, 1200))
sleep(1)
axs_esc_E = [Axis(f_esc_E[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    #limits = (nothing, (1e-4, 1e4)),
    #yscale = log10,
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:5, col in 1:2] 
[ax.xtickformat = values -> ["$(value)" for value in values] for ax in axs_esc_E[end, :]]
[ax.ytickformat = values -> ["$(value)" for value in values] for ax in axs_esc_E[:, 1]]

for (i1, r_) in enumerate(runs_)
    println("Processing run: ", r_)

    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, df_escaping = [0, 0, 0, 0, 0 , 0, 0]
    io = open(joinpath(dir, "escaping", r_), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
    df = deserialize(io)
    close(io)
    df.E = E_ev.(norm.(df.v))
    his_E = fit(Histogram, df.E, collect(0:100:E0))
    E_middle = his_E.edges[1][1:end-1] .+ diff(his_E.edges[1])

    """
    f, ax, h = hist(df.E, bins = 0:20:E0, 
        axis = (xlabel = "Energy [eV]", 
            ylabel = "Counts [1]",
            title = "(E0)eV_(lim_pitch_deg)_deg"
            ),
        )
    save(joinpath(dir, "plots", "escaping_hist_E_$(E0)eV_$(lim_pitch_deg)_deg.png"), f)
    """

    stephist!(axs_esc_E[:][i1], df.E, bins = 0:20:E0, label="$(E0)eV_$(lim_pitch_deg)_deg")
    axislegend(axs_esc_E[:][i1])


end

#display(f_esc_E_all)
#save(joinpath(dir, "plots", "escaping_hist_allE_pitch.png"), f_esc_E, px_per_unit = 3)


#axs_esc_E[1, 1].title = "20 deg"
#axs_esc_E[1, 2].title = "90 deg"

[ax.xlabel = "Energy [eV]" for ax in axs_esc_E[end, :]]
[ax.ylabel = "Counts [1]" for ax in axs_esc_E[:, 1]]

#axislegend(axs_esc_E[1, 1])

linkyaxes!(axs_esc_E[:])
linkxaxes!(axs_esc_E[:])

display(f_esc_E)

save(joinpath(dir, "plots", "panel_escaping_hist_allE_pitch.png"), f_esc_E, px_per_unit = 3)
"""

