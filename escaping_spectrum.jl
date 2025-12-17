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

get_escaping_df = false

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


##

dir_con = readdir(joinpath(dir, "escaping"))
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)
runs_ = unique([d for d in dir_con_raw])
runs_ = [
 "escaping_res_500.0eV_20.0deg_.bin",
 "escaping_res_500.0eV_90.0deg_.bin",
 "escaping_res_1000.0eV_20.0deg_.bin",
 "escaping_res_1000.0eV_90.0deg_.bin",
 "escaping_res_2000.0eV_20.0deg_.bin",
 "escaping_res_2000.0eV_90.0deg_.bin",
 "escaping_res_4000.0eV_20.0deg_.bin",
 "escaping_res_4000.0eV_90.0deg_.bin",
 "escaping_res_8000.0eV_20.0deg_.bin",
 "escaping_res_8000.0eV_90.0deg_.bin"
]


E_edges = 2 .^collect(0:1/4:log2(8000*1.5))./1.024
E_edges = 2 .^ 2 .^collect(0:1/128:log2(log2(8000*1.5))) / 2545.4561526280904 * 2000

##
WGLMakie.activate!()
f_esc_E_all = Figure()
sleep(1)
ax_esc = Axis(f_esc_E_all[1, 1], 
    xscale = log10,
    yscale = log10,
    limits = ((10, nothing), (1e0, 1e5)),
    xlabel = "Energy [eV]", 
    ylabel = "Density [1/log(E)]")

for (i1, r_) in enumerate(runs_)
    println("Processing run: ", r_)

    io = open(joinpath(dir, "escaping", r_), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
    df = deserialize(io)
    close(io)
    df.E = E_ev.(norm.(df.v))
    his_E = fit(Histogram, df.E, E_edges)
    E_middle = his_E.edges[1][1:end-1] .+ diff(his_E.edges[1])
    norm_his = his_E.weights .* diff(log2.(his_E.edges[1]))
    #open(joinpath(dir, "escaping", "escaping_" * r * ".bin"), "w") do io
    #    serialize(io, [E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals])
    #    serialize(io, df_escaping)
    #end

    #stephist!(ax_esc, df.E, bins = 0:20:E0, label="$(E0)eV_$(lim_pitch_deg)_deg")
    #stephist!(ax_esc, df.E, bins = 10 .^collect(0:0.1:log10(8000*1.5)), label="$(E0)eV_$(lim_pitch_deg)_deg")
    color = Makie.wong_colors()[ceil(Int, i1/2)];
    if lim_pitch_deg == 20.0
        s = stairs!(ax_esc, [his_E.edges[1]; E0], [1e-30; norm_his .+ 1e-30; 1e-30], label="$(E0)eV", color = color)
    else
        s = stairs!(ax_esc, [his_E.edges[1]; E0], [1e-30; norm_his .+ 1e-30; 1e-30], linestyle = :dot, color = color)
    end
end

lines!(ax_esc, 1, 1, label = "20 deg", color = "black")
lines!(ax_esc, 1, 1, label = "90 deg", color = "black", linestyle = :dot)
axislegend(ax_esc, position=:lt)
save(joinpath(dir, "plots", "escaping_hist_allE_pitch.png"), f_esc_E_all, px_per_unit = 30)

##

"""
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

