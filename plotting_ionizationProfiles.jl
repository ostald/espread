using Serialization
using DataFrames
using LinearAlgebra
using CSV
using ImageFiltering

using MathTeXEngine
set_texfont_family!(FontFamily("TeXGyreHeros"))

using Bonito    
Bonito.set_cleanup_time!(1)
# ssh -L 9384:localhost:9384 user@server

#using WGLMakie
#WGLMakie.activate!()
using CairoMakie
CairoMakie.activate!()

include("analysis_util.jl")
include("constants.jl")
include("magnetic_field.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
#dir = "results/r8_conicB_He_500eV_2026-01-21T19:21:08.258/"
dir = "results/r13_pitchAngle_2026-05-06T17:41:50.299/"
dir = "results/r14_pitchAngle_2026-05-07T14:46:40.781/"
dir = "results/r16_reevaluate-correctedMagField_2026-07-08T17:41:57.774/"

if !isdir(joinpath(dir, "plots"))
    mkdir(joinpath(dir, "plots"))
end

dir_con = readdir(joinpath(dir, "hist_summed"))
dir_con_raw = filter(x-> contains(x, ".hist"), dir_con)
runs = unique(dir_con_raw)
runs = [
    "h_hrp_500.0eV_20.0deg_summed.hist",
    "h_hrp_500.0eV_90.0deg_summed.hist",
    "h_hrp_1000.0eV_20.0deg_summed.hist",
    "h_hrp_1000.0eV_90.0deg_summed.hist",
    "h_hrp_2000.0eV_20.0deg_summed.hist",
    "h_hrp_2000.0eV_90.0deg_summed.hist",
    "h_hrp_4000.0eV_20.0deg_summed.hist",
    "h_hrp_4000.0eV_90.0deg_summed.hist",
    "h_hrp_8000.0eV_20.0deg_summed.hist",
    "h_hrp_8000.0eV_90.0deg_summed.hist",
    "h_hrp_16000.0eV_20.0deg_summed.hist",
    "h_hrp_16000.0eV_90.0deg_summed.hist",
    "h_hrp_32000.0eV_20.0deg_summed.hist",
    "h_hrp_32000.0eV_90.0deg_summed.hist",
    "h_xyz_500.0eV_20.0deg_summed.hist",
    "h_xyz_500.0eV_90.0deg_summed.hist",
    "h_xyz_1000.0eV_20.0deg_summed.hist",
    "h_xyz_1000.0eV_90.0deg_summed.hist",
    "h_xyz_2000.0eV_20.0deg_summed.hist",
    "h_xyz_2000.0eV_90.0deg_summed.hist",
    "h_xyz_4000.0eV_20.0deg_summed.hist",
    "h_xyz_4000.0eV_90.0deg_summed.hist",
    "h_xyz_8000.0eV_20.0deg_summed.hist",
    "h_xyz_8000.0eV_90.0deg_summed.hist",
    "h_xyz_16000.0eV_20.0deg_summed.hist",
    "h_xyz_16000.0eV_90.0deg_summed.hist",
    "h_xyz_32000.0eV_20.0deg_summed.hist",
    "h_xyz_32000.0eV_90.0deg_summed.hist",
]
runs_xyz = filter(x-> contains(x, "xyz"), runs)
runs_xyz_20 = filter(x-> contains(x, "20.0deg"), runs_xyz)
runs_xyz_90 = filter(x-> contains(x, "90.0deg"), runs_xyz)

runs_hrp = filter(x-> contains(x, "hrp"), runs)
runs_hrp_20 = filter(x-> contains(x, "20.0deg"), runs_hrp)
runs_hrp_90 = filter(x-> contains(x, "90.0deg"), runs_hrp)


#____________________________________________________________________________________________________
#ionization profiles
## Plot ionization vs height for different energies and pitch angles
#check atmosphere
df = CSV.read("/nfs/revontuli/data/oliver/espread/results/r4_conicB_2025-09-05T14:19:27.566/atmosphere.txt",
            DataFrame;
            header = false,
            delim = " ", 
            types=Float64,
            ignorerepeated=true)
rename!(df, :Column1 => :height, :Column2 => :nN2, :Column3 => :nO2, :Column4 => :nO)
fig, ax, lin = lines(log10.(df.nN2), df.height, label = "N2 Aurora")
lines!(ax, log10.(df.nO2), df.height, label = "O2 Aurora")
lines!(ax, log10.(df.nO), df.height, label = "O Aurora")
include("get_msis.jl")
hmin = 80e3
hmax = 600e3
hintervals = 1e3
z = hmin+hintervals/2:hintervals:600e3
densityf = make_densityf(hmin+hintervals/2, hmax, hintervals, [69.58, 19.23])
atm = stack(densityf.(z))'
lines!(ax, log10.(atm[:, 1]), z ./1e3, label = "N2 julia_msis", linestyle = :dash)
lines!(ax, log10.(atm[:, 2]), z ./1e3, label = "O2 julia_msis", linestyle = :dash)
lines!(ax, log10.(atm[:, 3]), z ./1e3, label = "O julia_msis", linestyle = :dash)
lines!(ax, log10.(atm[:, 4]), z ./1e3, label = "He julia_msis", linestyle = :dash)
axislegend()
fig
using DelimitedFiles
writedlm( "atmosphere_julia.csv",  atm, ',')
##


#get production models from matlab:
prod_dir = "/nfs/revontuli/data/oliver/espread/results/r4_conicB_2025-09-05T14:19:27.566/plots/"
prod_f = [
"prod500.txt",
"prod1000.txt",
"prod2000.txt",
"prod4000.txt",
"prod8000.txt"]

reference_prod = []
for (i1, f) in enumerate(prod_f)
    df = CSV.read("/nfs/revontuli/data/oliver/espread/results/r4_conicB_2025-09-05T14:19:27.566/plots/"*f,
            DataFrame;
            header = false,
            delim = " ", 
            types=Float64,
            ignorerepeated=true)
    rename!(df, :Column1 => :height, :Column2 => :fang, :Column3 => :serg, :Column4 => :rees)
    #df.height = map(x -> parse(Float64, x), df.height)
    append!(reference_prod, [df])
end

#compare to other production formulas (from matlab elspec)
f = Figure(size = (600, 1000))
sleep(1)
axs = [Axis(f[i, 1], 
        #xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e1, 1e5), (80, 250)),
        xscale = log10,
        xtickformat = values -> ["" for value in values]
        #title = "Production vs Height 20 deg"
        ) for i in 1:5]
linkxaxes!(axs)
for (i,r) in enumerate(runs_xyz_90[1:5])
    #if occursin("500",r) continue end
    #if occursin("1000",r) continue end
    #if occursin("2000",r) continue end
    #if occursin("8000",r) continue end
    
    ax = axs[i]
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    # normalise by density is the same as dividing by bin volume
    # check:
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   

    max_montecarlo = maximum(data)

    color = Makie.wong_colors()[i]

    lines!(ax, data, z_middle/1e3, label = "$E0 eV", color = color)

    df = reference_prod[i]

    lines!(ax, df.fang/maximum(df.fang) * max_montecarlo, df.height, linestyle = :dash, color = color)
    lines!(ax, df.serg/maximum(df.serg) * max_montecarlo, df.height, linestyle = :dot, color = color)
    lines!(ax, df.rees/maximum(df.rees[.!isnan.(df.rees)]) * max_montecarlo, df.height, linestyle = :dashdot, color = color)

    #lines!(ax, df.fang*1e6, df.height, label = "Fang", linestyle = :dash, color = color)
    #lines!(ax, df.serg*1e6, df.height, label = "Sergienko", linestyle = :dot, color = color)
    #lines!(ax, df.rees*1e6, df.height, label = "Rees", linestyle = :dashdot, color = color)
    if i != 1 axislegend(ax) end

end

lines!(axs[1], 0, 0, label = "MC", linestyle = :solid, color = :black)
lines!(axs[1], 0, 0, label = "Fang", linestyle = :dash, color = :black)
lines!(axs[1], 0, 0, label = "Sergienko", linestyle = :dot, color = :black)
lines!(axs[1], 0, 0, label = "Rees", linestyle = :dashdot, color = :black)

axislegend(axs[1], position =:rt)
axs[end].xtickformat = values -> ["$(value)" for value in values]
axs[1].title = "Production vs Height isotropic"
axs[end].xlabel = "Production [m⁻³]"
[ylims!(ax, 80, 200) for ax in axs[4:end]]
#save(joinpath(dir, "plots", "comp_prod_hist_height_xyz_90deg_4000ev.png"), f, px_per_unit = 3.3)
f
##

#comparison of models 4000ev only

f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-4, 1e-2), (100, 180)),
        xscale = log10,
        xticks = LogTicks(-4:-2),
        #title = "Production vs Height isotropic"
        )
for r in runs_xyz_90[4:4]
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    # normalise by density is the same as dividing by bin volume
    # check:
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data/1e6, z_middle/1e3, label = "MC $E0 eV")
    
    max_montecarlo = maximum(data)/1e6

    df = reference_prod[4]

    lines!(ax, df.fang/maximum(df.fang) * max_montecarlo, df.height, label = "Fang", linestyle = :dash)
    lines!(ax, df.serg/maximum(df.serg) * max_montecarlo, df.height, label = "Sergienko", linestyle = :dot)
    lines!(ax, df.rees/maximum(df.rees[.!isnan.(df.rees)]) * max_montecarlo, df.height, label = "Rees", linestyle = :dashdot)
    
    #lines!(ax, df.fang, df.height, label = "Fang", linestyle = :dash)
    #lines!(ax, df.serg, df.height, label = "Sergienko", linestyle = :dot)
    #lines!(ax, df.rees, df.height, label = "Rees", linestyle = :dashdot)


end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_90deg_4000ev.png"), f, px_per_unit = 3.3)
f

##
r = "h_xyz_4000.0eV_90.0deg_summed.hist"
#r = "h_hrp_4000.0eV_63.0deg_summed.hist"
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
close(io)
new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))
his_xyz = rebin(his_xyz, new_edges);
his_xyz = normalize(his_xyz, mode=:density)
data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
df = reference_prod[4]
f, a, l = lines(df.fang*1e6./data, df.height)

f = Figure()
ax = Axis(f[1, 1])
for i1 in 1:5
    r = runs_xyz_90[i1]
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))
    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    df = reference_prod[i1]
    lines!(ax, df.serg*1e6./data, df.height)
    #display(f)
end
xlims!(0, 4)
f

##
#field aligned
f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-1, 1e5), (80, 600)),
        xscale = log10,
        title = "Production vs Height 20 deg"
        )
#for r in runs_xyz_20
for r in filter(x -> contains(x, "4000.0"), runs_xyz)[2:2]
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    #his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    # normalise by density is the same as dividing by bin volume
    # check:
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
axislegend(ax)
#save(joinpath(dir, "plots", "hist_height_xyz_20deg_allE.png"), f, px_per_unit = 3.3)
f


##
# check mirror altitudes!!!
include("constants.jl")
p_angles = 65:1:75
r0 = 600e3 + c.re
r1 = ((sin.(deg2rad.(p_angles)) .^2 * r0^3) .^(1/3) .- c.re) ./1e3


##

f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-1, 1e5), (80, 600)),
        xscale = log10,
        title = "Production vs Height 90 deg"
        )
for r in runs_xyz_90
    println(r)
    #dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    
    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_90deg_allE.png"), f, px_per_unit = 3.3)
f
##

using MathTeXEngine
set_texfont_family!(FontFamily("TeXGyreHeros"))
CairoMakie.activate!()
f = Figure()
sleep(2)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e1, 1e5), (80, 300)),
        xscale = log10,
        title = "Production vs Height"
        )

r = "h_xyz_8000.0eV_90.0deg_summed.hist"
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz_8 = deserialize(io)
close(io)

for (i, r) in enumerate(runs_xyz)
    if mod(ceil(i/2), 2) == 0 continue end
    println(r)
    #dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))
    new_edges = (his_xyz_8.edges[1], his_xyz_8.edges[1], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    
    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    
    color = Makie.wong_colors()[ceil(Int, ceil(Int, i/2)/2)]
    #color = Makie.to_colormap(:batlow)[(ceil(Int, i/2)- 1) * 42 + 1 ]
    if lim_pitch_deg == 20.0
        lines!(ax, data, z_middle/1e3, label = "$E0 eV", color = color)
    else
        lines!(ax, data, z_middle/1e3, linestyle = :dashdot, color = color)
    end
end
#lines!(ax, [0, 0], [0, 0], color = "black", label = "field-aligned θₗᵢₘ = 20 deg")# lim \theta_{lim} = 20^\circ")
lines!(ax, [0, 0], [0, 0], color = "black", label = L"\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}")
lines!(ax, [0, 0], [0, 0], color = "black", linestyle = :dashdot, label = "isotropic")
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_allE.png"), f, px_per_unit = 3.3)
f

ylims!(80, 120)
