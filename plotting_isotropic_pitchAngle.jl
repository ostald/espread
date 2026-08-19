using Serialization
using DataFrames
using LinearAlgebra
#using CSV
#using ImageFiltering

using MathTeXEngine
set_texfont_family!(FontFamily("TeXGyreHeros"))

using Bonito    
Bonito.set_cleanup_time!(1)
# ssh -L 9384:localhost:9384 user@server

using WGLMakie
WGLMakie.activate!()
using CairoMakie
CairoMakie.activate!()

include("analysis_util.jl")

dir = "results/r17_reevaluate-correctedMagField_2026-08-15T20:00:37.960/"
#dir_con = readdir(joinpath(dir, "hist_summed"))
dir_con = readdir(joinpath(dir, "hist_pitch_summed"))
dir_con_raw = filter(x-> contains(x, ".hist"), dir_con)
runs = unique(dir_con_raw)

runs_xyz = filter(x-> contains(x, "xyz"), runs)
runs_hrp = filter(x-> contains(x, "hrp"), runs)

runs_xyz_4kev = filter(x-> contains(x, "4000.0eV"), runs_xyz)
runs_xyz_10kev = filter(x-> contains(x, "10000.0eV"), runs_xyz)
runs_xyz_40kev = filter(x-> contains(x, "40000.0eV"), runs_xyz)

runs_xyz_10deg = filter(x-> contains(x, "10.0deg"), runs_xyz)
runs_xyz_20deg = filter(x-> contains(x, "20.0deg"), runs_xyz)
runs_xyz_30deg = filter(x-> contains(x, "30.0deg"), runs_xyz)
runs_xyz_40deg = filter(x-> contains(x, "40.0deg"), runs_xyz)
runs_xyz_50deg = filter(x-> contains(x, "50.0deg"), runs_xyz)
runs_xyz_60deg = filter(x-> contains(x, "60.0deg"), runs_xyz)
runs_xyz_70deg = filter(x-> contains(x, "70.0deg"), runs_xyz)
runs_xyz_80deg = filter(x-> contains(x, "80.0deg"), runs_xyz)

#runs_xyz_10deg = filter(x-> contains(x, "10.0deg"), runs_xyz)

"""
r = "results/r4_conicB_2025-09-05T14:19:27.566/hist_summed/h_xyz_32000.0eV_90.0deg_summed.hist"
dir = ""
io = open(joinpath(dir, r), "r")
"""
#include("Magne")


if !isdir(joinpath(dir, "plots"))
    mkdir(joinpath(dir, "plots"))
end


##
r = runs[end]
println(r)
io = open(joinpath(dir, "hist_pitch_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_pitch, ne_pitch_summed = deserialize(io)
#E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_pitch = deserialize(io)
close(io)

#his_pitch.weights = (his_pitch.weights' ./ ne_pitch_summed)'
# replace his_pitch with a copy, to be able to replace weights (Int) with
# n_primary normalized weights

h = zero(his_pitch)
h = normalize(h)
h.weights = (his_pitch.weights' ./ ne_pitch_summed)'
h.isdensity = false
his_pitch = h
his_pitch.weights[isnan.(his_pitch.weights)] .= 0

pitch_edges = his_pitch.edges[2]
pitch_middle = pitch_edges[1:end-1] + diff(pitch_edges)/2
z_edges = his_pitch.edges[1]
z_middle = z_edges[1:end-1] + diff(z_edges)/2

#=
fig, ax, hm = heatmap(rad2deg.(pitch_middle), z_middle./1e3, his_pitch.weights',
    axis = (xlabel = "Pitch Angle [deg]", ylabel = "Height [km]"),
    colorscale = log10,
    colorrange = (1e-5, 1e0)
    )
Colorbar(fig[1, 2], hm)
fig
xlims!(0, 90)
=#

id = 68
fig, ax, lin = lines(his_pitch.weights[:, id], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id]))) - $(round(rad2deg(pitch_edges[id+1])))",
    axis = (xscale = log10, limits = ((1e-3, 1e1), nothing)),)
lines!(ax, his_pitch.weights[:, id-1], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id-1]))) - $(round(rad2deg(pitch_edges[id])))")
axislegend(ax)
ylims!(60, 300)
fig

#=
fig = Figure()
ax = Axis(fig[1, 1], xscale = log10, limits = ((1e-4, 1e2), nothing))
for id in 60:70
    lines!(ax, his_pitch.weights[:, id], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id]))) - $(round(rad2deg(pitch_edges[id+1])))",)
end
axislegend(ax)
fig
=#


id = 68
data = dropdims(sum(his_pitch.weights[:, id:end], dims = 2), dims = 2)
fig, ax, lin = lines(his_pitch.weights[:, id], z_middle./1e3, 
    label = "$(round(Int, rad2deg(pitch_edges[id])))° - $(round(Int, rad2deg(pitch_edges[id+1])))°",
    axis = (xscale = log10, 
        limits = ((1e-3, 1e1), (60, 400)),
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        title = "Extended Pitch Angle Source"
        ),
    )
lines!(ax, data, z_middle./1e3, label = "$(round(Int, rad2deg(pitch_edges[id])))° - 90°",)
axislegend(ax)
ylims!(60, 400)

xlims!(1e-3, 1e1)

fig
save(joinpath(dir, "plots", "production_profile_sum90.png"), fig, px_per_unit = 3.3)



