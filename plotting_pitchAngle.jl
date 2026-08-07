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
using CairoMakie
CairoMakie.activate!()
#WGLMakie.activate!()

include("analysis_util.jl")

dir = "results/r9_pitchAngle2026-01-26T11:09:00.654/"
#dir = "results/r8_conicB_He_500eV_2026-01-21T19:21:08.258/"
dir = "results/r12_pitchAngle_2026-03-06T16:58:01.010/"
dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir = "results/r13_pitchAngle_2026-05-06T17:41:50.299/"
dir = "results/r16_reevaluate-correctedMagField_2026-07-08T17:41:57.774/"
dir_con = readdir(joinpath(dir, "hist_summed"))
#dir_con = readdir(joinpath(dir, "hist_pitch_summed"))
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



##
r = runs[2]
println(r)
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_pitch, ne_pitch_summed = deserialize(io)
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_pitch = deserialize(io)
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

fig, ax, hm = heatmap(rad2deg.(pitch_middle), z_middle./1e3, his_pitch.weights',
    axis = (xlabel = "Pitch Angle [deg]", ylabel = "Height [km]"),
    colorscale = log10,
    colorrange = (1e-5, 1e0)
    )
Colorbar(fig[1, 2], hm)
fig
xlims!(0, 90)

id = 67
fig, ax, lin = lines(his_pitch.weights[:, id], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id]))) - $(round(rad2deg(pitch_edges[id+1])))",
    axis = (xscale = log10, limits = ((1e-5, 1e0), nothing)),)
lines!(ax, his_pitch.weights[:, id-1], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id-1]))) - $(round(rad2deg(pitch_edges[id])))")
axislegend(ax)
fig

fig = Figure()
ax = Axis(fig[1, 1], xscale = log10, limits = ((1e-4, 1e2), nothing))
for id in 60:70
    lines!(ax, his_pitch.weights[:, id], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id]))) - $(round(rad2deg(pitch_edges[id+1])))",)
end
axislegend(ax)
fig



id = 68
data = dropdims(sum(his_pitch.weights[:, id:end], dims = 2), dims = 2)
fig, ax, lin = lines(data, z_middle./1e3, label = "Sum $(rad2deg(pitch_edges[id])) deg to 90 deg",
    axis = (xscale = log10, limits = ((1e-5, 1e1), nothing)),)
lines!(ax, his_pitch.weights[:, id], z_middle./1e3, label = "$(round(rad2deg(pitch_edges[id]))) deg - $(round(rad2deg(pitch_edges[id+1]))) deg",)
axislegend(ax)
xlims!(1e-4, 1e0)
fig
save(joinpath(dir, "plots", "production_profile_sum90.png"), fig, px_per_unit = 3.3)





##

dir = "results/r12_pitchAngle_2026-03-06T16:58:01.010/"
dir_con = readdir(joinpath(dir, "hist_summed"))

dir_con_raw = filter(x-> contains(x, ".hist"), dir_con)
runs = unique(dir_con_raw)

runs_xyz = filter(x-> contains(x, "xyz"), runs)
runs_hrp = filter(x-> contains(x, "hrp"), runs)

runs_xyz_4kev = filter(x-> contains(x, "4000.0eV"), runs_xyz)
runs_xyz_10kev = filter(x-> contains(x, "10000.0eV"), runs_xyz)
runs_xyz_40kev = filter(x-> contains(x, "40000.0eV"), runs_xyz)

E0 = 0

f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-6, 1e2), (80, 600)),
        xscale = log10,
        #xticks = LogTicks(-4:-2),
        #title = "Production vs Height isotropic"
        )
#for r in runs_xyz_40kev[1:5]
for r in runs_xyz
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    #dA = diff(his_xyz.edges[1])[1]*diff(his_xyz.edges[2])[1]
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data/1e6, z_middle/1e3, label = "MC $E0 eV, $lim_pitch_deg")
end
axislegend(ax)
xlims!(1e-5, 1e0)
ylims!(80, 400)
f
save(joinpath(dir, "plots", "production_profile_pitch_angle_scan_$(E0)eV.png"), f, px_per_unit = 3.3)
##


draw_plim=string(63)
f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-6, 1e3), (80, 600)),
        xscale = log10,
        #xticks = LogTicks(-4:-2),
        #title = "Production vs Height isotropic"
        )
r =  "h_hrp_10000.0eV_$draw_plim.0deg_summed.hist"
r =  "h_pitch_4000.0eV_$draw_plim.0deg_summed.hist"
println(r)
io = open(joinpath(dir, "hist_summed", r), "r")
io = open(joinpath(dir, "hist_pitch_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)

data = dropdims(sum(his_hrp.weights, dims = (2, 3)), dims = (2, 3))
h_edges = his_hrp.edges[1]
h_middle = h_edges[1:end-1] + diff(h_edges)/2   
lines!(ax, data/1e6, h_middle/1e3, label = "MC $E0 eV, $lim_pitch_deg")

r =  "h_hrp_40000.0eV_$draw_plim.0deg_summed.hist"
println(r)
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)

data = dropdims(sum(his_hrp.weights, dims = (2, 3)), dims = (2, 3))
h_edges = his_hrp.edges[1]
h_middle = h_edges[1:end-1] + diff(h_edges)/2   
lines!(ax, data/1e6, h_middle/1e3, label = "MC $E0 eV, $lim_pitch_deg")
#lines(data/1e6, h_middle/1e3, label = "MC $E0 eV", axis = (xscale = log10,),)
axislegend(ax)
ylims!(80, 400)
xlims!(1e-5, 1e0)
f
save(joinpath(dir, "plots", "production_profile_pitch_angle_$lim_pitch_deg.png"), f, px_per_unit = 3.3)
