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

using WGLMakie
WGLMakie.activate!()
#using CairoMakie
#CairoMakie.activate!()

include("analysis_util.jl")

##
# check mirror altitudes!!!
include("constants.jl")
p_angles = 65:1:75
r0 = 600e3 + c.re
r1 = ((sin.(deg2rad.(p_angles)) .^2 * r0^3) .^(1/3) .- c.re) ./1e3

##

dir = "results/r14_pitchAngle_2026-05-07T14:46:40.781/"
dir_con = readdir(joinpath(dir, "hist_summed"))

dir_con_raw = filter(x-> contains(x, ".hist"), dir_con)
runs = unique(dir_con_raw)

runs_xyz = filter(x-> contains(x, "xyz"), runs)
runs_hrp = filter(x-> contains(x, "hrp"), runs)

runs_xyz_4kev = filter(x-> contains(x, "4000.0eV"), runs_xyz)
runs_xyz_10kev = filter(x-> contains(x, "10000.0eV"), runs_xyz)
runs_xyz_40kev = filter(x-> contains(x, "40000.0eV"), runs_xyz)

E0 = 0

##
# Pitch angle multiplot
# showing production profile for several pitch angles

f = Figure()
sleep(2)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-6, 1e2), (80, 600)),
        xscale = log10,
        #xticks = LogTicks(-4:-2),
        #title = "Production vs Height isotropic"
        )
for r in runs_xyz_10kev[1:4]
#for r in runs_xyz
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    #dA = diff(his_xyz.edges[1])[1]*diff(his_xyz.edges[2])[1]
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data/1e6, z_middle/1e3, label = "$(round(Int, E0/1000)) keV, $(round(Int, lim_pitch_deg)) deg")
end
for r in r1[1:4]
    hlines!(ax, r, label="mirror" * string(r))
end
axislegend(ax)
xlims!(1e-5, 1e0)
xlims!(1e-6, 1e-2)
ylims!(80, 400)
#f
#save(joinpath(dir, "plots", "production_profile_pitch_angle_scan_$(E0)eV.png"), f, px_per_unit = 3.3)
##
#__________________________________________________________________________________________________________________



#energy mulitpplit
#showing the effecto of primary electron energy at the same pitch angle

angle = 66
runs_xyz_deg = filter(x-> contains(x, string(angle)*".0deg"), runs_xyz)

f = Figure()
sleep(2)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-6, 1e2), (80, 600)),
        xscale = log10,
        #xticks = LogTicks(-4:-2),
        title = "Pitch angle "*string(angle)* "deg"
        )
for r in runs_xyz_deg
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    #dA = diff(his_xyz.edges[1])[1]*diff(his_xyz.edges[2])[1]
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data/1e6, z_middle/1e3, label = "$(round(Int, E0/1000)) keV")
end
axislegend(ax)
xlims!(1e-5, 1e0)
xlims!(1e-6, 1e-2)
ylims!(80, 400)
#f
#save(joinpath(dir, "plots", "production_profile_pitch_angle_scan_$(E0)eV.png"), f, px_per_unit = 3.3)
