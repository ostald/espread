using Serialization
using DataFrames
using WGLMakie

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir_con = readdir(joinpath(dir, "hist_summed"))
dir_con_raw = filter(x-> contains(x, ".hist"), dir_con)
runs = unique(dir_con_raw)
runs = [
    "h_hrp_500.0eV_20.0deg_summed.hist",
    "h_hrp_500.0eV_90.0deg_summed.hist",
    "h_hrp_1000.0eV_90.0deg_summed.hist",
    "h_hrp_2000.0eV_20.0deg_summed.hist",
    "h_hrp_2000.0eV_90.0deg_summed.hist",
    "h_hrp_4000.0eV_20.0deg_summed.hist",
    "h_hrp_4000.0eV_90.0deg_summed.hist",
    "h_hrp_8000.0eV_20.0deg_summed.hist",
    "h_hrp_8000.0eV_90.0deg_summed.hist",
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
]
runs_xyz = filter(x-> contains(x, "xyz"), runs)
runs_xyz_20 = filter(x-> contains(x, "20.0deg"), runs_xyz)
runs_xyz_90 = filter(x-> contains(x, "90.0deg"), runs_xyz)

runs_hrp = filter(x-> contains(x, "hrp"), runs)
runs_hrp_20 = filter(x-> contains(x, "20.0deg"), runs_hrp)
runs_hrp_90 = filter(x-> contains(x, "90.0deg"), runs_hrp)

# to do:
# - rebinning
# - normalization from counts to production density: be careful about radial bins!


#for normalization!
normalize(h, mode=:density)


f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Ionizations [1]",
        ylabel = "Height [km]",
        limits = ((1, 1e9), (80, 600)),
        xscale = log10,
        title = "Ionizations vs Height 20 deg"
        )
for r in runs_xyz_20
    println(r)
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    x_edges = his_xyz.edges[1]
    y_edges = his_xyz.edges[2]
    z_edges = his_xyz.edges[3]

    x_middle = x_edges[1:end-1] + diff(x_edges)/2
    y_middle = y_edges[1:end-1] + diff(y_edges)/2
    z_middle = z_edges[1:end-1] + diff(z_edges)/2

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_20deg_allE.png"), f)


f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Ionizations [1]",
        ylabel = "Height [km]",
        limits = ((1, 1e9), (80, 600)),
        xscale = log10,
        title = "Ionizations vs Height 90 deg"
        )
for r in runs_xyz_90
    println(r)
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    x_edges = his_xyz.edges[1]
    y_edges = his_xyz.edges[2]
    z_edges = his_xyz.edges[3]

    x_middle = x_edges[1:end-1] + diff(x_edges)/2
    y_middle = y_edges[1:end-1] + diff(y_edges)/2
    z_middle = z_edges[1:end-1] + diff(z_edges)/2

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_90deg_allE.png"), f)


f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Ionizations [1]",
        ylabel = "Height [km]",
        limits = ((1, 1e9), (80, 600)),
        xscale = log10,
        title = "Ionizations vs Height"
        )
for (i, r) in enumerate(runs_xyz)
    println(r)
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    x_edges = his_xyz.edges[1]
    y_edges = his_xyz.edges[2]
    z_edges = his_xyz.edges[3]

    x_middle = x_edges[1:end-1] + diff(x_edges)/2
    y_middle = y_edges[1:end-1] + diff(y_edges)/2
    z_middle = z_edges[1:end-1] + diff(z_edges)/2

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    color = Makie.wong_colors()[ceil(Int, i/2)]
    if lim_pitch_deg == 20.0
        lines!(ax, data, z_middle/1e3, label = "$E0 eV, 20 deg",  color = color)
    else
        lines!(ax, data, z_middle/1e3, label = "$E0 eV, 90 deg", linestyle = :dash, color = color)
    end
end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_allE.png"), f)




#load histogram
io = open(joinpath(dir, "hist_summed", "h_hrp_1000.0eV_90.0deg_summed.hist"), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
io = open(joinpath(dir, "hist_summed", "h_xyz_1000.0eV_90.0deg_summed.hist"), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
close(io)

#extract edges
h_edges = his_hrp.edges[1]
r_edges = his_hrp.edges[2]
p_edges = his_hrp.edges[3]
h_middle = h_edges[1:end-1] + diff(h_edges)/2
r_middle = r_edges[1:end-1] + diff(r_edges)/2
p_middle = p_edges[1:end-1] + diff(p_edges)/2

# compute gyroradius
include("magnetic_field.jl")
B = norm.([convergent_vertical_field([0, 0, h+c.re]) for h in h_middle])
r_gyro_h_max = v_abs(E0) * c.me ./ (c.qe * B)
lim_pitch = lim_pitch_deg /180 *pi
mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
B_top = convergent_vertical_field([0, 0, 600e3+c.re])
v_perp_mean = mean_vperp_theory * v_abs(E0)
r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)


h_ind = 7


using Bonito    
Bonito.set_cleanup_time!(12)
# ssh -L 9384:localhost:9384 user@server

function plot_histogram_phase_rdist(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
    data = his_hrp.weights[h_ind, :, :]' 
    i_max = [argmax(d) for d in eachrow(data)]

    fig = Figure()
    sleep(5)
    ax = PolarAxis(fig[1:2, 1], title = "Ionizations\n$E0 eV, $lim_pitch_deg deg, $(h_middle[h_ind]/1e3)km",)
    p = voronoiplot!(ax, p_middle, r_middle, data.+1e-3,
        colorscale = log10,
        colorrange = (1, max(1, maximum(data))), lowclip = ("white", 0),
        show_generators = false, strokewidth = 0,)
    rlims!(ax, 0.0, 20.0)
    l1 = lines!(ax, p_middle, r_middle[i_max],
        color= "red", label = "Max Counts")
    l2 = contour!(ax, p_middle, r_middle, data, 
        levels=[maximum(data)/exp(1)], color= "red", 
        linestyle = :dash, label = "1/e contour")
    l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), 
        color= "black", linestyle = :dot, label = "Mean Gyroradius")
    Colorbar(fig[2, 2], p, label="Counts [1]")
    Legend(fig[1, 2], [l1, l2, l4], ["Max Counts", "1/e contour", "Mean Gyroradius"])
    ax.titlegap = 20
    #fig
    save(joinpath(dir, "plots", "hist3d_phase_rd_$(E0)_$(lim_pitch_deg)_$(h_middle[h_ind]/1e3)km.png"), fig)
end

for h_ind in [3, 5, 8, 13, 23]
    plot_histogram_phase_rdist(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
end


function plot_histogram_phase_rdist_normalised(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
    r_edges = his_hrp.edges[2]
    p_edges = his_hrp.edges[3]
    area_pr = [dr * dp / 2  for dp in diff(p_edges), dr in (r_edges[2:end] .^2 - r_edges[1:end-1] .^2)]

    data = his_hrp.weights[h_ind, :, :]' ./ area_pr
    i_max = [argmax(d) for d in eachrow(data)]
    
    fig = Figure()
    sleep(5)
    ax = PolarAxis(fig[1:2, 1], title = "Production\n$E0 eV, $lim_pitch_deg deg, $(h_middle[h_ind]/1e3)km",)
    p = voronoiplot!(ax, p_middle, r_middle, data,
        colorscale = log10,
        colorrange = (10, max(10, maximum(data))), lowclip = ("white", 0),
        show_generators = false, strokewidth = 0,)
    rlims!(ax, 0.0, 20.0)
    l1 = lines!(ax, p_middle, r_middle[i_max],
        color= "red", label = "Max Production")
    l2 = contour!(ax, p_middle, r_middle, data, 
        levels=[maximum(data)/exp(1)], color= "red", 
        linestyle = :dash, label = "1/e contour")
    l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), 
        color= "black", linestyle = :dot, label = "Mean Gyroradius")
    Colorbar(fig[2, 2], p, label="Production [1/m²]")
    Legend(fig[1, 2], [l1, l2, l4], ["Max Production", "1/e contour", "Mean Gyroradius"])
    ax.titlegap = 20
    #fig
    save(joinpath(dir, "plots", "hist3d_phase_rd_norm$(E0)_$(lim_pitch_deg)_$(h_middle[h_ind]/1e3)km.png"), fig)
end

for h_ind in [3, 5, 8, 13, 23]
    plot_histogram_phase_rdist_normalised(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
end


## Radial - Height plot
data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)'
data_h_normal = data ./ maximum(data, dims = 1)
data_h_normal[isnan.(data_h_normal)] .= 0.0
i_max = [argmax(d) for d in eachcol(data)]

fig, ax, hm = heatmap(r_edges, h_edges/1e3, data,
    colorscale = log10,
    colorrange = (1e-2, maximum(data)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 20), (50, 600)),)
    )
sleep(1)
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Counts [1]")
axislegend(ax)
#fig
save(joinpath(dir, "plots", "hist2d_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)



## Radial - Height plot normalised density
r_edges = his_hrp.edges[2]
p_edges = his_hrp.edges[3]
area_pr = [dr * dp / 2  for dp in diff(p_edges), dr in (r_edges[2:end] .^2 - r_edges[1:end-1] .^2)]
## Radial - Height plot normalised density
data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' ./ dropdims(sum(area_pr, dims = 1), dims = 1)
data_h_normal = data ./ maximum(data, dims = 1)
data_h_normal[isnan.(data_h_normal)] .= 0.0
i_max = [argmax(d) for d in eachcol(data)]

fig, ax, hm = heatmap(r_edges, h_edges/1e3, data,
    colorscale = log10,
    colorrange = (1e-2, maximum(data)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 20), (50, 600)),)
    )
sleep(1)
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
contour!(ax, r_middle, h_middle[1:25]/1e3, data_h_normal[:, 1:25], levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Density [1/m²]")
axislegend(ax)
#fig
save(joinpath(dir, "plots", "hist2d_rd_height_norm$(E0)_$(lim_pitch_deg).png"), fig)



x_edges = his_xyz.edges[1]
y_edges = his_xyz.edges[2]
z_edges = his_xyz.edges[3]
x_middle = x_edges[1:end-1] + diff(x_edges)/2
y_middle = y_edges[1:end-1] + diff(y_edges)/2
z_middle = z_edges[1:end-1] + diff(z_edges)/2