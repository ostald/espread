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
using CairoMakie
CairoMakie.activate!()
WGLMakie.activate!()

include("analysis_util.jl")

dir = dir = "results/r9_pitchAngle2026-01-26T11:09:00.654/"
#dir = "results/r8_conicB_He_500eV_2026-01-21T19:21:08.258/"
dir_con = readdir(joinpath(dir, "hist_summed"))
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

include("Magne")



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
for r in runs_xyz_60deg
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    #dA = diff(his_xyz.edges[1])[1]*diff(his_xyz.edges[2])[1]
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    lines!(ax, data/1e6, h_middle/1e3, label = "MC $E0 eV, $lim_pitch_deg")
end
axislegend(ax)
f


f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = (nothing, (80, 600)),
        xscale = log10,
        #xticks = LogTicks(-4:-2),
        #title = "Production vs Height isotropic"
        )
r =  "h_hrp_10000.0eV_60.0deg_summed.hist"
println(r)
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)

data = dropdims(sum(his_hrp.weights, dims = (2, 3)), dims = (2, 3))
h_edges = his_hrp.edges[1]
h_middle = h_edges[1:end-1] + diff(h_edges)/2   
#lines!(ax, data/1e6, h_middle/1e3, label = "MC $E0 eV")
lines(data/1e6, h_middle/1e3, label = "MC $E0 eV", axis = (xscale = log10,),)
axislegend(ax)
f

##

using MathTeXEngine
set_texfont_family!(FontFamily("TeXGyreHeros"))
CairoMakie.activate!()
f = Figure()
sleep(1)
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
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
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

##

using CairoMakie
CairoMakie.activate!()

#polar plot at different heights and energies
f = Figure(size = (1000, 800))
sleep(1)
axs = [PolarAxis(f[i, j],
    thetaticksvisible = false,
    alignmode=Mixed(bottom=0)
    ) for i in 1:4, j in 1:4]
hidedecorations!.(axs, grid=false)

# Add row labels (heights)
heights = [200, 150, 120, 105]
for (i2, h) in enumerate(heights)
    Label(f[i2, 0], "$(h) km", rotation = π/2, tellheight = false)
end

# Add row labels (heights)
#energies = [0.5, 1, 2, 4, 8]
#for (i1, e) in enumerate(energies)
#    Label(f[0, i1], "$(e) keV", tellwidth = false)
#end
Label(f[0, 1], " ", tellwidth = false)


for i in 1:4
    rowgap!(f.layout, i, Relative(-0.0))  # Minimal row gap
end
colgap!(f.layout, 1, Real(0.1))  # Minimal row gap
for j in 2:4
    colgap!(f.layout, j, Relative(-0.02))  # Minimal row gap
end
Label(f[-1, :], "Production Distribution in (r, θ) at Different Heights and Energies",)

p = nothing

to_evaluate = runs_hrp_90[1:2:end]

sleep(1)
for (i2, r) in enumerate(to_evaluate)
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
    close(io)

    Label(f[0, i2], "$E0 keV", tellwidth = false)


    new_edges = (his_hrp.edges[1], his_hrp.edges[2][1:10:end], his_hrp.edges[3])

    his_hrp = rebin(his_hrp, new_edges);

    h_edges = his_hrp.edges[1]
    r_edges = his_hrp.edges[2]
    p_edges = his_hrp.edges[3]
    h_middle = h_edges[1:end-1] + diff(h_edges)/2
    r_middle = r_edges[1:end-1] + diff(r_edges)/2
    p_middle = p_edges[1:end-1] + diff(p_edges)/2

    dv = [dh*dr*dp /2 for dh in diff(h_edges), dr in diff(r_edges.^2), dp in diff(p_edges)]
    his_hrp_norm = normalize_histogram_density(his_hrp, dv)
    println(minimum(his_hrp_norm.weights[his_hrp_norm.weights .!= 0]))

    h_inds = [findfirst(h -> h > height * 1e3, h_middle) for height in heights]

    for (i1, hi) in enumerate(h_inds)
        data = his_hrp_norm.weights[hi, :, :]' 
        i_max = [argmax(d) for d in eachrow(data)]
        println(maximum(data[:]))

        p = voronoiplot!(axs[i1, i2], p_middle, r_middle, data,
            colorscale = log10,
            colorrange = (1e-4, 1e2), lowclip = ("white", 0),
            show_generators = false, strokewidth = 0)
        rlims!(axs[i1, i2], 0.0, 30.0)
        #l1 = lines!(axs[i1, i2], p_middle, r_middle[i_max],
        #    color= "red", label = "Max Counts")
        #l2 = contour!(axs[i1, i2], p_middle, r_middle, data, 
        #    levels=[maximum(data)/exp(1)], color= "red", 
        #    linestyle = :dash, label = "1/e contour")
        #l4 = lines!(ax, p_middle, fill(r_gyro_mean[hi], size(p_middle)), 
        #    color= "black", linestyle = :dot, label = "Mean Gyroradius")
        
        #axislegend(ax)
    end
end
Colorbar(f[1:end, 5], p, label="Production [m⁻³]")
axs[4, 1].rticklabelsvisible = true
text!(axs[4, 1], 5.99, 12, text = "[m]", align = (:center, :top))
save(joinpath(dir, "plots", "panel_radial_hist_$(lim_pitch_deg)_all_E.png"), f, px_per_unit = 3.3)
display(f)


##
#WGLMakie.activate!()

#using CairoMakie
CairoMakie.activate!()


showlines = false
f_hr_ion = Figure(size = (900, 800))
sleep(1)
axs_ion = [Axis(f_hr_ion[row, col], 
    #xlabel = "Radial Distance [m]", 
    #ylabel = "Height [km]", 
    limits = ((nothing, 20), (50, 600)),
    ytickformat = values -> ["" for value in values], 
    xtickformat = values -> ["" for value in values]
    ) for row in 1:3, col in 1:2] 
[ax.xtickformat = values -> ["$(value)" for value in values] for ax in axs_ion[end, :]]
[ax.ytickformat = values -> ["$(value)" for value in values] for ax in axs_ion[:, 1]]
[ax.xminorticks = [1, 2, 3, 4] for ax in axs_ion[:]]
[ax.xminorticksvisible = true for ax in axs_ion[:]]
[ax.xminorgridvisible = true for ax in axs_ion[:]]


f_hr_prod = Figure(size = (900, 800))
sleep(1)
axs_prod = [Axis(f_hr_prod[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    limits = ((nothing, 20), (50, 600)),
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:3, col in 1:2] 
[ax.xtickformat = values -> ["$(Int(value))" for value in values] for ax in axs_prod[end, :]]
[ax.ytickformat = values -> ["$(Int(value))" for value in values] for ax in axs_prod[:, 1]]
[ax.xminorticks = [1, 2, 3, 4] for ax in axs_prod[:]]
[ax.xminorticksvisible = true for ax in axs_prod[:]]
[ax.xminorgridvisible = true for ax in axs_prod[:]]


f_prod_r_h = Figure(size = (900, 800))
sleep(1)
axs_r_h = [Axis(f_prod_r_h[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    limits = ((-20, 20), (2e-3, 2e3)),
    yscale = log10,
    xticks = [-30, -20, -10, 0, 10, 20, 30],
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:3, col in 1:2] 
[ax.xtickformat = values -> ["$(Int(value))" for value in values] for ax in axs_r_h[end, :]]
[ax.ytickformat = Makie.automatic for ax in axs_r_h[:, 1]]
[ax.xminorticks=IntervalsBetween(2) for ax in axs_r_h[:]]
[ax.xminorticksvisible = true for ax in axs_r_h[:]]
#+    limits = (nothing, (1e-4, 1e4)), 
#+    xlabel = "Radial Distance [m]", 
#+    ylabel = "Production [1/m³]",
 
heights = [95, 105, 120, 140, 200, 350]

hm_ion = nothing
hm_prod = nothing

for (i2, collection) in enumerate([runs_hrp_20[2:2:end], runs_hrp_90[2:2:end]])
    for (i1, r) in enumerate(collection)
        #io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", r), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        @time his_hrp = rebin(his_hrp, (
            filter(x -> mod(x, 1e4) == 0, his_hrp.edges[1]),
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        #dv = [dh*dr2*dp /2 for dh in diff(h_edges), dr2 in diff(r_edges.^2), dp in diff(p_edges)]
        #his_hrp_norm = normalize_histogram_density(his_hrp, dv)

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
        #mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
        mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (2 * sin(lim_pitch)^2) #higher precision by avoiding subtraction
        B_top = convergent_vertical_field([0, 0, 600e3+c.re])
        v_perp_mean = mean_vperp_theory * v_abs(E0)
        r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)

        mag_moment = c.me * (v_abs(E0)*sin(lim_pitch))^2 / (2 * norm(B_top))
        v_perp_max_pitch_lim = sqrt.(mag_moment * 2 * B / c.me)
        r_gyro_max_pitch_lim = v_perp_max_pitch_lim * c.me ./ (c.qe * B)

        if lim_pitch_deg == 90
            r_gyro_max_pitch_lim = v_abs(E0) * c.me ./ (c.qe * B)
        end

        if any(v_perp_max_pitch_lim .> v_abs(E0))
            println("Mirroring $lim_pitch_deg deg")
            #error()
        end

        # Radial - Height plot
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)'        
        data_h_normal = data ./ maximum(data, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        i_max = [argmax(d) for d in eachcol(data)]

        ax = axs_ion[i1, i2]
        hm_ion = heatmap!(ax, r_edges, h_edges/1e3, data,
            colorscale = log10,
            colorrange = (1, 1e6 #maximum(data)
            ), lowclip = ("white", 1e-2),
            )
        println( maximum(data))  
        sleep(1)
        if lim_pitch_deg == 90
            text!(ax, 15, 460, text = "$E0 eV", font =:bold)
        end

        data_h_normal_f = mapwindow(median!, data_h_normal, (3, 3))
        data_h_normal_f = data_h_normal

        lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = (:dot, 10), linewidth = 1, label = "Max gyroradius")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(1)], color= "red", linestyle = (:dash, :loose), label = "1/e contour")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(2)], color= "red", linestyle = (:dash, :loose), label = "1/2e contour")

        if showlines
            lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
            contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
            contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
            lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
            lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
        end

        # Radial - Height plot normalised density
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data = data ./ dA
        data_h_normal = data ./ maximum(data, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        i_max = [argmax(d) for d in eachcol(data)]

        ax = axs_prod[i1, i2]
        hm_prod = heatmap!(ax, r_edges, h_edges/1e3, data,
            colorscale = log10,
            colorrange = (1e-4, 1e2 #maximum(data)
            ), lowclip = ("white", 1e-4),
            )
        println( maximum(data))
        sleep(1)
        if lim_pitch_deg == 90
            text!(ax, 15, 460, text = "$E0 eV", font =:bold)
        end


        data_h_normal_f = mapwindow(median!, data_h_normal, (3, 3))
        #data_h_normal_f = data_h_normal

        lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = (:dot, 10), linewidth = 1, label = "Max gyroradius")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(1)], color= "red", linestyle = (:dash, 10), label = "1/e contour")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(2)], color= "red", linestyle = (:dash, 10), label = "1/2e contour")

        if showlines
            lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
            contour!(ax, r_middle, h_middle[1:25]/1e3, data_h_normal[:, 1:25], levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
            lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
            lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
        end 

        ax = axs_r_h[i1, i2]
        d = dropdims(sum(data, dims = 2), dims = 2)
        i_min = argmin(abs.(d .- maximum(d)/exp(1)))
        #l = lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d], label= "Integral", color = Makie.to_colormap(:viridis)[7*36+1])
        l = lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d], label= "Integral", color = "black", linestyle = :dashdot)#Makie.to_colormap(:managua)[7*36+1])
        #s = scatter!(ax, r_middle[i_min], d[i_min])
        for (i3, h) in enumerate(heights)  
            hi = findfirst(x -> x > h * 1e3, h_middle)
            d = data[:, hi]
        #    lines!(ax, [-r_middle; r_middle], [d; d], label= "$h km", color = Makie.to_colormap(:viridis)[(i3-1)*36+1])
            lines!(ax, [-reverse(r_middle); r_middle], [reverse(d); d], label= "$h km", color = Makie.to_colormap(:managua)[(i3-1)*36+1])
            if lim_pitch_deg == 90
                text!(ax, 10, 460, text = "$E0 eV", font =:bold)
            end
        end
    end
end

Colorbar(f_hr_ion[1:3, 3], hm_ion, label="Counts [1]")
Colorbar(f_hr_prod[1:3, 3], hm_prod, label="Production [m⁻³]")

axs_ion[1, 1].title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
axs_ion[1, 2].title = "isotropic"

axs_prod[1, 1].title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
axs_prod[1, 2].title = "isotropic"

axs_r_h[1, 1].title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
axs_r_h[1, 2].title = "isotropic"

[ax.xlabel = "Radial Distance [m]" for ax in axs_ion[end, :]]
[ax.ylabel = "Height [km]" for ax in axs_ion[:, 1]]

[ax.xlabel = "Radial Distance [m]" for ax in axs_prod[end, :]]
[ax.ylabel = "Height [km]" for ax in axs_prod[:, 1]]

[ax.xlabel = "Radial Distance [m]" for ax in axs_r_h[end, :]]
[ax.ylabel = "Production [m⁻³]" for ax in axs_r_h[:, 1]]

axislegend(axs_ion[1, 1])
axislegend(axs_prod[1, 1])
axislegend(axs_r_h[1, 1])

linkyaxes!(axs_ion[:])
linkxaxes!(axs_ion[:])
linkyaxes!(axs_prod[:])
linkxaxes!(axs_prod[:])
linkyaxes!(axs_r_h[:])
linkxaxes!(axs_r_h[:])

display(f_hr_ion)
display(f_hr_prod)
display(f_prod_r_h)

save(joinpath(dir, "plots", "panel_hist_hr_all_E_pitch.png"), f_hr_ion, px_per_unit = 3.3)
save(joinpath(dir, "plots", "panel_hist_norm_hr_all_E_pitch.png"), f_hr_prod, px_per_unit = 3.3)
save(joinpath(dir, "plots", "panel_prod_r_h_all_E_pitch.png"), f_prod_r_h, px_per_unit = 3.3)
##


##
WGLMakie.activate!()

#using CairoMakie
CairoMakie.activate!()



f_prod_r_h = Figure(size = (900, 1800))
sleep(1)
axs_r_h = [Axis(f_prod_r_h[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    limits = ((-20, 20), (2e-3, 2e3)),
    yscale = log10,
    xticks = [-30, -20, -10, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30],
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:7, col in 1:1] 
[ax.xtickformat = values -> ["$(Int(value))" for value in values] for ax in axs_r_h[end, :]]
[ax.ytickformat = Makie.automatic for ax in axs_r_h[:, 1]]
[ax.xminorticks=IntervalsBetween(2) for ax in axs_r_h[:]]
[ax.xminorticksvisible = true for ax in axs_r_h[:]]
#+    limits = (nothing, (1e-4, 1e4)), 
#+    xlabel = "Radial Distance [m]", 
#+    ylabel = "Production [1/m³]",
 
heights = [105, 120, 140, 180, 250, 350, 500]

lines_r_h = []


for (i2, collection) in enumerate([runs_hrp_20, runs_hrp_90])
    for (i1, r) in enumerate(collection)
        #io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", r), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        @time his_hrp = rebin(his_hrp, (
            filter(x -> mod(x, 1e4) == 0, his_hrp.edges[1]),
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        #dv = [dh*dr2*dp /2 for dh in diff(h_edges), dr2 in diff(r_edges.^2), dp in diff(p_edges)]
        #his_hrp_norm = normalize_histogram_density(his_hrp, dv)

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
        #mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
        mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (2 * sin(lim_pitch)^2) #higher precision by avoiding subtraction
        B_top = convergent_vertical_field([0, 0, 600e3+c.re])
        v_perp_mean = mean_vperp_theory * v_abs(E0)
        r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)

        mag_moment = c.me * (v_abs(E0)*sin(lim_pitch))^2 / (2 * norm(B_top))
        v_perp_max_pitch_lim = sqrt.(mag_moment * 2 * B / c.me)
        r_gyro_max_pitch_lim = v_perp_max_pitch_lim * c.me ./ (c.qe * B)

        if lim_pitch_deg == 90
            r_gyro_max_pitch_lim = v_abs(E0) * c.me ./ (c.qe * B)
        end

        if any(v_perp_max_pitch_lim .> v_abs(E0))
            println("Mirroring $lim_pitch_deg deg")
            #error()
        end


        # Radial - Height plot normalised density
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data = data ./ dA
        data_h_normal = data ./ maximum(data, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        i_max = [argmax(d) for d in eachcol(data)]
 

        ax = axs_r_h[i1, 1]
        d = dropdims(sum(data, dims = 2), dims = 2)
        i_min = argmin(abs.(d .- maximum(d)/exp(1)))
        lines_int = []
        if lim_pitch_deg == 20
            l = lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d], label= "Integral", color = Makie.to_colormap(:turbo)[7*36+1])
            append!(lines_int, [l])
            #hlines!(ax, maximum(d)/exp(1))
            l = vlines!(ax, [-r_middle[i_min], r_middle[i_min]])
            append!(lines_int, [l])
        else
            l = lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d],
                label= "Integral", color = Makie.to_colormap(:turbo)[7*36+1],
                linestyle = :dot)
            append!(lines_int, [l])
            #hlines!(ax, maximum(d)/exp(1), linestyle = :dot)
            l = vlines!(ax, [-r_middle[i_min], r_middle[i_min]], linestyle = :dot)
            append!(lines_int, [l])
        end
        s = scatter!(ax, r_middle[i_min], d[i_min])
        append!(lines_int, [s])
        append!(lines_r_h, [lines_int])
        for (i3, h) in enumerate(heights)  
            hi = findfirst(x -> x > h * 1e3, h_middle)
            d = data[:, hi]
            #lines!(ax, [-r_middle; r_middle], [d; d], label= "$h km", color = Makie.to_colormap(:turbo10)[i3+3])
            if lim_pitch_deg == 20
                l =  lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d], label= "$h km", color = Makie.to_colormap(:turbo)[(i3-1)*36+1])
                append!(lines_r_h, [l])
            else
                l = lines!(ax, [reverse(r_middle).*(-1); r_middle], [reverse(d); d], label= "$h km", color = Makie.to_colormap(:turbo)[(i3-1)*36+1], linestyle = :dash)
                append!(lines_r_h, [l])
            end
            if lim_pitch_deg == 90
                text!(ax, 15, 460, text = "$E0 eV", font =:bold)
            end
        end
    end
end

[ax.xlabel = "Radial Distance [m]" for ax in axs_r_h[end, :]]
[ax.ylabel = "Production [m⁻³]" for ax in axs_r_h[:, 1]]

a = vcat([[i8 == 1 ? lines_r_h[i8 + ii*8] : lines_r_h[i8 + ii*8] for ii in 0:9] for i8 in 1:8][1]...)
leg = Legend(f_prod_r_h[:, 2],[i8 == 1 ? a : [lines_r_h[i8 + ii*8] for ii in 0:9] for i8 in 1:8], ["Integral"; ["$h km" for h in heights]])
#axislegend(axs_r_h[1, 1])
linkyaxes!(axs_r_h[:])
linkxaxes!(axs_r_h[:])

display(f_prod_r_h)


##
#CairoMakie.activate!()
fig = Figure()
sleep(1)
ax = Axis(fig[1, 1],
    limits = (sqrt.((3e2, 35e3)./1e3), (-0.8, 31)),
    xticks = sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]./1e3),
    #ax.title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
    xtickformat = values -> [rich(value, superscript("1/2")) for value in ["0.5", "1", "2", "4", "8", "16", "32"]], 
    #xscale = log10,
    xlabel = rich("(Energy [keV])", superscript("1/2")),
    yticks = 0:5:30,
    ylabel = "Width [m]"
    )
width_20 = zeros(size(runs_hrp_20, 1))
width_90 = zeros(size(runs_hrp_20, 1))
err_width_20 = zeros(size(runs_hrp_20, 1))
err_width_90 = zeros(size(runs_hrp_20, 1))
for (i2, collection) in enumerate([runs_hrp_20, runs_hrp_90])
    for (i1, r) in enumerate(collection)
        #io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        #println(r)
        io = open(joinpath(dir, "hist_summed", r), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
    
        #extract edges
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        h_middle = h_edges[1:end-1] + diff(h_edges)/2
        r_middle = r_edges[1:end-1] + diff(r_edges)/2
        p_middle = p_edges[1:end-1] + diff(p_edges)/2

        data = dropdims(sum(his_hrp.weights, dims = (1,3)), dims = (1,3))
        dh = h_edges[end] - h_edges[1]
        dA = [dh*dr2*pi for dr2 in diff(r_edges.^2)]
        d = data ./ dA
        #lines([reverse(r_middle).*(-1); r_middle], [reverse(d); d])
        #hlines!(ax, maximum(d)/exp(1))
        i_min = argmin(abs.(d .- maximum(d)/exp(1)))
        println(E0, " ", lim_pitch_deg, " ", r_middle[i_min])
        #scatter!(ax, r_middle[i_min], d[i_min])
        #vlines!(ax, [-r_middle[i_min], r_middle[i_min]])

        d = data ./ dA
deltad = sqrt.(data) ./dA
d_max = maximum(d)
deltad_max = deltad[argmax(d)]
d1e = maximum(d) / exp(1)
i_min = argmin(abs.(d .- maximum(d) / exp(1)))
using Interpolations
r_itp = linear_interpolation(d[i_min+1:-1:i_min-1], r_middle[i_min+1:-1:i_min-1])
w = r_itp(maximum(d) / exp(1))
dddr = (d[i_min-1] - d[i_min+1])/(r_middle[i_min-1] - r_middle[i_min+1])
deltaw = abs(deltad_max / dddr)

        println(E0, " ", lim_pitch_deg, " ", w, " ", deltaw)


        if lim_pitch_deg == 20
            #scatter!(ax, sqrt(E0/1e3), r_middle[i_min], color = Makie.wong_colors()[1])
            width_20[i1] =  r_middle[i_min]
            width_20[i1] =  w
            err_width_20[i1] =  deltaw
        else
            #scatter!(ax, sqrt(E0/1e3), r_middle[i_min], color = Makie.wong_colors()[2])
            width_90[i1] =  r_middle[i_min]
            width_90[i1] =  w
            err_width_90[i1] =  deltaw
        end
    end
end
scatterlines!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]./1e3), 2 .* width_20, color = Makie.wong_colors()[1],label = L"{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}")
errorbars!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]./1e3), 2 .* width_20, err_width_20, color = Makie.wong_colors()[1])
scatterlines!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]./1e3), 2 .* width_90, color = Makie.wong_colors()[2],label = "isotropic")
errorbars!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]./1e3), 2 .* width_90, err_width_90, color = Makie.wong_colors()[2])
scatter!(ax, sqrt.([5, 25]), [9.5, 22.4], color = Makie.wong_colors()[3],label = "Borovsky 1991", marker = :star5)
axislegend(ax, position = :lt)
#display(fig)
save(joinpath(dir, "plots", "width_h_int.png"), fig, px_per_unit = 3.3)

using LsqFit

function linear(x, a)
    return a[1] .* x .+ a[2]
end

fit_90 = curve_fit(linear, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]), 2*width_90, [1.0, 1])
fit_20 = curve_fit(linear, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3, 16e3, 32e3]), 2*width_20, [1.0, 1])

x = [0.1, 33]*1e3
lines!(ax, sqrt.(x), linear(sqrt.(x), fit_90.param))
lines!(ax, sqrt.(x), linear(sqrt.(x), fit_20.param))
##

r = "h_hrp_16000.0eV_20.0deg_summed.hist"
io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)

#extract edges
h_edges = his_hrp.edges[1]
r_edges = his_hrp.edges[2]
p_edges = his_hrp.edges[3]
h_middle = h_edges[1:end-1] + diff(h_edges)/2
r_middle = r_edges[1:end-1] + diff(r_edges)/2
p_middle = p_edges[1:end-1] + diff(p_edges)/2
data = dropdims(sum(his_hrp.weights, dims = (1,3)), dims = (1,3))
dh = h_edges[end] - h_edges[1]
dA = [dh*dr2*pi for dr2 in diff(r_edges.^2)]

d = data ./ dA
deltad = sqrt.(data) ./dA
d_max = maximum(d)
deltad_max = deltad[argmax(d)]
d1e = maximum(d) / exp(1)
i_min = argmin(abs.(d .- maximum(d) / exp(1)))
using Interpolations
r_itp = linear_interpolation(d[i_min+1:-1:i_min-1], r_middle[i_min+1:-1:i_min-1])
w = r_itp(maximum(d) / exp(1))
dddr = (d[i_min-1] - d[i_min+1])/(r_middle[i_min-1] - r_middle[i_min+1])
deltaw = abs(deltad_max / dddr)

lines(r_middle, data)
errorbars(r_middle, data, sqrt.(data), #axis = (limits = ((39, 40), (0, 5e2)),),
)

errorbars(r_middle, d, deltad, #axis = (limits = ((6, 7), (1.8, 1.9)),),
)


lines!(r_middle, abs.(d .- maximum(d) / exp(1)))

vlines!(r_middle[i_min])

using Interpolations
r_itp = linear_interpolation(d[i_min+1:-1:i_min-1], r_middle[i_min+1:-1:i_min-1])
w = r_itp(maximum(d) / exp(1))
dddr = (d[i_min-1] - d[i_min+1])/(r_middle[i_min-1] - r_middle[i_min+1])
deltaw = abs(deltad_max / dddr)
scatter!(r_itp(maximum(d) / exp(1)), maximum(d) / exp(1))
lines!(r_itp(d[i_min+1]:0.001:d[i_min-1]), d[i_min+1]:0.001:d[i_min-1], linestyle = :dash)
errorbars!([w], [d1e], [deltaw], direction =:x)

##
WGLMakie.activate!()
CairoMakie.activate!()

#for (i2, collection) in enumerate(zip(runs_hrp_20, runs_hrp_90))
for (i2, collection) in enumerate([runs_hrp_20, runs_hrp_90])
    println(collection)
    fig = Figure()
    ax = Axis(fig[1, 1])
    for (i1, r) in enumerate(collection)
        println(r)
#        io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", r), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        
        """
        @time his_hrp = rebin(his_hrp, (
#            filter(x -> mod(x, 1e3) == 0, his_hrp.edges[1]),
            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:249e3; 250e3:5e3:299e3; 300e3:5e3:600e3],
#            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:600e3; 600e3],
#            [80e3:1e3:359e3; 360e3:2e3:600e3],
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))
        """

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        
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
        #mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
        mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (2 * sin(lim_pitch)^2) #higher precision by avoiding subtraction
        B_top = convergent_vertical_field([0, 0, 600e3+c.re])
        v_perp_mean = mean_vperp_theory * v_abs(E0)
        r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)

        mag_moment = c.me * (v_abs(E0)*sin(lim_pitch))^2 / (2 * norm(B_top))
        v_perp_max_pitch_lim = sqrt.(mag_moment * 2 * B / c.me)
        r_gyro_max_pitch_lim = v_perp_max_pitch_lim * c.me ./ (c.qe * B)

        if lim_pitch_deg == 90
            r_gyro_max_pitch_lim = v_abs(E0) * c.me ./ (c.qe * B)
        end

        if any(v_perp_max_pitch_lim .> v_abs(E0))
            println("Mirroring $lim_pitch_deg deg")
            #error()
        end

        # Radial - Height plot
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data = data ./ dA
        data_f = mapwindow(median!, data, (3, 3))
        data_f = mapwindow(median!, data_f, (3, 3))
        #data_f = data
        data_h_normal = data_f ./ maximum(data_f, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        #data_f = mapwindow(median!, data_h_normal, (3, 3))

        i_max = [argmax(d) for d in eachcol(data_h_normal)]

        #ax = axs_prod[i1, i2]
        #fig, ax, hm_prod = heatmap(r_edges, h_edges/1e3, data,
        #    axis=(limits = ((0, 10), nothing),),
        #    colorscale = log10,
        #    colorrange = (1e-4, 1e2), lowclip = ("white", 1e-4),
        #    colormap = :turbo,
        #    alpha = 0.35
        #    )
        
        #sleep(1)
        #ax2 = Axis(fig[1, 2])
        #sleep(1)
        #heatmap!(ax2, r_edges, h_edges/1e3, data_h_normal,
        #    colorscale = log10,
        #    colorrange = (1e-4, 1 #maximum(data)
        #    ), lowclip = ("white", 1e-4),
        #    colormap = :turbo
        #    )
        #linkxaxes!(ax, ax2)
        #linkyaxes!(ax, ax2)

        println( maximum(data))
        sleep(1)
        #if lim_pitch_deg == 90
        #    text!(ax, 15, 460, text = "$E0 eV", font =:bold)
        #end
        #text!(ax, 6, 460, text = "$E0 eV, $lim_pitch_deg", font =:bold)
        
        #lines!(ax, r_middle[i_max], h_middle/1e3)
        #cont = contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal[:])/exp(n) for n in 1:15], 
        dc = copy(data_h_normal)
        [dc[1:im, i1] .= 1 for (i1, im) in enumerate(i_max)]

        tot_ion = dropdims(sum(data, dims = 1), dims = 1)
        #scatter!(ax,tot_ion, h_middle/1e3)

        mask = BitVector(zeros(size(tot_ion, 1)))# .> 1
        mask[1:100] = tot_ion[1:100] .< 3
        mask[250:end] = tot_ion[250:end] .< 0.2
        #dc[:, mask] .= undef

        color = Makie.to_colormap(:batlow)[256 - (i1-1)*41]
        cont = contour!(ax, r_middle, h_middle[.!mask]/1e3, dc[:, .!mask], levels=[maximum(dc)/exp(1)], color = color, linewidth = 0)
        #cont = contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], 
        #    color= "red", linestyle = :dash, label = "1/e contour")
        x = [d[1] for d in cont.contour_points.value.x]
        y = [d[2] for d in cont.contour_points.value.x]
        lines!(ax, x[.!isnan.(x)], y[.!isnan.(y)], label = "$E0 eV", color = color)
        lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color = color, linestyle = :dot)

        xlims!(0, 17)
        ylims!(80, 400)

    
        #expecv = data


    end

    if lim_pitch_deg ==20
        ax.title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
    else
        ax.title = "isotropic"
    end
    ax.xlabel = "Radial Distance [m]" 
    ax.ylabel = "Height [km]"
    lines!(ax, 0, 0, color = "black",  linestyle = :dot, label="Min. Gyroradius")
    ax.xticks = 0:2:20
    axislegend(ax)
    display(fig)
    if lim_pitch_deg ==20
        save(joinpath(dir, "plots", "contour_allE_20deg.png"), fig, px_per_unit = 3.3)
    else
        save(joinpath(dir, "plots", "contour_allE_90deg.png"), fig, px_per_unit = 3.3)
    end
end



##
fig = Figure(size = (800, 500))


for (i2, collection) in enumerate(zip(runs_hrp_20[1:2:end], runs_hrp_90[1:2:end]))
    #println(collection)
    if i2 == 1
        ax = Axis(fig[1, i2])
    else
        ax = Axis(fig[1, i2], ytickformat = values -> ["" for value in values],)
    end
    for (i1, r) in enumerate(collection)
        #println(r)
#        io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", r), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        
        """
        @time his_hrp = rebin(his_hrp, (
#            filter(x -> mod(x, 1e3) == 0, his_hrp.edges[1]),
            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:249e3; 250e3:5e3:299e3; 300e3:5e3:600e3],
#            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:600e3; 600e3],
#            [80e3:1e3:359e3; 360e3:2e3:600e3],
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))
        """

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        
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
        #mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
        mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (2 * sin(lim_pitch)^2) #higher precision by avoiding subtraction
        B_top = convergent_vertical_field([0, 0, 600e3+c.re])
        v_perp_mean = mean_vperp_theory * v_abs(E0)
        r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)

        mag_moment = c.me * (v_abs(E0)*sin(lim_pitch))^2 / (2 * norm(B_top))
        v_perp_max_pitch_lim = sqrt.(mag_moment * 2 * B / c.me)
        r_gyro_max_pitch_lim = v_perp_max_pitch_lim * c.me ./ (c.qe * B)

        if lim_pitch_deg == 90
            r_gyro_max_pitch_lim = v_abs(E0) * c.me ./ (c.qe * B)
        end

        if any(v_perp_max_pitch_lim .> v_abs(E0))
            #println("Mirroring $lim_pitch_deg deg")
            #error()
        end

        # Radial - Height plot
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data = data ./ dA
        data_f = mapwindow(median!, data, (3, 3))
        data_f = mapwindow(median!, data_f, (3, 3))
        #data_f = data
        data_h_normal = data_f ./ maximum(data_f, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        #data_f = mapwindow(median!, data_h_normal, (3, 3))

        i_max = [argmax(d) for d in eachcol(data_h_normal)]

        #ax = axs_prod[i1, i2]
        #fig, ax, hm_prod = heatmap(r_edges, h_edges/1e3, data,
        #    axis=(limits = ((0, 10), nothing),),
        #    colorscale = log10,
        #    colorrange = (1e-4, 1e2), lowclip = ("white", 1e-4),
        #    colormap = :turbo,
        #    alpha = 0.35
        #    )
        
        #sleep(1)
        #ax2 = Axis(fig[1, 2])
        #sleep(1)
        #heatmap!(ax2, r_edges, h_edges/1e3, data_h_normal,
        #    colorscale = log10,
        #    colorrange = (1e-4, 1 #maximum(data)
        #    ), lowclip = ("white", 1e-4),
        #    colormap = :turbo
        #    )
        #linkxaxes!(ax, ax2)
        #linkyaxes!(ax, ax2)

        #println( maximum(data))
        #sleep(1)
        #if lim_pitch_deg == 90
        #    text!(ax, 15, 460, text = "$E0 eV", font =:bold)
        #end
        #text!(ax, 6, 460, text = "$E0 eV, $lim_pitch_deg", font =:bold)
        
        #lines!(ax, r_middle[i_max], h_middle/1e3)
        #cont = contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal[:])/exp(n) for n in 1:15], 
        dc = copy(data_h_normal)
        [dc[1:im, i1] .= 1 for (i1, im) in enumerate(i_max)]

        tot_ion = dropdims(sum(data, dims = 1), dims = 1)
        i_max_totion = argmax(tot_ion)

        #scatter!(ax,tot_ion, h_middle/1e3)

        mask = BitVector(zeros(size(tot_ion, 1)))# .> 1
        mask[1:100] = tot_ion[1:100] .< 3
        mask[250:end] = tot_ion[250:end] .< 0.2
        #dc[:, mask] .= undef

        color = Makie.wong_colors()[i1]
        #color = Makie.to_colormap(:batlow)[(i1-1)*123+1]
        cont = contour!(ax, r_middle, h_middle[.!mask]/1e3, dc[:, .!mask], levels=[maximum(dc)/exp(1)], color = color, linewidth = 0)
        #cont = contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], 
        #    color= "red", linestyle = :dash, label = "1/e contour")
        x = [d[1] for d in cont.contour_points.value.x]
        y = [d[2] for d in cont.contour_points.value.x]

        x = x[.!isnan.(x)]
        y = y[.!isnan.(y)]

        y = y[x .> minimum(r_gyro_max_pitch_lim)]
        x = x[x .> minimum(r_gyro_max_pitch_lim)]
        #lines!(ax, x[.!isnan.(x)], y[.!isnan.(y)], label = "$E0 eV", color = color)
        lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color = color, linestyle = :dot)
        
        nearest_index = argmin(abs.(y[.!isnan.(y)] .- h_middle[i_max_totion]/1e3))
        println(E0, " ", lim_pitch_deg, " ", x[nearest_index])

        if i1 ==1
            lines!(ax, x[.!isnan.(x)], y[.!isnan.(y)], label = "f.a.", color = color)
            scatter!(ax, x[nearest_index], y[nearest_index], color = color, marker = :utriangle)
            #hlines!(ax, h_middle[i_max_totion]/1e3)
        else 
            lines!(ax, x[.!isnan.(x)], y[.!isnan.(y)], label = "iso.", color = color)
            scatter!(ax, x[nearest_index], y[nearest_index], color = color, marker = :dtriangle)
        end

        if i2 == 1
            ax.ylabel = "Height [km]"
        end
    
        #expecv = data


    end

    xlims!(0, 2 + i2^2 )
    ylims!(80, 400)

    ax.title = "$E0 eV"
    #
    #if lim_pitch_deg ==20
    #    ax.title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
    #else
    #    ax.title = "isotropic"
    #end
    
    if i2 == 1
        lines!(ax, 0, 0, color = "black",  linestyle = :dot, label="rₘᵢₙ")
        axislegend(ax)
    end
    ax.xticks = 0:i2:20    
end

Label(fig[2, 1:end], "Radial Distance [m]")
display(fig)
save(joinpath(dir, "plots", "panel_contour_allE.png"), fig, px_per_unit = 3.3)

##
WGLMakie.activate!()
CairoMakie.activate!()

#for (i2, collection) in enumerate(zip(runs_hrp_20, runs_hrp_90))
#for (i2, collection) in enumerate(zip(runs_hrp_20, runs_hrp_90))
collection = ("h_hrp_500.0eV_20.0deg_summed.hist", "h_hrp_500.0eV_90.0deg_summed.hist")
    println(collection)
    
    r_fa = collection[1]
#        io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", r_fa), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        
        """
        @time his_hrp = rebin(his_hrp, (
#            filter(x -> mod(x, 1e3) == 0, his_hrp.edges[1]),
            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:249e3; 250e3:5e3:299e3; 300e3:5e3:600e3],
#            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:600e3; 600e3],
#            [80e3:1e3:359e3; 360e3:2e3:600e3],
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))
        """

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        
        #extract edges
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        h_middle = h_edges[1:end-1] + diff(h_edges)/2
        r_middle = r_edges[1:end-1] + diff(r_edges)/2
        p_middle = p_edges[1:end-1] + diff(p_edges)/2

        # Radial - Height plot
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data_fa = data ./ dA

        riso = collection[2]
#        io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
        io = open(joinpath(dir, "hist_summed", riso), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
        close(io)
        
        """
        @time his_hrp = rebin(his_hrp, (
#            filter(x -> mod(x, 1e3) == 0, his_hrp.edges[1]),
            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:249e3; 250e3:5e3:299e3; 300e3:5e3:600e3],
#            [80e3:1e3:119e3; 120e3:2e3:139e3; 140e3:3e3:600e3; 600e3],
#            [80e3:1e3:359e3; 360e3:2e3:600e3],
            his_hrp.edges[2],
            his_hrp.edges[3],
            ))
        """

        # compute bin volumes for normalization
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        
        #extract edges
        h_edges = his_hrp.edges[1]
        r_edges = his_hrp.edges[2]
        p_edges = his_hrp.edges[3]
        h_middle = h_edges[1:end-1] + diff(h_edges)/2
        r_middle = r_edges[1:end-1] + diff(r_edges)/2
        p_middle = p_edges[1:end-1] + diff(p_edges)/2

        # Radial - Height plot
        data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
        data_iso = data ./ dA

        data_fa = mapwindow(median!, data_fa, (3, 3))
        data_fa = mapwindow(median!, data_fa, (3, 3))
        data_iso = mapwindow(median!, data_iso, (3, 3))
        data_iso = mapwindow(median!, data_iso, (3, 3))
        ratio_fa_iso = data_fa ./ data_iso
        diff_fa_iso = data_fa .- data_iso

        fratio, ax, hm_rat = heatmap(r_edges, h_edges/1e3, ratio_fa_iso,
            axis=(limits = ((0, 20), nothing),),
            colorscale = log10,
            colorrange = (1/30,30), #highclip = ("white", 1),
            colormap = :BrBG_5,
            )
        Colorbar(fratio[1, 2], hm_rat, label="Ratio prod [1]")


        fdiff, ax, hm_diff = heatmap(r_edges, h_edges/1e3, asinh.(diff_fa_iso),
            axis=(limits = ((0, 20), nothing),),
            #colorscale = :Symlog10,
            colorrange = asinh.((-1, 1).*maximum(abs.(diff_fa_iso[:]))), lowclip = ("white", 1/30), highclip = ("white", 1/30),
            colormap = :balance,
            )
        Colorbar(fdiff[1, 2], hm_diff, label="Diff prod [m-3]")

        display(fratio)
        display(fdiff)
end

 
##
        #contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "black", linestyle = :dash, label = "1/e contour")
        x = [d[1] for d in cont.contour_points.value.x]
        y = [d[2] for d in cont.contour_points.value.x]
        inan = findall(x-> isnan(x), x)
        i1 = 0
        for ina in axes(inan, 1)
            if y[inan[ina]+1] < 190 && y[inan[ina+1]-1]>190
                i1 = ina
                break
            end
        end

        lines!(ax, x[inan[i1]+1:inan[i1+1]-1], y[inan[i1]+1:inan[i1+1]-1])
        lines!(ax, x[1:inan[i1]-1], y[1:inan[i1]-1])

        x[isnan.(x)] .= 0
        y[isnan.(y)] .= 0
        lines!(ax, x, y)

        ih = 30
         r_75p = zeros(size(h_middle))
        for ih in axes(h_middle, 1)
            cdf = cumsum(data_h_normal[:, ih])
            cdf = cdf ./ cdf[end]
        #    lines(r_middle, cdf)
            ir_75p = findfirst(x -> x>0.75, cdf)
            if isnothing(ir_75p) ir_75p = 1 end
        #    scatter!(r_middle[ir_75p], cdf[ir_75p])
            r_75p[ih] = r_middle[ir_75p]
        end
        lines!(ax, r_75p, h_middle/1e3)


##



#dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
io = open(joinpath(dir, "hist_summed", "h_xyz_1000.0eV_90.0deg_summed.hist"), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
close(io)
x_edges = his_xyz.edges[1]
y_edges = his_xyz.edges[2]
z_edges = his_xyz.edges[3]
dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
his_xyz_norm = normalize_histogram_density(his_xyz, dv)


x_edges = his_xyz.edges[1]
y_edges = his_xyz.edges[2]
z_edges = his_xyz.edges[3]
x_middle = x_edges[1:end-1] + diff(x_edges)/2
y_middle = y_edges[1:end-1] + diff(y_edges)/2
z_middle = z_edges[1:end-1] + diff(z_edges)/2




##
#single plots

#for r in runs_hrp
io = open(joinpath(dir, "hist_summed", "h_hrp_8000.0eV_20.0deg_summed.hist"), "r")
#io = open(joinpath(dir, "hist_summed", r), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)
@time his_hrp = rebin(his_hrp, (
    filter(x -> mod(x, 1e4) == 0, his_hrp.edges[1]),
    his_hrp.edges[2],
    his_hrp.edges[3],
    ))

# compute bin volumes for normalization
h_edges = his_hrp.edges[1]
r_edges = his_hrp.edges[2]
p_edges = his_hrp.edges[3]
dv = [dh*dr*dp /2 for dh in diff(h_edges), dr in diff(r_edges.^2), dp in diff(p_edges)]
his_hrp_norm = normalize_histogram_density(his_hrp, dv)




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
#mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))
mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (2 * sin(lim_pitch)^2) #higher precision by avoiding subtraction
B_top = convergent_vertical_field([0, 0, 600e3+c.re])
v_perp_mean = mean_vperp_theory * v_abs(E0)
r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)

r_gyro_pitch_lim = v_abs(E0) * sin(lim_pitch) * c.me ./ (c.qe * B)

h_ind = 7
#change h_ind to heights, and find correspdng h_ind


## production at selected altitudes:
fig = Figure()
sleep(1)
ax = Axis(fig[1, 1], 
    title = "$E0 eV, $lim_pitch_deg deg", 
    limits = (nothing, (1e-4, 1e4)), 
    xlabel = "Radial Distance [m]", 
    ylabel = "Production [1/m³]",
    yscale = log10,
    )


heights = [105, 120, 140, 180, 250, 350, 500]
h_inds = [findfirst(h -> h > height * 1e3, h_middle) for height in heights]
for (i1, h) in enumerate(heights)  
    hi = findfirst(x -> x > h * 1e3, h_middle)
    d = dropdims(sum(his_hrp_norm.weights[hi, :, :], dims = 2), dims = 2)
    lines!(ax, [-r_middle; r_middle], [d; d], label= "$h km")
end
axislegend()
fig


fig = Figure()

 



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
    rlims!(ax, 0.0, 5.0)
    #l1 = lines!(ax, p_middle, r_middle[i_max],
    #    color= "red", label = "Max Counts")
    #l2 = contour!(ax, p_middle, r_middle, data, 
    #    levels=[maximum(data)/exp(1)], color= "red", 
    #    linestyle = :dash, label = "1/e contour")
    #l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), 
    #    color= "black", linestyle = :dot, label = "Mean Gyroradius")
    Colorbar(fig[2, 2], p, label="Counts [1]")
    #Legend(fig[1, 2], [l1, l2, l4], ["Max Counts", "1/e contour", "Mean Gyroradius"])
    ax.titlegap = 20
    #fig
    #save(joinpath(dir, "plots", "hist3d_phase_rd_$(E0)_$(lim_pitch_deg)_$(h_middle[h_ind]/1e3)km.png"), fig)
    display(fig)
end

#find h_inds for heights 100, 120, 150, 200 km
h_inds = [findfirst(h -> h > height * 1e3, h_middle) for height in [100, 120, 150, 200]]

for h_ind in [8]#, 5, 8, 13, 23]
    plot_histogram_phase_rdist(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
end


function plot_histogram_phase_rdist_normalised(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp, r_gyro_mean, h_ind)
    #r_edges = his_hrp.edges[2]
    #p_edges = his_hrp.edges[3]
    #area_pr = [dr * dp / 2  for dp in diff(p_edges), dr in (r_edges[2:end] .^2 - r_edges[1:end-1] .^2)]

    data = his_hrp.weights[h_ind, :, :]' #./ area_pr
    i_max = [argmax(d) for d in eachrow(data)]
    min = minimum(his_hrp.weights[his_hrp.weights .!= 0])
    
    fig = Figure()
    sleep(5)
    ax = PolarAxis(fig[1:2, 1], title = "Production\n$E0 eV, $lim_pitch_deg deg, $(h_middle[h_ind]/1e3)km",)
    p = voronoiplot!(ax, p_middle, r_middle, data,
        colorscale = log10,
        colorrange = (min, max(min, maximum(data))), lowclip = ("white", 0),
        show_generators = false, strokewidth = 0,)
    rlims!(ax, 0.0, 10.0)
    #l1 = lines!(ax, p_middle, r_middle[i_max],
    #    color= "red", label = "Max Production")
    #l2 = contour!(ax, p_middle, r_middle, data, 
    #    levels=[maximum(data)/exp(1)], color= "red", 
    #    linestyle = :dash, label = "1/e contour")
    #l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), 
    #    color= "black", linestyle = :dot, label = "Mean Gyroradius")
    Colorbar(fig[2, 2], p, label="Production [m⁻³]")
    #Legend(fig[1, 2], [l1, l2, l4], ["Max Production", "1/e contour", "Mean Gyroradius"])
    ax.titlegap = 20
    #fig
    save(joinpath(dir, "plots", "hist3d_phase_rd_norm$(E0)_$(lim_pitch_deg)_$(h_middle[h_ind]/1e3)km.png"), fig)
    display(fig)
end



for h_ind in [8]#, 5, 8, 13, 23]
    plot_histogram_phase_rdist_normalised(dir, E0, lim_pitch_deg, h_middle, r_middle, p_middle, his_hrp_norm, r_gyro_mean, h_ind)
end

##


## Radial - Height plot
data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)'
data_h_normal = data ./ maximum(data, dims = 1)
data_h_normal[isnan.(data_h_normal)] .= 0.0
i_max = [argmax(d) for d in eachcol(data)]

fig, ax, hm = heatmap(r_edges, h_edges/1e3, data,
    colorscale = log10,
    colorrange = (1e-2, maximum(data)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 5), (50, 600)),)
    )
sleep(1)
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
lines!(ax, r_gyro_pitch_lim, h_middle/1e3, color= "black", linestyle = :dot, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Counts [1]")
axislegend(ax)
display(fig)
save(joinpath(dir, "plots", "hist2d_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)



## Radial - Height plot normalised density
r_edges = his_hrp.edges[2]
p_edges = his_hrp.edges[3]
area_pr = [dr * dp / 2  for dp in diff(p_edges), dr in (r_edges[2:end] .^2 - r_edges[1:end-1] .^2)]
## Radial - Height plot normalised density
data = dropdims(sum(his_hrp_norm.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
data_h_normal = data ./ maximum(data, dims = 1)
data_h_normal[isnan.(data_h_normal)] .= 0.0
i_max = [argmax(d) for d in eachcol(data)]

fig, ax, hm = heatmap(r_edges, h_edges/1e3, data,
    colorscale = log10,
    colorrange = (1e-2, maximum(data)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 5), (50, 600)),)
    )
sleep(1)
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
contour!(ax, r_middle, h_middle[1:25]/1e3, data_h_normal[:, 1:25], levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
#lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
lines!(ax, r_gyro_pitch_lim, h_middle/1e3, color= "black", linestyle = :dot, label = "Max Gyroradius Pitch")
Colorbar(fig[1, 2], hm, label="Density [1/m²]")
axislegend(ax)
display(fig)
save(joinpath(dir, "plots", "hist2d_rd_height_norm$(E0)_$(lim_pitch_deg).png"), fig)

end



#dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
io = open(joinpath(dir, "hist_summed", "h_xyz_1000.0eV_90.0deg_summed.hist"), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
close(io)

his_xyz = rebin(his_xyz, (
            his_xyz.edges[1],
            his_xyz.edges[2],
            filter(x -> mod(x, 1e4) == 0, his_xyz.edges[3]),
            ))

x_edges = his_xyz.edges[1]
y_edges = his_xyz.edges[2]
z_edges = his_xyz.edges[3]
heatmap(x_middle, y_middle, his_xyz.weights[:, :, 8], colorscale = log10,colorrange = (1e-2, maximum(his_xyz.weights[:, :, 8])), lowclip = ("white", 1e-2),)

dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
his_xyz_norm = normalize_histogram_density(his_xyz, dv)


x_edges = his_xyz.edges[1]
y_edges = his_xyz.edges[2]
z_edges = his_xyz.edges[3]
x_middle = x_edges[1:end-1] + diff(x_edges)/2
y_middle = y_edges[1:end-1] + diff(y_edges)/2
z_middle = z_edges[1:end-1] + diff(z_edges)/2

fig, ax, l = lines(0, 0, axis=(yscale=log10,),)
[lines!(ax, x_middle, his_xyz_norm.weights[:, 21, i3]) for i3 in 1:50]
display(fig)