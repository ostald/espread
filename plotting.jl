using Serialization
using DataFrames
using LinearAlgebra


using Bonito    
Bonito.set_cleanup_time!(1)
# ssh -L 9384:localhost:9384 user@server

using WGLMakie
using CairoMakie
CairoMakie.activate!()
WGLMakie.activate!()

include("analysis_util.jl")

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
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



#_______________________________________________________________________________________
## Plot ionization vs height for different energies and pitch angles
#check atmosphere
df = CSV.read("/nfs/revontuli/data/oliver/espread/results/r4_conicB_2025-09-05T14:19:27.566/atmosphere.txt",
            DataFrame;
            header = false,
            delim = " ", 
            types=Float64,
            ignorerepeated=true)
rename!(df, :Column1 => :height, :Column2 => :nN2, :Column3 => :nO2, :Column4 => :nO)
fig, ax, lin = lines(log10.(df.nN2), df.height)
lines!(ax, log10.(df.nO2), df.height)
lines!(ax, log10.(df.nO), df.height)
include("get_msis.jl")
z = hmin+hintervals/2:hintervals:600e3
densityf = make_densityf(hmin+hintervals/2, hmax, hintervals, [69.58, 19.23])
atm = stack(densityf.(z))'
lines!(ax, log10.(atm[:, 1]), z)
lines!(ax, log10.(atm[:, 2]), z)
lines!(ax, log10.(atm[:, 3]), z)
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

using CSV
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
for (i,r) in enumerate(runs_xyz_90)
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
save(joinpath(dir, "plots", "comp_prod_hist_height_xyz_90deg_4000ev.png"), f, px_per_unit = 3.3)
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
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))
    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    df = reference_prod[i1]
    lines!(ax, df.fang*1e6./data, df.height)
    #display(f)
end

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
for r in runs_xyz_20
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
    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_20deg_allE.png"), f, px_per_unit = 3.3)
f

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
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
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

f = Figure()
sleep(1)
ax = Axis(f[1, 1], 
        xlabel = "Production [m⁻³]",
        ylabel = "Height [km]",
        limits = ((1e-1, 1e5), (80, 600)),
        xscale = log10,
        title = "Production vs Height"
        )
for (i, r) in enumerate(runs_xyz)
    println(r)
    dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    
    new_edges = (his_xyz.edges[1], his_xyz.edges[2], filter(x -> mod(x, 1000) == 0, his_xyz.edges[3]))

    his_xyz = rebin(his_xyz, new_edges);
    his_xyz = normalize(his_xyz, mode=:density)
    
    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   
    
    color = Makie.wong_colors()[ceil(Int, i/2)]
    if lim_pitch_deg == 20.0
        lines!(ax, data, z_middle/1e3, label = "$E0 eV", color = color)
    else
        lines!(ax, data, z_middle/1e3, linestyle = :dot, color = color)
    end
end
lines!(ax, [0, 0], [0, 0], color = "black", label = "field-aligned θₗᵢₘ = 20 deg")# lim \theta_{lim} = 20^\circ")
lines!(ax, [0, 0], [0, 0], color = "black", linestyle = :dot, label = "isotropic")
axislegend(ax)
save(joinpath(dir, "plots", "hist_height_xyz_allE.png"), f, px_per_unit = 3.3)
f


##

using CairoMakie
CairoMakie.activate!()

#polar plot at different heights and energies
f = Figure(size = (1000, 700))
sleep(1)
axs = [PolarAxis(f[i, j],
    thetaticksvisible = false,
    alignmode=Mixed(bottom=0)
    ) for i in 1:4, j in 1:5]
hidedecorations!.(axs, grid=false)

# Add row labels (heights)
heights = [200, 150, 120, 105]
for (i2, h) in enumerate(heights)
    Label(f[i2, 0], "$(h) km", rotation = π/2, tellheight = false)
end

# Add row labels (heights)
energies = [0.5, 1, 2, 4, 8]
for (i1, e) in enumerate(energies)
    Label(f[0, i1], "$(e) keV", tellwidth = false)
end

for i in 1:4
    rowgap!(f.layout, i, Relative(-0.0))  # Minimal row gap
end
colgap!(f.layout, 1, Real(0.1))  # Minimal row gap
for j in 2:5
    colgap!(f.layout, j, Relative(-0.02))  # Minimal row gap
end
Label(f[-1, :], "Production Distribution in (r, θ) at Different Heights and Energies",)

p = nothing

sleep(1)
for (i2, r) in enumerate(runs_hrp_90)
    println(r)
    io = open(joinpath(dir, "hist_summed", r), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
    close(io)

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
            show_generators = false, strokewidth = 0,)
        rlims!(axs[i1, i2], 0.0, 20.0)
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
Colorbar(f[1:end, 6], p, label="Production [m⁻³]")
axs[4, 1].rticklabelsvisible = true
text!(axs[4, 1], 5.99, 12, text = "[m]", align = (:center, :top))
save(joinpath(dir, "plots", "panel_radial_hist_$(lim_pitch_deg)_all_E.png"), f, px_per_unit = 3.3)
display(f)




##

WGLMakie.activate!()

using CairoMakie
CairoMakie.activate!()


showlines = false
f_hr_ion = Figure(size = (900, 1200))
sleep(1)
axs_ion = [Axis(f_hr_ion[row, col], 
    #xlabel = "Radial Distance [m]", 
    #ylabel = "Height [km]", 
    limits = ((nothing, 20), (50, 600)),
    ytickformat = values -> ["" for value in values], 
    xtickformat = values -> ["" for value in values]
    ) for row in 1:5, col in 1:2] 
[ax.xtickformat = values -> ["$(value)" for value in values] for ax in axs_ion[end, :]]
[ax.ytickformat = values -> ["$(value)" for value in values] for ax in axs_ion[:, 1]]

f_hr_prod = Figure(size = (900, 1200))
sleep(1)
axs_prod = [Axis(f_hr_prod[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    limits = ((nothing, 20), (50, 600)),
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:5, col in 1:2] 
[ax.xtickformat = values -> ["$(value)" for value in values] for ax in axs_prod[end, :]]
[ax.ytickformat = values -> ["$(value)" for value in values] for ax in axs_prod[:, 1]]

f_prod_r_h = Figure(size = (900, 1200))
sleep(1)
axs_r_h = [Axis(f_prod_r_h[row, col], 
#    xlabel = "Radial Distance [m]", 
#    ylabel = "Height [km]", 
    limits = (nothing, (1e-4, 1e4)),
    yscale = log10,
    ytickformat = values -> ["" for value in values],
    xtickformat = values -> ["" for value in values]
    ) for row in 1:5, col in 1:2] 
[ax.xtickformat = values -> ["$(value)" for value in values] for ax in axs_r_h[end, :]]
[ax.ytickformat = Makie.automatic for ax in axs_r_h[:, 1]]

#+    limits = (nothing, (1e-4, 1e4)), 
#+    xlabel = "Radial Distance [m]", 
#+    ylabel = "Production [1/m³]",
 
heights = [105, 120, 140, 180, 250, 350, 500]

hm_ion = nothing
hm_prod = nothing

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
        dv = [dh*dr2*dp /2 for dh in diff(h_edges), dr2 in diff(r_edges.^2), dp in diff(p_edges)]
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

        lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = :dot, label = "Mean Gyroradius Pitch")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dash, label = "1/2e contour")

        if showlines
            lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
            contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
            contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
            lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
            lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
        end


        # Radial - Height plot normalised density
        data = dropdims(sum(his_hrp_norm.weights, dims = 3), dims = 3)' #./ dropdims(sum(area_pr, dims = 1), dims = 1)
        data_h_normal = data ./ maximum(data, dims = 1)
        data_h_normal[isnan.(data_h_normal)] .= 0.0
        i_max = [argmax(d) for d in eachcol(data)]

        ax = axs_prod[i1, i2]
        hm_prod = heatmap!(ax, r_edges, h_edges/1e3, data,
            colorscale = log10,
            colorrange = (1e-2, 1e4 #maximum(data)
            ), lowclip = ("white", 1e-2),
            )
        println( maximum(data))
        sleep(1)
        if lim_pitch_deg == 90
            text!(ax, 15, 460, text = "$E0 eV", font =:bold)
        end

        #lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = :dot, label = "Max Gyroradius Pitch")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(1)], color= "red", linestyle = :dash, label = "1/e contour")
        contour!(ax, r_middle, h_middle/1e3, data_h_normal, levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dash, label = "1/2e contour")

        if showlines
            lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
            contour!(ax, r_middle, h_middle[1:25]/1e3, data_h_normal[:, 1:25], levels=[maximum(data_h_normal)/exp(2)], color= "red", linestyle = :dashdot, label = "1/2e contour")
            lines!(ax, r_gyro_h_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
            lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
        end 

        ax = axs_r_h[i1, i2]
        for (i3, h) in enumerate(heights)  
            hi = findfirst(x -> x > h * 1e3, h_middle)
            d = dropdims(sum(his_hrp_norm.weights[hi, :, :], dims = 2), dims = 2)
            #lines!(ax, [-r_middle; r_middle], [d; d], label= "$h km", color = Makie.to_colormap(:batlow10)[i3+3])
            lines!(ax, [-r_middle; r_middle], [d; d], label= "$h km", color = Makie.to_colormap(:batlow)[i3*36])
            if lim_pitch_deg == 90
                text!(ax, 15, 460, text = "$E0 eV", font =:bold)
            end
        end
    end
end

Colorbar(f_hr_ion[2:4, 3], hm_ion, label="Counts [1]")
Colorbar(f_hr_prod[2:4, 3], hm_prod, label="Production [m⁻³]")

axs_ion[1, 1].title = "20 deg"
axs_ion[1, 2].title = "isotropic"

axs_prod[1, 1].title = "20 deg"
axs_prod[1, 2].title = "isotropic"

axs_r_h[1, 1].title = "20 deg"
axs_r_h[1, 2].title = "isotropic"

[ax.xlabel = "Radial Distance [m]" for ax in axs_ion[end, :]]
[ax.ylabel = "Height [km]" for ax in axs_ion[:, 1]]

[ax.xlabel = "Radial Distance [m]" for ax in axs_prod[end, :]]
[ax.ylabel = "Height [km]" for ax in axs_prod[:, 1]]

[ax.xlabel = "Radial Distance [m]" for ax in axs_r_h[end, :]]
[ax.ylabel = "Production [m⁻³]" for ax in axs_r_h[:, 1]]

axislegend(axs_ion[1, 1])
#axislegend(axs_prod[1, 1])
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

#save(joinpath(dir, "plots", "panel_hist_hr_all_E_pitch.png"), f_hr_ion, px_per_unit = 3.3)
#save(joinpath(dir, "plots", "panel_hist_norm_hr_all_E_pitch.png"), f_hr_prod, px_per_unit = 3.3)
save(joinpath(dir, "plots", "panel_prod_r_h_all_E_pitch.png"), f_prod_r_h, px_per_unit = 3.3)


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
    l1 = lines!(ax, p_middle, r_middle[i_max],
        color= "red", label = "Max Production")
    l2 = contour!(ax, p_middle, r_middle, data, 
        levels=[maximum(data)/exp(1)], color= "red", 
        linestyle = :dash, label = "1/e contour")
    l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), 
        color= "black", linestyle = :dot, label = "Mean Gyroradius")
    Colorbar(fig[2, 2], p, label="Production [m⁻³]")
    Legend(fig[1, 2], [l1, l2, l4], ["Max Production", "1/e contour", "Mean Gyroradius"])
    ax.titlegap = 20
    #fig
    save(joinpath(dir, "plots", "hist3d_phase_rd_norm$(E0)_$(lim_pitch_deg)_$(h_middle[h_ind]/1e3)km.png"), fig)
    display(fig)
end



for h_ind in [18]#, 5, 8, 13, 23]
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
