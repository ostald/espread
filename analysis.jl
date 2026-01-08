using StatsBase
#using CSV
using DataFrames
using LinearAlgebra
using JLD2
using Serialization
using WGLMakie
using LsqFit
include("constants.jl")
include("magnetic_field.jl")

#dir = "results/run1_2025-07-19T21:53:59.887/"
dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
if !isdir(joinpath(dir, "plots"))
    mkdir(joinpath(dir, "plots"))
end

dir_con = readdir(dir)
dir_con_raw = filter(x-> contains(x, ".bin"), dir_con)

runs = unique([d[1:end-8] for d in dir_con_raw])


#for r in runs
    #filter_crit = r
    filter_crit = "res_8000.0eV_20.0deg"

    files = filter(x-> contains(x, filter_crit), dir_con_raw)
    files = files[1:10]

#    res = Vector{Any}(undef, length(files))
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = [0, 0, 0, 0, 0, 0]

    df = DataFrame(generation = Int[],
            idx_scatter = Int[], 
            r0 = Vector{Float64}[], 
            v0 = Vector{Float64}[], 
            status = Int[], 
            r = Vector{Float64}[], 
            v = Vector{Float64}[])

    #f = files[2]
    @time for (id, f) in enumerate(files)
        io = open(joinpath(dir,  f), "r")
        E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
        data = [deserialize(io) for _ in 1:Int(1e6) if !eof(io)]

        for d in data
            push!(df, d)
        end
        #df_arr[id] = df

        #while !eof(io)
        #    push!(df, deserialize(io))
        #end
        close(io)


        #data = load(dir * f)
        #all_keys = collect(keys(data))
        #df_keys = filter(x-> !contains(x, "setup"), all_keys)

        #df = vcat([data[k] for k in sort(df_keys)]...)
        #unique(df)

        """
        df = CSV.read(open(joinpath(dir, f)),
            header=false,
            delim = "\t",
            skipto=8,
            DataFrame,
            )
        rename!(df, [:r0, :v0, :status, :r, :v])
        df.r0 = [eval(Meta.parse(value)) for value in df.r0]
        df.v0 = [eval(Meta.parse(value)) for value in df.v0]
        df.r  = [eval(Meta.parse(value)) for value in df.r]
        df.v  = [eval(Meta.parse(value)) for value in df.v]
        """
        #res[id] = df
    end
    #df = vcat(df_arr...)
    #df_comb = unique(vcat(res...))
    
    #@time jldsave(joinpath(dir, r * ".jld2"); df, E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, compress = false)

#end
#exit()

#df_dict = load(joinpath(dir, "res_1000.0eV_20.0deg.jld2"))
#df = unique(df_dict["df_comb"])

"""
res_file = "results/conicB_run2_2025-08-26T12:50:35.891/res_8000.0eV_20.0deg_001.bin"
df = DataFrame(Generation = Int[],
    idx_scatter = Int[], 
    r0 = Vector{Float64}[], 
    v0 = Vector{Float64}[], 
    status = Int[], 
    r = Vector{Float64}[], 
    v = Vector{Float64}[])
io = open(res_file, "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
while !eof(io)
    push!(df, deserialize(io))
end
close(io)
"""

## sanity checks
# check for fails in Boris mover:
if size(filter(row -> :status .== 0, df))[1] > 0
    error("Failure of Boris mover!")
end

if any(nonunique(df))
    error("Dublicates found.")
end


## 
df.E0 = E_ev.(norm.(df.v0))
df.E_end = E_ev.(norm.(df.v))

#check for status -1 and Energy larger than lowest ionization energy
filter(row -> row.E0 < 12.072 && row.status != -1, df)


## primary electrons:
df.alt0 = altitude.(df.r0)
df.alt_end = altitude.(df.r)
df_alt0 = filter(:alt0 => x -> x > 599e3, df)

#injection points of primaries:
fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax, Point3.(df_alt0.r0), markersize = 2)#, 
#display(fig)
save(joinpath(dir, "plots", "r0_primary_e_$(E0)_$(lim_pitch_deg).png"), fig)


##
# select endpoint for primary electrons, starting point for secondary electrons
df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)
df.pos_ce = ifelse.(df.alt0 .> 599e3, df.r, df.r0)
df.pos = [p - [0, 0, c.re] for p in df.pos_ce]

hmsis = hmin:hintervals:hmax

using WGLMakie
WGLMakie.activate!()

#Altitude histogram
#using GLMakie
fig, ax, his = hist(df.alt./1e3,
    bins = hmsis./1e3, 
    direction=:x,
    axis = (xlabel = "Eelctron Count [1]",
        ylabel = "Height [km]",
        limits = ((1, nothing), nothing),
        xscale = log10,
        ),
    )
display(fig)
save(joinpath(dir, "plots", "hist_alt_$(E0)_$(lim_pitch_deg).png"), fig)



# 3d point plot of endpoints (earth centered coordinates)
fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1),)
sc = scatter!(ax,Point3.(df.pos), markersize = 2)#,
    #axis=(limits=(nothing, nothing, nothing),),)
#zlims!(ax, 5.99e6, 6.01e6)
#ylims!(ax, 2.465e6, 2.48e6)
#xlims!(ax, -1.5e4/2, 1.5e4/2)
display(fig)
save(joinpath(dir, "plots", "position_$(E0)_$(lim_pitch_deg).png"), fig)


# zoom in by using smaller dataset
filter_hmin = 150e3
filter_hmax = 161e3
df_filtered = filter(row -> row.alt > filter_hmin && row.alt < filter_hmax, df)

# meshscatter keeps dimensions const
fig, ax, ms = scatter(Point3.(df_filtered.pos), markersize = 1e1)
display(fig)


fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax,Point3.(df_filtered.pos), markersize = 1e1)
scatter!(ax,Point3(sum(df_filtered.pos)./size(df_filtered.pos, 1)), markersize = 1e1)
display(fig)
save(joinpath(dir, "plots", "position_$(filter_hmin)-$(filter_hmax)km_$(E0)_$(lim_pitch_deg).png"), fig)





# investigate phase angle
# first make sure the magntic field line tracing from the gyrocenter is centered at [0, 0]
p0 = [0, 0, 600e3 + c.re]
p_along_B = zeros(600_000, 3)

p = p0
for i in axes(p_along_B, 1)
    p_along_B[i, :] = p
    B = convergent_vertical_field(p)
    p = p + B/norm(B)
end

if !iszero(p_along_B[:, 1:2][:])
    error("path along fieldine is not centered to [0, 0]")
end

function phase_angle(p)
    return atan(p[1], p[2])
end

df.phase = phase_angle.(df.pos)

fig, ax, sc = scatter([Point2(p[1:2]) for p in df.pos], markersize = 2, axis=(aspect = 1,),)
save(joinpath(dir, "plots", "horizontal_scatter_$(E0)_$(lim_pitch_deg).png"), fig)


function PolarHist(anglesRad, limits, title, colour, transparancy, graphposition)
    binsRad = collect(-pi:pi/18:pi)

    h = fit(Histogram, anglesRad, binsRad)
    x = collect(h.edges[1])
    y = convert.(Float64, h.weights)

    f = Figure()
    ax = PolarAxis(f[1, graphposition],
        title="$title",
        thetaticks=((collect(0:15:limits)) .* 2*pi ./ limits, string.(string.(collect(0:15:limits)), "°")),
        #rticks=collect(0:10:ceil(maximum(y) / 10)*10)
    )
    thetalims!(ax, 0, 2*pi )
    #rlims!(ax, 0, ceil(maximum(y) / 10) * 10)
    for i in eachindex(y)
        if y[i] > 0
            Makie.poly!(ax, [(0.0, 0.0), (x[i], y[i]), (x[i+1], y[i])], [1, 2, 3], strokewidth=1.5, strokecolor=:black, color=Makie.wong_colors()[colour], alpha=transparancy)
        end
    end
    return f
end
fig = PolarHist(df.phase, 360, "$E0 eV, $lim_pitch_deg deg", 1, 0.6, 1)
fig
save(joinpath(dir, "plots", "hist_phase_$(E0)_$(lim_pitch_deg).png"), fig)


df.rdist = [norm(p[1:2]) for p in df.pos]

maximum(df.alt - [p[3] for p in df.pos])

minimum(df.alt)

##

h_bins = collect(80e3:10e3:600e3)
r_bins = collect(0:0.1:40)
p_bins = collect(-pi:pi/18:pi)

h_middle = h_bins[1:end-1] + diff(h_bins)/2
r_middle = r_bins[1:end-1] + diff(r_bins)/2
p_middle = p_bins[1:end-1] + diff(p_bins)/2

B = norm.([convergent_vertical_field([0, 0, h+c.re]) for h in h_middle])
r_gyro_max = v_abs(E0) * c.me ./ (c.qe * B)
lim_pitch = lim_pitch_deg /180 *pi

mean_vperp_theory = (1/2 * lim_pitch - 1/4*sin(2*lim_pitch)) / (1-cos(lim_pitch))

B_top = convergent_vertical_field([0, 0, 600e3+c.re])

mean_vperp_data = mean(norm.(df_alt0.v0 .- (dot.(df_alt0.v0, [B_top]) .* [B_top] ./ norm(B_top)^2)))/v_abs(E0)
mean_vperp_data2 = mean(sin.(atan.(norm.(cross.(df_alt0.v0,[B_top])),dot.(df_alt0.v0,[B_top]))))

if mean_vperp_data-mean_vperp_data2 > 1e-5
    error("invistigate")
end


if mean_vperp_data-mean_vperp_theory> 1e-2
    error("invistigate")
end

v_perp_mean = mean_vperp_theory * v_abs(E0)
r_gyro_mean = v_perp_mean * c.me ./ (c.qe * B)


# Radial Distance - Height Histogram
his_rh = fit(Histogram, (df.rdist, df.alt), (r_bins, h_bins))
i_max_r_h = [argmax(h) for h in eachcol(his_rh.weights)]
fig, ax, hm = heatmap(r_bins, h_bins/1e3, his_rh.weights,
    colorscale = log10,
    colorrange = (1, maximum(his_rh.weights)), lowclip = ("white", 0),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = (nothing, (50, 600)),)
    )
lines!(ax, r_middle[i_max_r_h], h_middle/1e3, color= "red", label = "Max Count")
lines!(ax, r_gyro_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Counts [1]")
axislegend(ax)
fig
save(joinpath(dir, "plots", "hist2d_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)

fig, ax, sc = scatter(r_middle, 
    his_rh.weights[:, 3], 
    axis=(xlabel = "Radial Distance [m]", 
        ylabel = "Counts [1]", 
        title = "$E0 eV, $lim_pitch_deg deg, $(h_middle[3]/1e3)km",
        ),
    )


# Radial Distance - Height Histogram, normalized to areal density of radial bins
his_rh = fit(Histogram, (df.rdist, df.alt), (r_bins, h_bins))
area_r = (r_bins[2:end] .^2 - r_bins[1:end-1] .^2) * pi
density_rh = his_rh.weights ./ area_r

his_rh = fit(Histogram, (df.rdist, df.alt), weights(1 ./df.rdist), (r_bins, h_bins))
density_rh2 = his_rh.weights

i_max_r_h = [argmax(h) for h in eachcol(density_rh)]
fig, ax, hm = heatmap(r_bins, h_bins/1e3, density_rh2,
    colorscale = log10,
    colorrange = (1e-2, maximum(density_rh)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = (nothing, (50, 600)),)
    )
lines!(ax, r_middle[i_max_r_h], h_middle/1e3, color= "red", label = "Max Count")
lines!(ax, r_gyro_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Density [1/m²]")
axislegend(ax)
fig
save(joinpath(dir, "plots", "hist2d_areaDensity_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)

fig, ax, sc = scatter(r_middle, 
    density_rh[:, 3], 
    axis=(xlabel = "Radial Distance [m]", 
        ylabel = "Density [1/m²]", 
        title = "$E0 eV, $lim_pitch_deg deg, $(h_middle[3]/1e3)km",
        ),
    )
scatter!(ax, r_middle, density_rh2[:, 3], color = "red")
fig

fig, ax, sc = scatter(r_middle, 
    density_rh[:, 3] ./ density_rh2[:, 3], )


# Phase - Heigth Histogram
his_ph = fit(Histogram, (df.phase, df.alt), (p_bins, h_bins))
fig, ax, hm = heatmap(p_bins, h_bins/1e3, his_ph.weights,
    colorscale = log10,
    colorrange = (1, maximum(his_ph.weights)), lowclip = ("white", 0),
    axis = (xlabel = "Phase [rad]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = (nothing, (50, 600)),)
    )
Colorbar(fig[1, 2], hm, label="Counts [1]")
fig
save(joinpath(dir, "plots", "hist2d_phase_height_$(E0)_$(lim_pitch_deg).png"), fig)


# Phase - Radial Distance Histogram, Polar plot
his_pr = fit(Histogram, (df.phase, df.rdist), (p_bins, r_bins))
i_max_r_p = [argmax(h) for h in eachrow(his_pr.weights)]

fig = Figure()
ax = PolarAxis(fig[1, 2], title = "$E0 eV, $lim_pitch_deg deg",)
p = voronoiplot!(ax, p_middle, r_middle, his_pr.weights,
    colorscale = log10,
    colorrange = (1, maximum(his_pr.weights)), lowclip = ("white", 0),
    show_generators = false, strokewidth = 0)
#rlims!(ax, 0.0, 10.5)
l1 = lines!(ax, p_middle, r_middle[i_max_r_p], color= "red", label = "Max Count")
Colorbar(fig[1, 3], p, label="Counts [1]")
Legend(fig[1, 1], [l1], ["Max Count"])
fig
save(joinpath(dir, "plots", "hist2d_phase_rd_$(E0)_$(lim_pitch_deg).png"), fig)

# Phase - Radial Distance Histogram, Cartesian
fig, ax, hm = heatmap(p_bins, r_bins, his_pr.weights,
    colorscale = log10,
    colorrange = (1, maximum(his_pr.weights)), lowclip = ("white", 0),
    axis = (xlabel = "Phase [rad]", ylabel = "Radial Distance [m]", title = "$E0 eV, $lim_pitch_deg deg",),)
Colorbar(fig[1, 2], hm, label="Counts [1]")
fig
#save(joinpath(dir, "plots", "hist2d_v2_phase_rd_$(E0)_$(lim_pitch_deg).png"), fig)

# Phase - Radial Distance Histogram, Polar plot, normalized to areal density 
area_pr = [dr * dp / 2  for dp in diff(p_bins), dr in (r_bins[2:end] .^2 - r_bins[1:end-1] .^2)]
density_pr = his_pr.weights ./ area_pr
i_max_r_p = [argmax(h) for h in eachrow(density_pr)]

fig = Figure()
ax = PolarAxis(fig[1, 2], title = "$E0 eV, $lim_pitch_deg deg",)
p = voronoiplot!(ax, p_middle, r_middle, density_pr,
    colorscale = log10,
    colorrange = (1e-2, maximum(density_pr)), lowclip = ("white", 1e-2),
    show_generators = false, strokewidth = 0)
#rlims!(ax, 0.0, 10.5)
l1 = lines!(ax, p_middle, r_middle[i_max_r_p], color= "red", label = "Max Count")
Colorbar(fig[1, 3], p, label="Density [1/m²]")
Legend(fig[1, 1], [l1], ["Max Count"])
fig
save(joinpath(dir, "plots", "hist2d_areaDensity_phase_rd_$(E0)_$(lim_pitch_deg).png"), fig)


#also do densities (divide by circumference) done
#scattering depth instead of height? => no
# countour of 1/e decay from max


#3D histogram in height, phase and radius

his_hrp = fit(Histogram, (df.alt, df.rdist, df.phase), (h_bins, r_bins, p_bins))


# use slices to visualize, or sum along phase
## Radial Phase plot
h_ind = 4
data = his_hrp.weights[h_ind, :, :]' 
i_max = [argmax(d) for d in eachrow(data)]
i_exp_hs = [findmin(abs.(d[i:end].-d[i]/exp(1)))[2] + i - 1 for (d, i) in zip(eachrow(data), i_max)]
i_exp_ls = [findmin(abs.(d[1:i].-d[i]/exp(1)))[2] for (d, i) in zip(eachrow(data), i_max)]


fig = Figure()
ax = PolarAxis(fig[1:2, 1], title = "$E0 eV, $lim_pitch_deg deg, $(h_bins[4]/1e3)km",)
p = voronoiplot!(ax, p_middle, r_middle, data,
    colorscale = log10,
    colorrange = (1, maximum(data)), lowclip = ("white", 0),
    show_generators = false, strokewidth = 0)
rlims!(ax, 0.0, 20.0)
l1 = lines!(ax, p_middle, r_middle[i_max], color= "red", label = "Max Counts")
l2 = lines!(ax, p_middle, r_middle[i_exp_hs], color= "red", linestyle = :dash, label = "1/e contour")
l3 = lines!(ax, p_middle, r_middle[i_exp_ls], color= "red", linestyle = :dash, label = "1/e contour")
l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), color= "black", linestyle = :dot, label = "Mean Gyroradius")
Colorbar(fig[2, 2], p, label="Counts [1]")
Legend(fig[1, 2], [l1, l2, l4], ["Max Counts", "1/e contour", "Mean Gyroradius"])
fig
save(joinpath(dir, "plots", "hist3d_phase_rd_$(E0)_$(lim_pitch_deg)_$(h_bins[4]/1e3)km.png"), fig)


## Radial Phase plot normalised density
h_ind = 4
data = his_hrp.weights[h_ind, :, :]' ./ area_pr
i_max = [argmax(d) for d in eachrow(data)]
i_exp_hs = [findmin(abs.(d[i:end].-d[i]/exp(1)))[2] + i - 1 for (d, i) in zip(eachrow(data), i_max)]
i_exp_ls = [findmin(abs.(d[1:i].-d[i]/exp(1)))[2] for (d, i) in zip(eachrow(data), i_max)]


fig = Figure()
ax = PolarAxis(fig[1:2, 1], title = "$E0 eV, $lim_pitch_deg deg, $(h_bins[4]/1e3)km",)
p = voronoiplot!(ax, p_middle, r_middle, data,
    colorscale = log10,
    colorrange = (1e-1, maximum(data)), lowclip = ("white", 1e-1),
    show_generators = false, strokewidth = 0)
#rlims!(ax, 0.0, 20.0)
l1 = lines!(ax, p_middle, r_middle[i_max], color= "red", label = "Max Density")
l2 = lines!(ax, p_middle, r_middle[i_exp_hs], color= "red", linestyle = :dash, label = "Max Count")
l3 = lines!(ax, p_middle, r_middle[i_exp_ls], color= "red", linestyle = :dash, label = "Max Count")
l4 = lines!(ax, p_middle, fill(r_gyro_mean[h_ind], size(p_middle)), color= "black", linestyle = :dot, label = "Mean Gyroradius")
Colorbar(fig[2, 2], p, label="Density [1/m²]")
Legend(fig[1, 2], [l1, l2, l4], ["Max Density", "1/e contour", "Mean Gyroradius"])
fig
save(joinpath(dir, "plots", "hist3d_density_phase_rd_$(E0)_$(lim_pitch_deg)_$(h_bins[4]/1e3)km.png"), fig)

"""
function find_ind_dist_width(data)
    i_max = [argmax(d) for d in eachrow(data)]
    i_hs = zeros(size(i_max))
    i_hs = zeros(size(i_max))
    for (d, i1) in zip(eachrow(data), i_max)
        value, i2 = findmin(abs.(d[i1:end].-d[i1]/exp(1)))
        i_hsi1 + i2
    end
"""

## Radial - Height plot

data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)
i_max = [argmax(d) for d in eachrow(data)]
i_exp_hs = [findmin(abs.(d[i:end].-d[i]/exp(1)))[2] + i - 1 for (d, i) in zip(eachrow(data), i_max)]
i_exp_ls = [findmin(abs.(d[1:i].-d[i]/exp(1)))[2] for (d, i) in zip(eachrow(data), i_max)]

fig, ax, hm = heatmap(r_bins, h_bins/1e3, data',
    colorscale = log10,
    colorrange = (1, maximum(data)), lowclip = ("white", 0),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 20), (50, 600)),)
    )
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Count")
lines!(ax, r_middle[i_exp_hs], h_middle/1e3, color= "red", linestyle = :dash, label = "1/e contour")
lines!(ax, r_middle[i_exp_ls], h_middle/1e3, color= "red", linestyle = :dash)#, label = "1/e contour")
lines!(ax, r_gyro_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Counts [1]")
axislegend(ax)
fig
save(joinpath(dir, "plots", "hist2d_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)




## Radial - Height plot normalised density
area_r = (r_bins[2:end] .^2 - r_bins[1:end-1] .^2) * pi

data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' ./ area_r

i_max = [argmax(d) for d in eachcol(data)]
i_exp_hs = [findmin(abs.(d[i:end].-d[i]/exp(1)))[2] + i - 1 for (d, i) in zip(eachcol(data), i_max)]
i_exp_ls = [findmin(abs.(d[1:i].-d[i]/exp(1)))[2] for (d, i) in zip(eachcol(data), i_max)]

fig, ax, hm = heatmap(r_bins, h_bins/1e3, data,
    colorscale = log10,
    colorrange = (1e-2, maximum(data)), lowclip = ("white", 1e-2),
    axis = (xlabel = "Radial Distance [m]", ylabel = "Height [km]", title = "$E0 eV, $lim_pitch_deg deg",
        limits = ((nothing, 20), (50, 600)),)
    )
lines!(ax, r_middle[i_max], h_middle/1e3, color= "red", label = "Max Density")
lines!(ax, r_middle[i_exp_hs], h_middle/1e3, color= "red", linestyle = :dash, label = "1/e contour")
lines!(ax, r_middle[i_exp_ls], h_middle/1e3, color= "red", linestyle = :dash)#, label = "1/e contour")
lines!(ax, r_gyro_max , h_middle/1e3, color= "black", label = "Max Gyroradius")
lines!(ax, r_gyro_mean, h_middle/1e3, color= "black", linestyle = :dash, label = "Mean Gyroradius")
Colorbar(fig[1, 2], hm, label="Density [1/m²]")
axislegend(ax)
fig
save(joinpath(dir, "plots", "hist2d_rd_height_$(E0)_$(lim_pitch_deg).png"), fig)




# wideining at different heights
# is the radial distribution a gaussian? chekc 
# => ceck 1d (for statistics) but also 2D (use x, y coord instaed of polar)
# vertical scattering depthof primaries?
# comparison between runs
# animation going through height of phase, read 

using CairoMakie
CairoMakie.activate!()

function gaussian(x, p0)
    A, std= p0
    gamma = 2
    return A * exp.( -(0.5 .* ((x)./std).^2) .^(gamma/2))
end

p0 = [100.0, 5.0]

function sum_gaussian(x, p0)
    A1, std1, A2, std2 = p0
    gamma = 2
    return A1 * exp.( -(0.5 .* ((x)./std1).^2) .^(gamma/2)) .+ A2 * exp.( -(0.5 .* ((x)./std2).^2) .^(gamma/2))
end



for i in axes(h_middle, 1)
    d = data[:, i]
    if sum(d) < 10
        continue
    end

    fig = Figure()
    ax = Axis(fig[1, 1], 
        title = "$E0 eV, $lim_pitch_deg deg, $(h_middle[i]/1e3)km", 
        limits = (nothing, (1e-4, nothing)), 
        xlabel = "Radial Distance [m]", 
        ylabel = "Density [1/m²]",
        yscale = log10,
        )

    scatter!(ax, [-r_middle; r_middle], [d; d], label= "data")
    #display(fig)

    ff = curve_fit(gaussian, r_middle, d, p0)
    #lines!(ax, -40:0.1:40, (gaussian(-40:0.1:40, ff.param)), label = "fit")

    lb = [0.0, 0, 0, 0]
    ff2 = curve_fit(sum_gaussian, r_middle, d, [ff.param ; ff.param[1]/10; ff.param[2]*2]; lower=lb)
    lines!(ax, -40:0.1:40, (sum_gaussian(-40:0.1:40, ff2.param)), label = "fit")
    axislegend(ax)
    #xlims!(ax, -20, 20)
    #save(joinpath(dir, "plots", "hist2d_rd_height_gaussfit_$(E0)_$(lim_pitch_deg)_$(h_middle[i]/1e3)km.png"), fig)
    display(fig)

    dv = [dh*dr*dp /2 for dh in diff(h_bins), dr in diff(r_bins.^2), dp in diff(p_bins)]
    data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)' ./ area_r
    d = data[:, i]

    boundary = (-20, 20)
    npoints = 40
    #kernel = :normal
    bandwidth = 1.0
    #U = kde(d, boundary=boundary, npoints=npoints, bandwidth=bandwidth)
    #lines!(ax, U.x, U.density, color = "green", label = "KDE")

end



z_bins = collect(80e3:10e3:600e3)
x_bins = collect(-20.5:1:20.5)
y_bins = copy(x_bins)

z_middle = z_bins[1:end-1] + diff(z_bins)/2
x_middle = x_bins[1:end-1] + diff(x_bins)/2
y_middle = y_bins[1:end-1] + diff(y_bins)/2

his_xyz = fit(Histogram, ([p[1] for p in df.pos], [p[2] for p in df.pos], [p[3] for p in df.pos]), (x_bins, y_bins, z_bins))

function twoD_Gaussian(xy, p)
    amplitude, xo, yo, sigma_x, sigma_y, theta = p
    offset = 0
    a = (cos(theta)^2)/(2*sigma_x^2) + (sin(theta)^2)/(2*sigma_y^2)
    b = -(sin(2*theta))/(4*sigma_x^2) + (sin(2*theta))/(4*sigma_y^2)
    c = (sin(theta)^2)/(2*sigma_x^2) + (cos(theta)^2)/(2*sigma_y^2)

    # creating linear meshgrid from xy
    x = xy[:, 1]
    y = xy[:, 2]
    g = offset .+ amplitude .* exp.( - (a.*((x .- xo).^2) + 2 .* b .* (x .- xo) .* (y .- yo) + c * ((y .- yo).^2)))
    return g[:]
end

function gaussian2D(xy, p0)
    x = xy[:, 1]
    y = xy[:, 2]
    A, stdx, stdy = p0
    return A .* exp.(.-0.5 .* ((x)./stdx).^2 .-0.5 .* ((y)./stdy).^2) 
end

p0 = [100.0, 5.0, 5.0, 100.0]
p = Float64.([3000, 0, 0, 2, 2, 0])

xy_middle =  hcat([[x, y] for x in x_middle, y in y_middle][:]...)'

for i in axes(z_middle, 1)
    d = his_xyz.weights[:, :, i]
    if sum(d) < 10
        continue
    end

    fig, ax, hm = heatmap(x_bins, y_bins, d,
        colorscale = log10,
        colorrange = (1, maximum(d)), lowclip = ("white", 0),
        axis = (xlabel = "x [m]", ylabel = "y [m]", title = "Data",
            limits = ((-20, 20), (-20, 20)),
            aspect = 1,
            )
        )
    #Colorbar(fig[1, 2], hm, label="Counts [1]")
    #display(fig)
    
    ff = curve_fit(twoD_Gaussian, xy_middle, d[:], p)
    fit_d = twoD_Gaussian(xy_middle, ff.param) |> (res -> reshape(res, length(x_middle), length(y_middle)))
    ax, hm = heatmap(fig[1, 2], x_bins, y_bins, fit_d,
        colorscale = log10,
        colorrange = (1, maximum(d)), lowclip = ("white", 0),
        axis = (xlabel = "x [m]", ylabel = "y [m]", title = "Fit",
            limits = ((-20, 20), (-20, 20)),
            aspect = 1,
            )
        )
    #display(fig)


    ax, hm = heatmap(fig[2, 1],x_bins, y_bins, abs.(d-fit_d),
        colorscale = log10,
        colorrange = (1, maximum(d)), lowclip = ("white", 0),
        axis = (xlabel = "x [m]", ylabel = "y [m]", title = "Residuals",
            limits = ((-20, 20), (-20, 20)),
            aspect = 1,
            )
        )

    Colorbar(fig[1:2, 3], hm, label="Counts [1]") 

    ax, lin = lines(fig[2, 2], x_middle, fit_d[:, 20], label = "fit")
    scatter!(ax, x_middle, d[:, 21], label = "fit")
    supertitle = Label(fig[0, :], "$E0 eV, $lim_pitch_deg deg, $(h_middle[i]/1e3)km")
    display(fig)
end

using Bonito
Bonito.set_cleanup_time!(1)
