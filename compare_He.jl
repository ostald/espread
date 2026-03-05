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

dir = "results/r4_conicB_2025-09-05T14:19:27.566/"
dir_he = "results/r8_conicB_He_500eV_2026-01-21T19:21:08.258/"
dir_he = "results/r10_scatter45deg_2026-01-26T17:12:40.089/"


#compare production profiles

fig = Figure()
sleep(1)
ax = Axis(fig[1, 1],
    xscale = log10)

for d in [dir, dir_he]
    io = open(joinpath(d, "hist_summed", "h_xyz_4000.0eV_20.0deg_summed.hist"), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_xyz = deserialize(io)
    close(io)
    println(his_xyz.isdensity)

    new_edges = (filter(x -> mod(x, 1) == 0.5, his_xyz.edges[1]), 
        filter(x -> mod(x, 1) == 0.5, his_xyz.edges[2]), 
        filter(x -> mod(x, 1000) == 0, his_xyz.edges[3])
        )

    his_xyz = rebin(his_xyz, new_edges);
    #his_xyz = normalize(his_xyz, mode=:density)
    # normalise by density is the same as dividing by bin volume
    # check:
    #x_edges = his_xyz.edges[1]
    #y_edges = his_xyz.edges[2]
    #z_edges = his_xyz.edges[3]
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true

    data = dropdims(sum(his_xyz.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_xyz.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   

    lines!(ax, data, z_middle/1e3, label = "$E0 eV")
end
xlims!(1, 1e7)
display(fig)




CairoMakie.activate!()
fig = Figure()
sleep(1)
ax = Axis(fig[1, 1],
    limits = (sqrt.((4e2, 1e4)./1e3), (-0.3, 13.3)),
    xticks = sqrt.([5e2, 1e3, 2e3, 4e3, 8e3]./1e3),
    #ax.title = L"\mathbf{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}"
    xtickformat = values -> [rich("$(Int(value))", superscript("1/2")) for value in [5e2, 1e3, 2e3, 4e3, 8e3]], 
    #xscale = log10,
    xlabel = rich("(Energy [keV])", superscript("1/2")),
    yticks = [0, 1, 2, 3, 4, 5, 6, 7].*2,
    ylabel = "Width [m]"
    )
width_20 = zeros(5)
width_90 = zeros(5)
#for (i2, collection) in enumerate([runs_hrp_20, runs_hrp_90])
#    for (i1, r) in enumerate(collection)
        #io = open(joinpath(dir, "hist_summed", "h_hrp_4000.0eV_20.0deg_summed.hist"), "r")
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
        if lim_pitch_deg == 20
            #scatter!(ax, sqrt(E0/1e3), r_middle[i_min], color = Makie.wong_colors()[1])
            width_20[i1] =  r_middle[i_min]
        else
            #scatter!(ax, sqrt(E0/1e3), r_middle[i_min], color = Makie.wong_colors()[2])
            width_90[i1] =  r_middle[i_min]
        end
    end
end
scatterlines!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3]./1e3), 2 .* width_20, color = Makie.wong_colors()[1],label = L"{\mathrm{field-aligned \, \theta_{lim} = 20^{\degree}}}")
scatterlines!(ax, sqrt.([5e2, 1e3, 2e3, 4e3, 8e3]./1e3), 2 .* width_90, color = Makie.wong_colors()[2],label = "isotropic")
axislegend(ax, position = :lt)
display(fig)
save(joinpath(dir, "plots", "width_h_int.png"), fig, px_per_unit = 3.3)

#for d in [dir, dir_he]
    io = open(joinpath(d, "hist_summed", "h_hrp_4000.0eV_90.0deg_summed.hist"), "r")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
    close(io)
    #println(his_hrp.isdensity)

    new_edges = (filter(x -> mod(x, 1) == 0.5, his_hrp.edges[1]), 
        filter(x -> mod(x, 1) == 0.5, his_hrp.edges[2]), 
        filter(x -> mod(x, 1000) == 0, his_hrp.edges[3])
        )

    his_xyz = rebin(his_xyz, new_edges);
    #his_xyz = normalize(his_xyz, mode=:density)
    # normalise by density is the same as dividing by bin volume
    # check:
    #x_edges = his_xyz.edges[1]
    #y_edges = his_xyz.edges[2]
    #z_edges = his_xyz.edges[3]
    #dv = [dx*dy*dz for dx in diff(x_edges), dy in diff(y_edges), dz in diff(z_edges)]
    #his_xyz.weights ./ dv == (normalize(his_xyz, mode=:density)).weights
    # >>> true

    data = dropdims(sum(his_hrp.weights, dims = (1, 2)), dims = (1, 2))
    z_edges = his_hrp.edges[3]
    z_middle = z_edges[1:end-1] + diff(z_edges)/2   

end
