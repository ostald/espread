

#loading the results:
using Serialization
dir = "results/r19_Yuri_pulsatingAurora_2026-08-25T16:10:33.477/" # directory to hist_summed folder
io = open(joinpath(dir, "hist_summed", "h_hrp_100000.0eV_90.0deg_summed.hist"), "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, his_hrp = deserialize(io)
close(io)
# E0 is the primary electron Energy
# lim_pitch_deg is the pitch angel distribution
# seed_value is the seed used for initiating the random number generator
# hmin is the lower altitude limit of the model atmosphere
# hmax is the upper altitude limit of the model atmosphere
# hinterval is the resolution of the model atmosphere
# his_hrp is the histogram in cylindrical coordinates
#   his_hrp.edges contains the edges of the histgram bins
#       with the first dimension being height
#       second dimension is radial Distance
#       third dimension is phase angle
#   his.weights contains the counts in each bin, with the same ordering as above.


using LinearAlgebra
using ImageFiltering
using StatsBase
include("analysis_util.jl")
#rebinning: 10 altitude bins together
his_hrp = rebin(his_hrp, (
    his_hrp.edges[1][1:10:end],
    his_hrp.edges[2],
    his_hrp.edges[3],
    ))


# make radial production plot 
function hist_radial_heihgt_production(ax, his_hrp)
    #extract edges
    h_edges = his_hrp.edges[1]
    r_edges = his_hrp.edges[2]
    p_edges = his_hrp.edges[3]

    #calculate middle points
    h_middle = h_edges[1:end-1] + diff(h_edges)/2
    r_middle = r_edges[1:end-1] + diff(r_edges)/2
    p_middle = p_edges[1:end-1] + diff(p_edges)/2

    # Radial - Height plot normalised density
    #sum over phase bins
    data = dropdims(sum(his_hrp.weights, dims = 3), dims = 3)'
    #calculate bin size
    dA = [dh*dr2*pi for dh in diff(h_edges), dr2 in diff(r_edges.^2)]'
    #normalize histogram
    data = data ./ dA

    #normalize with respect to the maximum at every height slice for contours
    data_h_normal = data ./ maximum(data, dims = 1)
    data_h_normal[isnan.(data_h_normal)] .= 0.0
    i_max = [argmax(d) for d in eachcol(data)]

    #plot
    hm_prod = heatmap!(ax, r_edges, h_edges/1e3, data,
        colorscale = log10,
        colorrange = (1e-4, 1e2 #maximum(data)
        ), lowclip = ("white", 1e-4),
        )
        
    text!(ax, 15, 460, text = "$E0 eV", font =:bold)
    
    data_h_normal_f = mapwindow(median!, data_h_normal, (3, 3))
    #lines!(ax, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = (:dot, 10), linewidth = 1, label = "Max gyroradius")
    contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(1)], color= "white", linestyle = (:dash, 10), label = "1/e contour")
    contour!(ax, r_middle, h_middle/1e3, data_h_normal_f, levels=[maximum(data_h_normal_f)/exp(2)], color= "white", linestyle = (:dot, 7), label = "1/e² contour")

    return hm_prod
end


#plot
using CairoMakie
f_hr_prod = Figure(size = (600, 400))
axs_prod = Axis(f_hr_prod[1, 1], 
    xlabel = "Radial Distance [m]", 
    ylabel = "Height [km]", 
    limits = ((nothing, 80), (70, 130)),
    title = "isotropic"
    )
hm_prod = hist_radial_heihgt_production(axs_prod, his_hrp)
      
Colorbar(f_hr_prod[1, 2], hm_prod, label="Production [m⁻³]")
axislegend(axs_prod)
dsiplay(f_hr_prod)

save(joinpath(dir, "plots", "hist_norm_hr_all_E_pitch.png"), f_hr_prod, px_per_unit = 3.3)




## OPTIONAL: plotting gyroradius
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


if lim_pitch_deg == 90
    r_gyro_max_pitch_lim = v_abs(E0) * c.me ./ (c.qe * B)
end

        
lines!(axs_prod, r_gyro_max_pitch_lim, h_middle/1e3, color= "black", linestyle = (:dot, 10), linewidth = 1, label = "Max gyroradius")
dsiplay(f_hr_prod)
