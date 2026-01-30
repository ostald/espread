#load cross sections
include("data_electron/e_N2_cross_sections.jl")
include("data_electron/e_O2_cross_sections.jl")
include("data_electron/e_O_cross_sections.jl")
include("scattering.jl")

using DelimitedFiles: readdlm

function get_scattering_parameters(species_name)
    # for every species name, load info on cross sections and corresponding energy loss 
    filename = joinpath(pwd(), "data_neutrals", species_name * "_levels.name")
    state_name = readdlm(filename, String, comments=true, comment_char='%')
    function_name = "e_" * species_name .* state_name
    function_handle = Symbol.(function_name)

    filename2 = joinpath(pwd(), "data_neutrals", species_name * "_levels.dat")
    energy_secondary_matrix = readdlm(filename2, comments=true, comment_char='%')
    energy_levels = energy_secondary_matrix[:, 1]
    n_secondary_el = Int.(energy_secondary_matrix[:, 2])

    scat_par = Matrix{Any}(undef, length(function_name), 6)
    scat_par[:, 1] = function_name
    scat_par[:, 2] = energy_levels
    scat_par[:, 3] = n_secondary_el
    scat_par[:, 4] = [getfield(Main, f) for f in function_handle] #assign functions from x_cross_sections.jl
    scat_par[:, 5] .= species_name
    
    
    if species_name == "N2"
        scat_par[:, 6] .= scatter_inelastic_N2
        scat_par[1, 6]  = scatter_elastic_N2
        scat_par[scat_par[:, 3] .> 0, 6] .= scatter_ion_N2
    elseif species_name == "O2"
        scat_par[:, 6] .= scatter_inelastic_O2
        scat_par[1, 6]  = scatter_elastic_O2
        scat_par[scat_par[:, 3] .> 0, 6] .= scatter_ion_O2
    elseif species_name == "O"
        scat_par[:, 6] .= scatter_inelastic_O
        scat_par[1, 6]  = scatter_elastic_O
        scat_par[scat_par[:, 3] .> 0, 6] .= scatter_ion_O
    end
    
    """
    if species_name == "N2"
        scat_par[:, 6] .= ["N2", "inelastic"] 
        scat_par[1, 6]  = ["N2", "elastic"]
        scat_par[scat_par[:, 3] .> 0, 6] .= ["N2", "ionizing"]]
    elseif species_name == "O2"
        scat_par[:, 6] .= ["O2", "inelastic"] 
        scat_par[1, 6]  = ["O2", "elastic"]
        scat_par[scat_par[:, 3] .> 0, 6] .= ["O2", "ionizing"]
    elseif species_name == "O"
        scat_par[:, 6] .= ["O", "inelastic"] 
        scat_par[1, 6]  = ["O", "elastic"]
        scat_par[scat_par[:, 3] .> 0, 6] .= ["O", "ionizing"]
    end
    """

    return scat_par
end

## for typestable code, avoid global variables!
# => solution is a function barrier

# get scattering parameters for all species:
sp_N2 = get_scattering_parameters("N2")
sp_O2 = get_scattering_parameters("O2")
sp_O  = get_scattering_parameters("O")

# cross section functions
cs_e_N2(Ep) = only.([cs([float(Ep)]) for cs in sp_N2[:, 4]])
cs_e_O2(Ep) = only.([cs([float(Ep)]) for cs in sp_O2[:, 4]])
cs_e_O(Ep)  = only.([cs([float(Ep)]) for cs in sp_O[:, 4]])

# stacked, so it can be multiplies by a vector of densities: 
# cs_all(Ep) .* [nN2, nO2, nO]
# be aware of ordering!
cs_all(Ep) = [cs_e_N2(Ep), cs_e_O2(Ep), cs_e_O(Ep)]::Vector{Vector{Float64}}
sp_all = [sp_N2; sp_O2; sp_O]

# summed for total cross section:
cs_all_sum(Ep) = sum.(cs_all(Ep))

#aliases, remove
#sigma_tot(Ep) = cs_all_sum(Ep)
#sigma_tot_list(Ep) = cs_all(Ep)    


"""
using CairoMakie

Ep = logrange(1e-1, 1e6, 700)

for sig in eachrow(sp_all)
    f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e5), (1e-25, 1e-18)),
        ),
    )
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
    axislegend(a)
    display(f)
end

f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e6), (1e-25, 1e-18)),
        ),
    )
for sig in eachrow(sp_all)
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
end
axlislegend(a)
display(f)
save("figures/cross_sections.png", f)


f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e5), (1e-25, 1e-18)),
        ),
    )
for sig in eachrow(sp_N2)
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
end
axislegend(a)
display(f)
save("figures/cross_sections_N2.png", f)


f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e5), (1e-25, 1e-18)),
        ),
    )
for sig in eachrow(sp_O2)
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
end
axislegend(a)
display(f)
save("figures/cross_sections_O2.png", f)


f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e5), (1e-25, 1e-18)),
        ),
    )
for sig in eachrow(sp_O)
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
end
axislegend(a)
display(f)
save("figures/cross_sections_O.png", f)


f, a, lin = scatter(1, 1, 
    axis= (xscale = log10,
        yscale = log10,
        limits = ((1e-1, 1e5), (1e-25, 1e-18)),
        ),
    )
for sig in eachrow(sp_He)
    if sig[2] == 0.0
        lines!(a, Ep, sig[4](Ep), label = "elastic", linewidth = 2)
    elseif sig[3] == 0
        lines!(a, Ep, sig[4](Ep), label = sig[1])
    else
        lines!(a, Ep, sig[4](Ep), label = sig[1], linestyle = :dashdot)
    end
end
axislegend(a)
display(f)
save("figures/cross_sections_He.png", f)
"""