#load cross sections
include("data_electron/e_N2_cross_sections.jl")
include("data_electron/e_O2_cross_sections.jl")
include("data_electron/e_O_cross_sections.jl")

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

    scat_par = Matrix{Any}(undef, length(function_name), 4)
    scat_par[:, 1] = function_name
    scat_par[:, 2] = energy_levels
    scat_par[:, 3] = n_secondary_el
    scat_par[:, 4] = [getfield(Main, f) for f in function_handle] #assign functions from x_cross_sections.jl
    
    return scat_par
end

# get scattering parameters for all species:
sp_N2 = get_scattering_parameters("N2")
sp_O2 = get_scattering_parameters("O2")
sp_O  = get_scattering_parameters("O")

# cross section functions
cs_e_N2(Ep) = only.([cs([float(Ep)]) for cs in sp_N2[:, 4]])
cs_e_O2(Ep) = only.([cs([float(Ep)]) for cs in sp_O2[:, 4]])
cs_e_O(Ep)  = only.([cs([float(Ep)]) for cs in sp_O[:, 4]])

#stacked, so it can be multiplies by a vector of densities: 
# cs_all(Ep) .* [nN2, nO2, nO]
# be aware of ordering!
cs_all(Ep) = [cs_e_N2(Ep), cs_e_O2(Ep), cs_e_O(Ep)]

# summed for total cross section:
cs_all_sum(Ep) = sum.(cs_all(Ep))

#aliases, remove
#sigma_tot(Ep) = cs_all_sum(Ep)
#sigma_tot_list(Ep) = cs_all(Ep)    
