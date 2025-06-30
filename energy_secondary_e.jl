using PCHIPInterpolation: Interpolator

function E_secondary_e_N2(E_primary, E_ionization)
    # returns a pdf of the secondary electron energy
    # limited by E_max
    # follows a Lorentz/Cauchy Distribution, but we can simply use Opal's approxiamtion
    # see Opal et al 1972
    E_hat = 11.4  #eV
    E_max = (E_primary - E_ionization)/2
    ps(Es) = 1/(Es^2 + E_hat^2)
    #normalization
    #trapeziodal integration: see https://se.mathworks.com/help/matlab/ref/trapz.html#bua4lsr
    # N+1 evenly spaced points, with N intervals inbetween
    # int_a^b f(x) dx ≈ (b-a)/2N * (f(x1) + 2f(x2) + ... + 2f(xN) + f(xN+1))
    ee = LinRange(0, E_max, Int(ceil(E_max))) # array [0, E_max] with approx 1eV steps
    n_ee = length(ee)-1
    int = E_max * (sum(ps.(ee)) + sum(ps.(ee[2:end-1]))) / 2 / n_ee
    normalization_factor = 1/int 
    p_secondary_e(E_secondary) = normalization_factor ./ (E_hat^2 .+ E_secondary.^2) .* (E_secondary .< E_max)

    return p_secondary_e, E_max
end


function E_secondary_e_O2(E_primary, E_ionization)
    # returns a pdf of the secondary electron energy
    # limited by E_max
    # follows a Lorentz/Cauchy Distribution, but we can simply use Opal's approxiamtion
    # see Opal et al 1972
    E_hat = 15.2 #eV
    E_max = (E_primary - E_ionization)/2
    ps(Es) = 1/(Es^2 + E_hat^2)
    #normalization
    #trapeziodal integration: see https://se.mathworks.com/help/matlab/ref/trapz.html#bua4lsr
    # N+1 evenly spaced points, with N intervals inbetween
    # int_a^b f(x) dx ≈ (b-a)/2N * (f(x1) + 2f(x2) + ... + 2f(xN) + f(xN+1))
    ee = LinRange(0, E_max, Int(ceil(E_max)))
    n_ee = length(ee)-1
    int = E_max * (sum(ps.(ee)) + sum(ps.(ee[2:end-1]))) / 2 / n_ee
    normalization_factor = 1/int 
    p_secondary_e(E_secondary) = normalization_factor ./ (E_hat^2 .+ E_secondary.^2) .* (E_secondary .< E_max)

    return p_secondary_e, E_max
end

"""
# to produce cdf, use fine a fine resolution when evaluating the pdf,
# use cumsum for simlpe quadrature integration,
# producing a discrete cdf at points Es_lim
# use cdf_discrete to randomly choose E_secondary with resolution dEs
pdf, E_max = E_secondary_e_O2(50, 12.8)
dEs = 0.1
Es_lim =  0:dEs:E_max+dEs
E_secondary = Es_lim[1:end-1] .+ dEs/2 #evaluate cdf in teh middle of a bin
pdf_discrete = pdf(E_secondary)

# cdf[1] == pdf[1] != 0  => cdf does not start with 0, instead first 
# element of cdf stands for the integrated pdf from 0:dEs, 
# i.e. first energy bin, ie. Es_lim[1] - Es_lim[2], 
# with mean energy E_secondary[1]. This allows the expression in the 
# histogram to succeed, as rand() > 0, so cdf > rand makes no sense if cdf has 0. 
# ALternatively, we could use Es_lim in histogram, and include 0 in cumsum.
cdf_discrete = cumsum(pdf_discrete) 
cdf_discrete = cdf_discrete ./ cdf_discrete[end] #cdf must start from 0

lines(E_secondary, pdf_discrete)
scatter!(Es_lim[2:end], cdf_discrete)

nsample = Int(1e6)
#specify correct binning!
hist([E_secondary[findfirst(cdf_discrete .> r)] for r in rand(nsample)],
    bins = Es_lim #0.5:1:E_max,
    )
lines!(E_secondary, pdf_discrete*nsample*dEs)
display(current_figure())
"""

# table 6.1 itikawa et al 1990
# [incident energy [eV]     A [eV]      B [1e-18 cm^2 / eV]]
secondary_e_O_table = [
        100     12.6    7.18
        200     13.7    4.97
        500     14.1    2.75
        1000    14.0    1.69
        2000    13.7    1.02
    ]

A_itp = Interpolator(secondary_e_O_table[:, 1], secondary_e_O_table[:, 2])


function E_secondary_e_O(E_primary, E_ionization)
    # returns a pdf of the secondary electron energy
    # limited by E_max
    # follows a Lorentz/Cauchy Distribution, but we can simply use Opal's approxiamtion
    # see Opal et al 1972
    # Opal is missing atomic oxigen, see Itikawa 1990 instead
    #table 6.1 itikawa et al 1990
    # [incident energy [eV]     A [eV]      B [1e-18 cm^2 / eV]]
    """
    secondary_e_O_table = [
        100     12.6    7.18
        200     13.7    4.97
        500     14.1    2.75
        1000    14.0    1.69
        2000    13.7    1.02
    ]

    E_incident = secondary_e_O_table[:, 1]
    A = secondary_e_O_table[:, 2]
    #B = secondary_e_O_table[:, 3]
    A_itp = Interpolator(E_incident, A)
    #B_itp = Interpolator(E_incident, B)

    fig, ax, lin = lines(E_incident, A, axis = (limits = (nothing, nothing),),)
    #lines!(E_incident, B)
    ee = 100:2000
    lines!(ee, A_itp.(ee)) 
    #lines!(ee, B_itp.(ee))
    display(current_figure())
    """

    if E_primary < 100
        A = 12.6  #eV
    elseif E_primary > 2000
        A = 13.7 #eV
    else
        A = A_itp(E_primary)
    end
    
    E_max = (E_primary - E_ionization)/2
    ps(Es) = B ./(1 .+(Es ./A).^(5/3))
    
    #normalization
    #trapeziodal integration: see https://se.mathworks.com/help/matlab/ref/trapz.html#bua4lsr
    # N+1 evenly spaced points, with N intervals inbetween
    # int_a^b f(x) dx ≈ (b-a)/2N * (f(x1) + 2f(x2) + ... + 2f(xN) + f(xN+1))
    ee = LinRange(0, E_max, Int(ceil(E_max)))
    n_ee = length(ee)-1
    int = E_max * (sum(ps.(ee)) + sum(ps.(ee[2:end-1]))) / 2 / n_ee
    normalization_factor = 1/int 
    
    p_secondary_e(E_secondary) = normalization_factor .* ps.(E_secondary) .* (E_secondary .< E_max)

    return p_secondary_e, E_max
end

"""
pdf, E_max = E_secondary_e_O(1e4, 13.16)
# energy resolution can be made dynamic. 
# for small energies, a resolution of 0.1 is recommended, but as the 
# primary electron energy gets large, so does Es_max. 
# for large secondary energy values (>100eV?), much bigger steps can be taken.
dEs = 0.1
Es_lim =  0:dEs:E_max+dEs
E_secondary = Es_lim[1:end-1] .+ dEs/2 #evaluate cdf in the middle of a bin
pdf_discrete = pdf(E_secondary) 
cdf_discrete = cumsum(pdf_discrete)
cdf_discrete = cdf_discrete ./ cdf_discrete[end] #cdf must start from 0

lines(E_secondary, pdf_discrete)
scatter!(Es_lim[2:end], cdf_discrete)

nsample = Int(1e4)
hist([E_secondary[findfirst(cdf_discrete .> r)] for r in rand(nsample)],
    bins = Es_lim#0:1:E_max,
    )
lines!(E_secondary, pdf_discrete*nsample*dEs)
display(current_figure())
"""