using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase
using PCHIPInterpolation: Interpolator
SpaceIndices.init()

function msis(times, heights, lats, longs)
    jdate = [date_to_jd(tt...) for tt in times]
    atm = [AtmosphericModels.nrlmsise00(
            jd,         # Julian Day
            h,          # Altitude [m]
            lat,         # Latitude [rad]
            long          # Longitude [rad]
            ) for jd in jdate, h in heights, lat in lats, long in longs];
    return atm
end

function atmospheric_model(times, heights, lats, longs)
    atm_matrix = msis(times, heights, lats, longs)[:]
    nN2 = [a.N2_number_density for a in atm_matrix]
    nO2 = [a.O2_number_density for a in atm_matrix]
    nO  = [a.O_number_density  for a in atm_matrix]
    """
    pn = propertynames(atm[1])
    atm = [getproperty(a, p) for a in atm_matrix, p in pn]  

    f = Figure()
    ax = Axis(f[1, 1], xscale=log10)
    for p in pn
        lines!([getproperty(a, p) for a in atm], hmsis, label=string(p))
    end
    axislegend()
    """ 

    # interpolate atmospheric desnities
    #    can we make a rather simple one, 
    #    fitting a low order polynom to the lig(density)? 
    #    should eb faster than pchipinterpolation
    nN2_ip(hh) = exp(Interpolator(heights, log.(nN2))(hh))
    nO2_ip(hh) = exp(Interpolator(heights, log.(nO2))(hh))
    nO_ip(hh)  = exp(Interpolator(heights, log.(nO))(hh))
    densities(hh) = [nN2_ip(hh), nO2_ip(hh), nO_ip(hh)]
    return densities
    """
    f, ax, lin = scatter(nO, heights, axis=(xscale=log10,),)
    lines!(nO_ip.(80e3:0.1e3:600e3), 80e3:0.1e3:600e3)
    """
end


function atmospheric_model_fast(times, heights, lats, longs)
    atm_matrix = msis(times, heights, lats, longs)[:]
    nN2 = [a.N2_number_density for a in atm_matrix]
    nO2 = [a.O2_number_density for a in atm_matrix]
    nO  = [a.O_number_density  for a in atm_matrix]
    return stack([nN2, nO2, nO])
end
"""
function consolidate_msis(atm)
    #function to consolidate results of msis() when called with lists as inputs
    # takes matrix of msis results (atm) and iterates through them, extracting all properties
    # and creates a new model, with the same properties, that contain matrices 
    pn = propertynames(atm[1])
    for a in atm
        for p in pn
            getproperty(a, p)
        end
    end
"""

function make_densityf(hmin, hmax, hintervals)
    hmsis = hmin:hintervals:hmax
    atm = Float64.(atmospheric_model_fast([[2020, 12, 12, 18, 0, 0]], hmsis, loc_geod[1], loc_geod[2]))
    return densityf_fast(alt) = atm[round(Int, (alt-hmin)/hintervals+1), :]
end
