using SatelliteToolboxAtmosphericModels
using SatelliteToolboxBase
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

