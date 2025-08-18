#using PhysicalConstants

c = (
# electron mass
me = 9.10938356e-31	#kg
,
# electron charge
qe = 1.6021766208e-19	#C
,
# earth radius
re = 6378e3  #m
,
vm = [0, 0, -8.22e15] # Tm3
,
)


function makeE_ev(c)
    function E_ev(v_abs::Float64)
        # energy in eV of an electron travelling with |v| = v_abs
        return v_abs^2 * c.me  / (2 * c.qe) 
    end
    return E_ev
end

E_ev = makeE_ev(c)

function make_v_abs(c)
    function v_abs(E_ev::Float64)
        # velocity of an electron with kinetic energy E_ev
        return sqrt(E_ev * c.qe * 2 / c.me)
    end
    return v_abs
end

v_abs = make_v_abs(c)

function make_altitude(c)
    function altitude(r::Vector{Float64})
        # altitude above earth's surface of a particle at coordinates r
        # with origin at geomagnetic center of earth
        return norm(r) - c.re
    end
    return altitude
end

"""
function v_abs(E_ev, c)
    # velocity of an electron with kinetic energy E_ev
    return sqrt(E_ev * c.qe * 2 / c.me)
end
"""

altitude = make_altitude(c)
