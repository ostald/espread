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
)

function E_ev(v_abs)
    # energy in eV of an electron travelling with |v| = v_abs
    return v_abs^2 * c.me  / (2 * c.qe) 
end

function v_abs(E_ev)
    # velocity of an electron with kinetic energy E_ev
    return sqrt(E_ev * c.qe * 2 / c.me)
end

function altitude(r)
    # altitude above earth's surface of a particle at coordinates r
    # with origin at geomagnetic center of earth
    return norm(r) - c.re
end