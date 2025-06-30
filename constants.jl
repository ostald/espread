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
    return v_abs^2 * c.me  / (2 * c.qe) 
end

function v_abs(E_ev)
    return sqrt(E_ev * c.qe * 2 / c.me)
end