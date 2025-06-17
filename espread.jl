using Distributions
using LinearAlgebra

include("get_msis.jl")
include("magnetic_field.jl")
include("constants.jl")

using GLMakie
hist([rand(Exponential(1)) for i in 1:100000], bins=1000)

#geomagnetic location of tromso:
# source https://eiscat.se/scientist/document/location-of-the-eiscat-facilities/
loc_gmag = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
lat_rad = loc_gmag[1] * pi /180
loc_geod = [69.58, 19.23]

#start point: gyrocenter at t= 0
gc0 = (re + 600e3) .* [0, cos(lat_rad), sin(lat_rad)]

#start velocity
#energy 
E = 1e3 #eV
v = sqrt(2*E*qe/me)         # convert to veloicty => relativistic?
pitch = rand()*pi/2         # pitch angle
v_par  = v*cos(pitch)       # calculate parallel component
# for parallel component, find normal vecotors, get random phase, distribute velocity accordingly
# velocity vector starts a gyroradius away from gyrocenter
# maybe do it the other way around, first declare v0, r0, find gyrocenter after?
v_perp = sqrt(v^2 - v_par^2)
B0 = dipole_field_earth(gc0)
r_gyro = me * v_perp / (qe * norm(B0))
#create orthogonal vectorsystem along B0
xu = [1, 0, 0]
n1 = cross(xu, B0)
n1u = n1/norm(n1)
n2 = cross(n1u, B0)
n2u =n2/norm(n2)
# phase angle
phase = rand()*2*pi
# calculate r0 and v0 in n1, n2 directions, and convert to original coordinate system
r_n1 =  cos(phase) * r_gyro .* n1u
r_n2 =  sin(phase) * r_gyro .* n2u
v_n1 = -sin(phase) * v_perp .* n1u
v_n2 =  cos(phase) * v_perp .* n2u
v_n3 =  v_par .* B0 / norm(B0)
r0 = r_gyro .+ r_n1 .+ r_n2
v0 = v_n1 .+ v_n2 .+ v_n3
# plotting:
fig = cm.arrows([Point3f([0, 0, 0])], [Vec3f(B0./norm(B0))], arrowcolor = :red, label = "adf")
axislegend()
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(n1u)], arrowcolor = :blue)
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(n2u)], arrowcolor = :green)
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(r_n1 .+ r_n2)])
cm.arrows!([Point3f(r_n1 .+ r_n2)], [Vec3f(v_n1 .+ v_n2)./norm(v0)])
cm.arrows!([Point3f(r_n1 .+ r_n2)], [Vec3f(v0)./norm(v0)])
display(GLMakie.Screen(), fig)

r0v0 = [r0; v0]
#magnetic field handle for boris mover: (static field)
Bin(t, p) = dipole_field_earth(p)

#mean free path calculations
#msis altitude resolution could be coarser??
altitude = norm(gc0) - re
h = 80e3:1e3:600e3 #km
atm_matrix = msis([[2020, 12, 12, 18, 0, 0]], h, loc_geod[1], loc_geod[2])[:]

atm = stack([[a.N2_number_density for a in atm_matrix],
             [a.O2_number_density for a in atm_matrix],
             [a.O_number_density for a in atm_matrix]])

nN2 = [a.N2_number_density for a in atm_matrix]
nO2 = [a.O2_number_density for a in atm_matrix]
nO  = [a.O_number_density  for a in atm_matrix]

"""
pn = propertynames(atm[1])
atm = [getproperty(a, p) for a in atm_matrix, p in pn]

f = Figure()
ax = Axis(f[1, 1], xscale=log10)
for p in pn
    lines!([getproperty(a, p) for a in atm], h, label=string(p))
end
axislegend()
"""

# interpolate atmospheric desnities
# function to find closest atmospheric density: 
#    can we make a rather simple one, 
#    fitting a low order polynom to the lig(density)? 
#    should eb faster than pchipinterpolation

using PCHIPInterpolation
nN2_ip(hh) = exp(Interpolator(h, log.(nN2))(hh))
nO2_itp = Interpolator(h, nO2)
nO_itp  = Interpolator(h, nO)

lines(nN2, h, axis=(xscale=log10,),)
scatter(nN2_ip.(80e3:0.1e3:600e3), 80e3:0.1e3:600e3, markerstyle="x")


lmfp = sum_i(sigma * ni) # *sqrt(2) #sqrt(2) for relative velocity, if electron and gas have compareable velocities. electron is much lighter, so no!
i = species

#selecting the path length until the next collision:
l_nc = rand(Exponential(1/lmfp))

t = l/v0
t ????? => replace by path length in boris mover
#initalize boris mover
ode_boris_mover(t, r0v0, qe, me, [0, 0, 0], Bin)



#to do:
#cross sections, mean free path
# ni depends on altitude, calculate it using norm(gc0) - re, get msis for neutral atmosphere
# cross section depends on energy, create lookup table?
# 