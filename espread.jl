using Distributions
using LinearAlgebra

include("magnetic_field.jl")
include("constants.jl")

using GLMakie
hist([rand(Exponential(1)) for i in 1:100000], bins=1000)

#geomagnetic location of tromso:
# source https://eiscat.se/scientist/document/location-of-the-eiscat-facilities/
loc_gmag = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
lat_rad = loc_gmag[1] * pi /180

#start point: gyrocenter at t= 0
gc0 = (6730e3 + 600e3) .* [0, cos(lat_rad), sin(lat_rad)]

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
#initalize boris mover

t ?????

ode_boris_mover(t, r0v0, qe, me, [0, 0, 0], Bin)



#selecting the path length until the next collision:
rand(Exponential(lam))


#to do:
#cross sections, mean free path