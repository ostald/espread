using Distributions
using LinearAlgebra

include("magnetic_field.jl")
include("constants.jl")
include("local_ortognal_basis.jl")
include("ode_boris_mover.jl")
#load cross sections
include("cross_sections.jl")
include("get_msis.jl")


using GLMakie
"""
cm = GLMakie
hist([rand(Exponential(1)) for i in 1:100000], bins=1000)
"""

#geomagnetic location of tromso:
# source https://eiscat.se/scientist/document/location-of-the-eiscat-facilities/
loc_gmag = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
lat_rad = loc_gmag[1] * pi /180
loc_geod = [69.58, 19.23]

#start point: gyrocenter at t= 0
gc0 = (c.re + 500e3) .* [0, cos(lat_rad), sin(lat_rad)]
altitude = norm(gc0) - c.re


#start velocity
#energy must be larger than 15.6
E = 10000.0 #eV must be float!
v = sqrt(2*E*c.qe/c.me)         # convert to veloicty => relativistic?
pitch = rand()*pi/2         # pitch angle
v_par  = v*cos(pitch)       # calculate parallel component
# for parallel component, find normal vecotors, get random phase, distribute velocity accordingly
# velocity vector starts a gyroradius away from gyrocenter
# maybe do it the other way around, first declare v0, r0, find gyrocenter after?
v_perp = sqrt(v^2 - v_par^2)
B0 = dipole_field_earth(gc0)
r_gyro = c.me * v_perp / (c.qe * norm(B0))
#create orthogonal vectorsystem along B0
u1, u2, u3 = local_orthogonal_basis(B0)
# phase angle
phase = rand()*2*pi
# calculate r0 and v0 in n1, n2 directions, and convert to original coordinate system
r_n1 =  cos(phase) * r_gyro .* u1
r_n2 =  sin(phase) * r_gyro .* u2
v_n1 = -sin(phase) * v_perp .* u1
v_n2 =  cos(phase) * v_perp .* u2
v_n3 =  v_par .* u3
r0 = gc0 .+ r_n1 .+ r_n2
v0 = v_n1 .+ v_n2 .+ v_n3
altitude = norm(r0) - c.re

"""
# plotting start positions, velocity and magnetic fieød
fig = cm.arrows([Point3f([0, 0, 0])], [Vec3f(B0./norm(B0))], arrowcolor = :red, label = "adf")
axislegend()
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(n1u)], arrowcolor = :blue)
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(n2u)], arrowcolor = :green)
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(r_n1 .+ r_n2)])
cm.arrows!([Point3f(r_n1 .+ r_n2)], [Vec3f(v_n1 .+ v_n2)./norm(v0)])
cm.arrows!([Point3f(r_n1 .+ r_n2)], [Vec3f(v0)./norm(v0)])
display(GLMakie.Screen(), fig)
"""

r0v0 = [r0; v0]
#magnetic field handle for boris mover: (static field)
Bin(p) = dipole_field_earth(p)

#mean free path calculations
#msis altitude resolution could be coarser??
hmsis = 80e3:1e3:600e3 #km
densityf = atmospheric_model([[2020, 12, 12, 18, 0, 0]], hmsis, loc_geod[1], loc_geod[2])


# sample number of mean free paths travelled:
n_mfp = rand(Exponential(1))
r, v, t = ode_boris_mover_mfp(n_mfp, r0v0, -c.qe, c.me, Bin, cs_all_sum, densityf)


# cumulative sum of cross section * density to decide scattering partner:
cross_sections = cs_all_sum(E)
ns = densityf(norm(r) - c.re)
# sample uniformly between scattering partners:
r_scatter = rand(1)[1]*sum(cross_sections .* ns)
#species index: according to densityf [N2, O2, O]
idx_species = findfirst(cumsum(cross_sections .* ns) .> r_scatter)


# alternatively, use cross sections of ionization processes directly, not summed for each species
# thereby one radnom number decides directly which scattering process is triggered
ns = densityf(norm(r) - c.re)
# 
scatter_p = vcat((cs_all(E) .* ns)...)
r_scatter = rand(1)[1]*sum(scatter_p)
idx_scatter = findfirst(cumsum(scatter_p) .> r_scatter)

"""
nsample = Int(1e7)
hist([findfirst(cumsum(scatter_p) .> r_scatter) for r_scatter in rand(nsample).*sum(scatter_p)],
    bins = 0.5:1:63.5,
    scale_to=scatter_p[1]
    )
"""

# scattering:
# scattering phase is uniform:
scatter_phase = rand(1)*2*pi
# scattering angle along flightpath:
# depends on secondary electron eenergy


#do collision:
#1. decide which particle: sum densities, normalize,

