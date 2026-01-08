# Setup file for individual runs.

using Dates

#magnetic field model
b_model = "vertical"
#b_model = "dipole"

# Define initial energy
#E0 = 1e4 #eV
#e_energy = [500, 1e3, 2e3, 4e3, 8e3]
e_energy = [16e3]

# Define number of particles
N_electrons = 1e6

# starting altitude
alt0 = 600e3 #m

# Define pitch angle limits
#lim_pitch_deg = 20 #deg
#lim_pitch = lim_pitch_deg/180*pi
#Bin! = dipole_field_earth!
#Bin! = convergent_vertical_field!

pitch_limits_deg = [20, 90]
#pitch_limits_deg = [60]

nPerGyro = 20

# initiate seed
#Random.seed!(Int(E0)) # use energy as seed so that runs of the same energy are reproducible


#geomagnetic location of tromso:
# source https://eiscat.se/scientist/document/location-of-the-eiscat-facilities/
loc_gmag_deg = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
loc_gmag = loc_gmag_deg ./ 180 * pi

loc_geod_deg = [69.58, 19.23]
loc_geod = loc_geod_deg ./ 180 * pi

name = "r7_conicB_16kev"
res_dir = joinpath("results", name * string(now()))

nprocesses = 200



"""
TO DO

- Double ion not handled correctly yet
- some tracing?
"""