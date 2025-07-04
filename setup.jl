# Setup file for individual runs.

using Dates

# Define initial energy
E0 = 50 #eV

# Define number of particles
N_electrons = 1e6

# starting altitude
alt0 = 600e3 #m

# Define pitch angle limits
lim_pitch_deg = 90 #deg
lim_pitch = 90/180*pi

# initiate seed
Random.seed!(E0) # use energy as seed so that runs of the same energy are reproducible


# Results directory
res_dir = joinpath("results", string(now()))


#geomagnetic location of tromso:
# source https://eiscat.se/scientist/document/location-of-the-eiscat-facilities/
loc_gmag_deg = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
loc_gmag = loc_gmag_deg ./ 180 * pi

loc_geod_deg = [69.58, 19.23]
loc_geod = loc_geod_deg ./ 180 * pi