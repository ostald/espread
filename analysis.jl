using CSV
using DataFrames
using LinearAlgebra
include("constants.jl")


dir = readdir("results/")
filter_crit = "res"
folders = filter(x-> contains(x, filter_crit), dir)

res = Vector{Any}(undef, length(folders))

#for (id, f) in enumerate(folders)
    df = CSV.read(open(joinpath("results", f, "results.txt")),
        header=false,
        delim = "\t",
        skipto=4,
        DataFrame
        )
    rename!(df, [:r0, :v0, :status, :r, :v])
    df.r0 = [eval(Meta.parse(value)) for value in df.r0]
    df.v0 = [eval(Meta.parse(value)) for value in df.v0]
    df.r = [eval(Meta.parse(value)) for value in df.r]
    df.v = [eval(Meta.parse(value)) for value in df.v]
    res[id] = df
end

# check for fails in Boris mover:
if size(filter(row -> :status .== 0, df))[1] > 0
    error("Failure of Boris mover!")
end

# primary electrons:
altitude.(df.r0)
filter(:r0 => x -> altitude(x) > 550e3, df)

df.E0 = E_ev.(norm.(df.v0))
df.E_end = E_ev.(norm.(df.v))

#check for status -1 and Energy larger than lowest ionization energy
filter(row -> row.E0 < 12.072 && row.status != -1, df)


df.alt0 = altitude.(df.r0)
df.alt_end = altitude.(df.r)
filter(:alt0 => x -> x > 599e3, df)

# select endpoint for primary electrons, starting point for secondary electrons
df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)
df.pos = ifelse.(df.alt0 .> 599e3, df.r, df.r0)



hmin = 80e3     #m
hmax = alt0+1e4 #m
intervals = 1e3 #m
hmsis = hmin:intervals:hmax

using GLMakie
fig, ax, his = hist(df.alt./1e3,
    bins = hmsis./1e3, 
    direction=:x,
    axis = (xlabel = "Eelctron Count [1]",
        ylabel = "Height [km]",
        limits = ((1, nothing), nothing),
        xscale = log10,
        ),
    )
display(fig)
save(joinpath("results", folders[1], "hist_ionizations_height.png"), fig)



fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1))
sc = meshscatter!(ax,Point3.(df.pos), markersize = 1e2)#, 
    #axis=(limits=(nothing, nothing, nothing),),)
zlims!(ax, 5.99e6, 6.01e6)
ylims!(ax, 2.465e6, 2.48e6)
xlims!(ax, -1.5e4/2, 1.5e4/2)


meshscatter(Point3.((filter(row -> row.alt > 105e3 && row.alt < 106e3, df).pos)), markersize = 1e1)

# transform into field-aligned coordinates
include("magnetic_field.jl")
loc_gmag_deg = [66.73, 102.18]   # degrees lat, long, geomagnetic coord.
loc_gmag = loc_gmag_deg ./ 180 * pi

loc_geod_deg = [69.58, 19.23]
loc_geod = loc_geod_deg ./ 180 * pi

lat_gmag = loc_gmag[1]

#straight up:
gc = [(c.re + h) * [0, cos(lat_gmag), sin(lat_gmag)] for h in hmsis]
    
B0 = dipole_field_earth.(gc)

or = [0, 0, 0]
fig = arrows3d(Point3(or), Vec3(B0[1]), alpha = 0.3)#, tipcolor = :red, label = "B/|B|")
arrows3d!(Point3(or), Vec3(B0[end]))
#magnetic field cunrvature can be neglected

#create orthogonal vectorsystem along B0
u1, u2, u3 = local_orthogonal_basis(B0[1])

origin = (c.re + hmsis[1]) * [0, cos(lat_gmag), sin(lat_gmag)]

#convert into magnetic field aligned coordinates
mat = inv([u1 -u2 -u3])
df.pos_fa = [mat * (p-origin) for p in  df.pos]

fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1))
sc = meshscatter!(ax, Point3.((filter(row -> row.alt > 105e3 && row.alt < 106e3, df).pos_fa))./1e3, markersize = 1e-3)#, 
    #axis=(limits=(nothing, nothing, nothing),),)
zlims!(ax, -2e5, -1e5)
ylims!(ax, -1.5e5, -0.5e5)
ylims!(ax, -3, 3)

meshscatter(Point3.((filter(row -> row.alt > 105e3 && row.alt < 106e3, df).pos_fa)), markersize = 1e1)
