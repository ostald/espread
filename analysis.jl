using CSV
using DataFrames
using LinearAlgebra
using JLD2
include("constants.jl")

#dir = "results/run1_2025-07-19T21:53:59.887/"
dir = "results/conicB/"
dir_con = readdir(dir)

runs = unique([d[1:end-8] for d in dir_con])

for r in runs
    #filter_crit = "res_1000.0eV_20.0deg"
    filter_crit = r
    files = filter(x-> contains(x, filter_crit), dir_con)
    
    res = Vector{Any}(undef, length(files))
    
    #f = files[2]
    for (id, f) in enumerate(files)
        data = load(dir * f)
        all_keys = collect(keys(data))
        df_keys = filter(x-> !contains(x, "setup"), all_keys)

        df = vcat([data[k] for k in sort(df_keys)]...)
        unique(df)

        """
        df = CSV.read(open(joinpath(dir, f)),
            header=false,
            delim = "\t",
            skipto=8,
            DataFrame,
            )
        rename!(df, [:r0, :v0, :status, :r, :v])
        df.r0 = [eval(Meta.parse(value)) for value in df.r0]
        df.v0 = [eval(Meta.parse(value)) for value in df.v0]
        df.r  = [eval(Meta.parse(value)) for value in df.r]
        df.v  = [eval(Meta.parse(value)) for value in df.v]
        """
        res[id] = df
    end
    df_comb = unique(vcat(res...))
    
    #jldsave(joinpath(dir, r * ".jld2"); df_comb)

end

df_dict = load(joinpath(dir, "res_1000.0eV_20.0deg.jld2"))
df = unique(df_dict["df_comb"])

res_file = "results/conicB_run2_2025-08-26T11:19:42.742/res_10000.0eV_20.0deg_000.bin"
df = DataFrame(Generation = Int[], idx_scatter = Int[], r0 = Vector{Float64}[], v0 = Vector{Float64}[], status = Int[], r = Vector{Float64}[], v = Vector{Float64}[])
io = open(res_file, "r")
E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals = deserialize(io)
while !eof(io)
    push!(df, deserialize(io))
end
close(io)


# check for fails in Boris mover:
if size(filter(row -> :status .== 0, df))[1] > 0
    error("Failure of Boris mover!")
end

df.E0 = E_ev.(norm.(df.v0))
df.E_end = E_ev.(norm.(df.v))

#check for status -1 and Energy larger than lowest ionization energy
filter(row -> row.E0 < 12.072 && row.status != -1, df)


# primary electrons:
df.alt0 = altitude.(df.r0)
df.alt_end = altitude.(df.r)
df_alt0 = filter(:alt0 => x -> x > 599e3, df)

fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax, Point3.(df_alt0.r0), markersize = 1)#, 


# select endpoint for primary electrons, starting point for secondary electrons
df.alt = ifelse.(df.alt0 .> 599e3, df.alt_end, df.alt0)
df.pos = ifelse.(df.alt0 .> 599e3, df.r, df.r0)


include("setup.jl")
hmin = 80e3     #m
hmax = alt0+1e4 #m
intervals = 1e3 #m
hmsis = hmin:intervals:hmax

#Altitude histogram
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
save(joinpath("results", "hist_ionizations_height.png"), fig)


# 3d point plot of endpoints (earth centered coordinates)
fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1),)
sc = scatter!(ax,Point3.(df.pos), markersize = 2)#,
    #axis=(limits=(nothing, nothing, nothing),),)
#zlims!(ax, 5.99e6, 6.01e6)
#ylims!(ax, 2.465e6, 2.48e6)
#xlims!(ax, -1.5e4/2, 1.5e4/2)

display(fig)

# zoom in by using smaller dataset
filter_hmin = 150e3
filter_hmax = 151e3
df_filtered = filter(row -> row.alt > filter_hmin && row.alt < filter_hmax, df)

# meshscatter keeps dimensions const
fig, ax, ms = scatter(Point3.(df_filtered.pos), markersize = 1e1)
display(fig)

fig = Figure()
ax = Axis3(fig[1, 1])
scatter!(ax,Point3.(df_filtered.pos), markersize = 1e1)
scatter!(ax,Point3(sum(df_filtered.pos)./size(df_filtered.pos, 1)), markersize = 1e1)

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
[b/norm(b) for b in B0]

or = [0, 0, 0.0]
#fig = arrow(Point3(or), Vec3(B0[1]), alpha = 0.3)#, tipcolor = :red, label = "B/|B|")
#arrows!(Point3(or), Vec3(B0[end]))
#magnetic field cunrvature can be neglected

#create orthogonal vectorsystem along B0, using B0 on top: B0[end]
include("local_ortognal_basis.jl")
u1, u2, u3 = local_orthogonal_basis(B0[end])

origin = (c.re + hmsis[1]) * [0, cos(lat_gmag), sin(lat_gmag)]

#convert into magnetic field aligned coordinates, using magnetic field on top
mat = inv([u2 -u1 -u3])
df.pos_fa = [mat * (p-origin) for p in  df.pos]
df_filtered.pos_fa = [mat * (p-origin) for p in  df_filtered.pos]

# figure shows 400m drift in y
fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1))
sc = scatter!(ax, Point3.(df_filtered.pos_fa)./1e3, markersize = 1)#, 
    #axis=(limits=(nothing, nothing, nothing),),)
#zlims!(ax, -2e5, -1e5)
#ylims!(ax, -1.5e5, -0.5e5)
#ylims!(ax, -3, 3)
display(fig)

meshscatter(Point3.(df_filtered.pos_fa), markersize = 1)
#clearly slanted

#trace entire fiedline, for better correction:
#start with B0 at starting point:
lat_gmag = loc_gmag[1]
gc0 = (c.re + alt0) .* [0, cos(lat_gmag), sin(lat_gmag)]
B0 = dipole_field_earth(gc0)
#create orthogonal vectorsystem along B0
u1, u2, u3 = local_orthogonal_basis(B0)

function trace_fieldline(p0, fB)
    p = p0
    B = fB(p0)
    p_next = Vector(undef, Int(1e6))
    p_next[1] = p
    B_save = Vector(undef, Int(1e6))
    B_save[1] = B
    ii = 1
    while altitude(p) > 0
        ii = ii + 1
        p = p + B/norm(B)*1 #m steps
        B = fB(p)
        p_next[ii] = p
        B_save[ii] = B
    end
    return p_next[1:ii], B_save[1:ii]
end
trace, B = trace_fieldline(gc0, dipole_field_earth)
origin = trace[end]

#visualize trace in original coord
scatter(Point3.(trace))
altitude.(trace)


#df.nearest_trp = [findfirst(x -> x < altitude(p)-10000, altitude.(trace)) for p in df.pos]

#some functions to align every position along fieldine
function align_to_B(B0, pos, origin)
    u1, u2, u3 = local_orthogonal_basis(B0)
    mat = inv([u2 -u1 -u3])
    pos_fa = mat * (pos-origin)
    return pos_fa
end
 

df.pos
trace

df.nearest_idx = [argmin([norm(t - pos) for t in trace]) for pos in df.pos]

function project_on_line(pos, a, b)
    #projects a vector pos on the closest point on a line,
    #given by vecotrs a, b: line = a + c*b
    # b must be unit vector
    # with c a free parameter
    return a - (dot((a - pos), b) * b)
end

function find_fa_pos(pos, idx)
    a = trace[idx]
    b = B[idx]/norm(B[idx])
    nearest_pos_on_line = project_on_line(pos, a, b)
    B0 = dipole_field_earth(nearest_pos_on_line)
    return align_to_B(B0, pos, origin)
end

df.pos_fa = [find_fa_pos(p, i) for (p, i) in zip(df.pos, df.nearest_idx)]
df_filtered = filter(row -> row.alt > filter_hmin && row.alt < filter_hmax, df)
meshscatter(Point3.(df_filtered.pos_fa), markersize = 1e1)
display(current_figure())

fig = Figure()
ax = Axis3(fig[1, 1], aspect=(1, 1, 1))
sc = scatter!(ax, Point3.(df_filtered.pos_fa)./1e3, markersize = 1e1)#, 
    #axis=(limits=(nothing, nothing, nothing),),)
#zlims!(ax, -2e5, -1e5)
#ylims!(ax, -1.5e5, -0.5e5)
#ylims!(ax, -3, 3)

display(fig)
# slant reduced to 75m !!!
# what does that mean????


#trace electron to surface:
include("ode_boris_mover.jl")
include("espread.jl")
n_mfp = 1
E0 = 10000.0
lim_pitch = 20/180*pi
densityf(alt) = 0
r0, v0 = initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c)
status, r, v, t = ode_boris_mover_mfp(n_mfp, [r0; v0], -c.qe, c.me, dipole_field_earth!, cs_all_sum, densityf)#; OPS = [])

rp = [copy(c) for c in eachcol(r)]
r_filtered = filter(x -> altitude(x) > filter_hmin && altitude(x) < filter_hmax, rp)


nearest_idx = [argmin([norm(t - pos) for t in trace]) for pos in r_filtered]

B0 = [dipole_field_earth(pos) for pos in eachcol(r)]
r_fa = [align_to_B(B[nearest_idx[i]], r[:, i], origin) for i in axes(r, 2)]
sc = scatter!(ax, Point3.(r_fa)./1e3, markersize = 1e1)#, 
zlims!(ax, 0, 500)
ylims!(ax, 0, 3)
xlims!(ax, -0.1, 0.1)


altitude(r[:, end])r