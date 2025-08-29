





# for dipole field!!
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