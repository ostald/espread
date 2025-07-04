using Distributions
using LinearAlgebra
using Revise 
using Random

include("magnetic_field.jl")
include("constants.jl")
include("local_ortognal_basis.jl")
include("ode_boris_mover.jl")
#load cross sections
include("cross_sections.jl")
include("get_msis.jl")
include("setup.jl")


#using GLMakie

n_e_sim = 0
while n_e_sim < N_electrons
    ##
    # initialize new electron
    
    # gyrocenter 
    lat_gmag = loc_gmag[1]
    gc0 = (c.re + alt0) .* [0, cos(lat_gmag), sin(lat_gmag)]
    #check altidute
    @assert abs(alt0 - altitude(gc0)) < 10 #m        10m within target altitude is ok.
    
    B0 = dipole_field_earth(gc0)
    #create orthogonal vectorsystem along B0
    u1, u2, u3 = local_orthogonal_basis(B0)

    pitch = rand()*lim_pitch    # random pitch angle within limits
    phase = rand()*2*pi         # random phase angle

    # starting velocity
    # energy must be float
    E0 = float(E0)
    v0_mag = v_abs(E0)           # convert to veloicty => relativistic?
    v0_par = v0_mag*cos(pitch)   # calculate parallel component

    # for parallel component, find normal vecotors, get random phase, distribute velocity accordingly
    # velocity vector starts a gyroradius away from gyrocenter
    v0_perp = sqrt(v0_mag^2 - v0_par^2)

    r0_gyro = c.me * v0_perp / (c.qe * norm(B0))

    # calculate r0 and v0 in n1, n2 directions, and convert to original coordinate system
    r_n1 =  cos(phase) * r0_gyro .* u1
    r_n2 =  sin(phase) * r0_gyro .* u2
    v_n1 = -sin(phase) * v0_perp .* u1
    v_n2 =  cos(phase) * v0_perp .* u2
    v_n3 =  v0_par .* u3
    r0 = gc0 .+ r_n1 .+ r_n2
    v0 = v_n1 .+ v_n2 .+ v_n3

    # alternatively, the offset from the gyrocenter can be calculated using
    # the cross product of v0, B0:
    #r0_2 = gc0 .- c.me * cross(v0, B0) / (c.qe * norm(B0)^2)
    #r0_2 ≈ r0
    #   > true
    #r0 ./ r0_2
    #   > 3-element Vector{Float64}:
    #   >   1.0000000000000002
    #   >   1.0
    #   >   1.0

    @assert abs(alt0 - altitude(r0)) < r0_gyro * 1.1 #m   # must be 1.1*gyroradius within target altitude.

    # this could be a function, eg. initialize_electron(E0, pitch_lim, lat_gmag, alt0)
    ##

    propagate_electron(v0, r0)
   
    ##
    n_e_sim = n_e_sim +1
end


##
function propagate_electron(v0, r0)
    #list of secondary electrons
    secondary_e = []

    #propagate electron until it runs out of energy to ionize
    while E_ev(norm(v0)) > 12.072
        error("notimplemented")
        #move and scatter
    end

    # save end point and velocity of parent electron
    save_result(res_dir, r, v)

    # go through secondary electrons
    for se in secondary_e
        r0 = se[1]
        v0 = se[2]
        initialize_electron(v0, r0)
    end
end



# r0 is by definition a vector away from the center of earth, i.e. pointing outwards
# projecting B0 on r0 gives us the orientation of B0 (negative for earthwards, i.e. down)
B0_proj_r0 = dot(B0, r0)
# projecting v0 on r0 gives tells us whether the particle is flying up or down (negative for down)
# but that can be depending on the position of the particle along the gyroradius
v0_proj_r0 = dot(v0, r0)
# therefore, first calculate v_par = v proj. on B0,
# and determine direction of B0 by projecting it on r
# if the product is negative, v_par is pointing down
v0_proj_B0 = dot(v0, B0)
v_par_vec = v0_proj_B0*B0 / norm(B0)^2
B0_proj_r0 = dot(B0, r0)
v0_proj_B0 * B0_proj_r0

# plotting start positions, velocity and magnetic field

or = [0, 0, 0]
pp = r0-gc0

fig = arrows3d(Point3(or), Vec3(u3), tipcolor = :red, label = "B/|B|")
arrows3d!(Point3.([or, or]), Vec3.([u1, u2]))

arrows3d!(Point3(or), Vec3((r_n1 .+ r_n2)), tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
arrows3d!(Point3(pp), Vec3((v_n1 .+ v_n2)./1e7), tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
arrows3d!(Point3(pp), Vec3(v_n3./1e7), tipcolor = :red, tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
arrows3d!(Point3(pp), Vec3(v0)./1e7, tipcolor = :blue, label = "v", tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
axislegend()


r0v0 = [r0; v0]
#magnetic field handle for boris mover: (static field)
Bin(p) = dipole_field_earth(p)

#mean free path calculations
#msis altitude resolution could be coarser??
hmsis = 80e3:1e3:700e3 #km
densityf = atmospheric_model([[2020, 12, 12, 18, 0, 0]], hmsis, loc_geod[1], loc_geod[2])


# sample number of mean free paths travelled:
n_mfp = rand(Exponential(1))
status, r, v, t = ode_boris_mover_mfp(n_mfp, r0v0, -c.qe, c.me, Bin, cs_all_sum, densityf)

if status == 1
    nothing
elseif status == 2
    #record? do something?
    println("Particle lost. should we record that?")
elseif status == 0
    error("Boris mover failed. Investigate!")
end


"""
OPS = (trace = true,)
status, r, v, t = ode_boris_mover_mfp(n_mfp, r0v0, -c.qe, c.me, Bin, cs_all_sum, densityf, OPS=OPS)
rp = r .- r0
zs = LinRange(0, 3, 200)
meshscatter([tuple(p...) for p in eachcol(rp[:, 1321001:1321200])], markersize = 1, color = zs)


tt = rvt[7, :]
altitude = [norm(p) - c.re for p in eachcol(r)]
lines(t, altitude)
"""

function plot_local_v(r, v, Bin)
    B0 = Bin(r)
    u1, u2, u3 = local_orthogonal_basis(B0)
    v_par = dot(v, u3) * u3
    v_u1  = dot(v, u1) * u1
    v_u2  = dot(v, u2) * u2
    v_perp = v_u1 .+ v_u2
    r_gyro = c.me * cross(v_perp, B0) / (c.qe * norm(B0)^2)
    
    or = [0, 0, 0] #r + r_gyro
    pp = -r_gyro

    fig = arrows3d(Point3(or), Vec3(u3), tipcolor = :red, label = "B/|B|")
    arrows3d!(Point3.([or, or]), Vec3.([u1, u2]))

    arrows3d!(Point3(or), Vec3(-r_gyro), tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
    arrows3d!(Point3(pp), Vec3((v_u1 .+ v_u2)./1e7), tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
    arrows3d!(Point3(pp), Vec3(v_par./1e7), tipcolor = :red, tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
    arrows3d!(Point3(pp), Vec3(v)./1e7, tipcolor = :blue, label = "v", tipradius = 0.034, tiplength = 0.1, shaftradius = 0.015)
    axislegend()
    display(current_figure())
end

plot_local_v(r, v, Bin)

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
r_scatter = rand()
idx_scatter = findfirst(cumsum(scatter_p) .> r_scatter * sum(scatter_p))

sp_all[idx_scatter, :]
E_loss = sp_all[idx_scatter, 2]
sc_f = sp_all[idx_scatter, 6]
v

for idx_scatter in 1:36
    E_loss = sp_all[idx_scatter, 2]
    sc_f = sp_all[idx_scatter, 6]
    out = sc_f(v, E_loss)
    println(idx_scatter, out)
    println(" ") 
end


out = sc_f(v, E_loss)
if typeof(out) == Vector{Float64}
    v = out
else
    vp_out, vs_out = out
    v = vp_out
    secondary_e = record_secondary(r, vs_out, secondary_e)
end



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
include("energy_secondary_e.jl")
# there is a lower bound on the E_secondary, given by the resolution (usually 0.1 eV)
# could be resolved by first selectong E_secondary bin, and assume an even distribution within the bin,
# and select the exact E_secondary by second random number.

#do collision:
#1. decide which particle: sum densities, normalize,

