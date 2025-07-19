using Pkg
Pkg.activate(".")

using Distributions
using LinearAlgebra
using Revise 
using Random

include("energy_secondary_e.jl")
include("magnetic_field.jl")
include("constants.jl")
include("local_ortognal_basis.jl")
include("ode_boris_mover.jl")
#load cross sections
include("cross_sections.jl")
include("get_msis.jl")
include("setup.jl")


#using GLMakie


function initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c)
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

    return r0, v0
end
    

function save_endpoint(res_dir, r0, v0, status, r, v) #save status too?
    open(joinpath(res_dir, "results.txt"), "a") do file
        write(file, "$r0\t$v0\t$status\t$r\t$v\n")
    end
end


##
function propagate_electron(v0, r0, densityf, res_dir, c)
    #list of secondary electrons
    secondary_e = []
    
    v = v0
    r = r0
    E = E_ev(norm(v))
    status = -1 #undef
    #status -1 will be retained if particle has low energy, such that the boris mover is not initiated

    #propagate electron until it runs out of energy to ionize
    while E > 12.072
        r0v0 = [r; v]
        #magnetic field handle for boris mover: (static field)
        Bin!(B, p) = dipole_field_earth!(B, p)
        # sample number of mean free paths travelled:
        n_mfp = rand(Exponential())
        status, r, v, t = ode_boris_mover_mfp(n_mfp, r0v0, -c.qe, c.me, Bin!, cs_all_sum, densityf)
        E = E_ev(norm(v))

        if status == 1
            nothing
        elseif status == 2
            #record? do something?
            println("Particle lost, recorded.\n")
            save_endpoint(res_dir, r0, v0, status, r, v) #save status too?
            break
            #alternatively, check in while statement if state != 2
        elseif status == 0
            error("Boris mover failed. Investigate!")
        else 
            pritnln("status :", status)
            error("This should not happen. Investigate!")
        end
        
        #plot_local_v(r, v, Bin)

        ##
        # after moving electron, scatter:

        # use cross sections of all scatter processes directly, not summed for each species
        # thereby one random number decides directly which scattering process is triggered 
        #try
        ns = densityf(altitude(r))
        #catch
        #    print(r)
        #    error()
        #end       
        scatter_p = vcat((cs_all(E) .* ns)...)
        idx_scatter = findfirst(cumsum(scatter_p) .> rand() * sum(scatter_p))

        """
        print("           Scattering Process = ", sp_all[idx_scatter, 1], 
            "\n                       E_loss = ", sp_all[idx_scatter, 2], 
            "\n                Ions Produced = ", sp_all[idx_scatter, 3], 
            #"\nCross Section function handle = ", sp_all[idx_scatter, 4], 
            "\n           Scattering Partner = ", sp_all[idx_scatter, 5], 
            "\n   Scattering function handle = ", sp_all[idx_scatter, 6],
            "\n\n")
        """

        # all scatterign functions have the same input:
        # v_in, E_loss
        # get paramters for scattering:
        E_loss = sp_all[idx_scatter, 2]
        sc_f   = sp_all[idx_scatter, 6]

        #ensure E_loss is not bigger than the kinetic energy
        while E_loss > E_ev(norm(v))
            idx_scatter = findfirst(cumsum(scatter_p) .> rand() * sum(scatter_p))
            E_loss = sp_all[idx_scatter, 2]
            sc_f   = sp_all[idx_scatter, 6]
        end

        # scatter:
        out = 0
        try
            out = sc_f(v, E_loss)
        catch
            print("r = ", r,"; v = ", v, "; idx_scatter = ", idx_scatter, "\n")
            print("Ekin - E_loss = ", E_ev(norm(v)) - E_loss)
            error("boom")
        end
        """
        # testing all scattering functions:
        for idx_scatter in 1:63
            E_loss = sp_all[idx_scatter, 2]
            sc_f = sp_all[idx_scatter, 6]
            out = sc_f(v, E_loss)
            println(idx_scatter, out)
            println(" ") 
        end
        """

        if sp_all[idx_scatter, 3] == 0 #non-ionising collisions
        #if typeof(out) == Vector{Float64}
            v = out
        else #elseif sp_all[idx_scatter, 3] == 1 #ionising collisions
            vp_out, vs_out = out
            v = vp_out
            secondary_e = record_secondary(r, vs_out, secondary_e)
        #else #double-ionising collisions not permitted as of now
        #    vp_out, vs1_out, vs2_out = out
        #    v = vp_out
        #    secondary_e = record_secondary(r, vs1_out, secondary_e)
        #    secondary_e = record_secondary(r, vs2_out, secondary_e)
        end

        #update energy after collision!
        E = E_ev(norm(v))

    end

    # save end point and velocity of parent electron
    save_endpoint(res_dir, r0, v0, status, r, v) #save status too?

    # go through secondary electrons
    for se in secondary_e
        r0 = se[1]
        v0 = se[2]
        propagate_electron(v0, r0, densityf, res_dir, c)
    end
end



function main(E0, N_electrons, alt0, lim_pitch_deg, loc_gmag, loc_geod, c)

    # make sure all values are floats
    E0 = Float64(E0)
    alt0 = Float64(alt0)
    lim_pitch_deg = Float64(lim_pitch_deg)
    loc_gmag = Float64.(loc_gmag)
    loc_geod = Float64.(loc_geod)
    
    #seed_value = round(Int, E0 + lim_pitch_deg)
    seed_value = rand(Int)
    Random.seed!(seed_value)
    
    hmin = 80e3     #m
    hmax = alt0+1e4 #m
    hintervals = 1e3 #m
    densityf = make_densityf(hmin, hmax, hintervals)
    #stack(densityf_fast.(hmsis))' == atm
    #    > true

    #hmsis = 80e3:1e3: #km
    #densityf = atmospheric_model([[2020, 12, 12, 18, 0, 0]], hmsis, loc_geod[1], loc_geod[2])


    # Results directory
    res_dir = joinpath("results", "res_$(E0)eV_$(lim_pitch_deg)deg_" * string(now()))
    mkdir(res_dir)
    open(joinpath(res_dir, "results.txt"), "w") do file
        write(file, "E0 = $E0\n")
        write(file, "lim_pitch_deg = $lim_pitch_deg\n")
        write(file, "seed_value = $seed_value\n")
        write(file, "hmin = $hmin\n")
        write(file, "hmax = $hmax\n")
        write(file, "hintervals = $hintervals\n\n")
    end
    
    lim_pitch = lim_pitch_deg/180*pi

    n_e_sim = 1
    
    while n_e_sim <= N_electrons
        println("Electron number: ", n_e_sim)
        #try
            r0, v0 = initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c)
            propagate_electron(v0, r0, densityf, res_dir, c)
            n_e_sim = n_e_sim +1
        #catch
        #    println("re-initiating electron")
        #end
    end
    return nothing
end

#main(E0, 10, alt0, pitch_lim, loc_gmag, loc_geod, c)


#E0 = 4000 #eV
#@time main(Int(E0), 10, alt0, lim_pitch_deg, loc_gmag, loc_geod, c)
"""
using Distributed
procs = addprocs(2)

@everywhere procs begin Base.MainInclude.eval(include("espread.jl")) end
f = @spawnat 2  main(E0, 2, alt0, lim_pitch, loc_gmag, loc_geod)
ff = fetch(f)


for E0 in e_energy
    for lim_pitch_deg in pitch_limits
        lim_pitch = lim_pitch_deg/180*pi
        main(E0, N_electrons, alt0, lim_pitch, res_dir, loc_gmag, loc_geod)
    end
end
"""


"""
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
"""


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
"""

"""
# cumulative sum of cross section * density to decide scattering partner:
cross_sections = cs_all_sum(E)
ns = densityf(altitude(r))
# sample uniformly between scattering partners:
r_scatter = rand()*sum(cross_sections .* ns)
#species index: according to densityf [N2, O2, O]
idx_species = findfirst(cumsum(cross_sections .* ns) .> r_scatter)
"""


"""
nsample = Int(1e7)
hist([findfirst(cumsum(scatter_p) .> r_scatter) for r_scatter in rand(nsample).*sum(scatter_p)],
    bins = 0.5:1:63.5,
    scale_to=scatter_p[1]
    )
"""


"""

using HDF5
res_dir = joinpath("results")
num_results = floor(Int, N_electrons * E0 / 35 * 1.1 ) #approximate number of electrons
mkdir(res_dir)
h5open(joinpath(res_dir, "results.h5"), "w") do file
    file["E0"] = E0
    file["lim_pitch_deg"] = pitch_limits_deg
    file["seed_value"] = seed_value
    file["results"] = zeros(0, num_results)  # Preallocate an array of size 10
end

data = ["lala", [1, 2, 3], 0.1]

function save_data(res_dir, data)
    i = 2
    h5open(joinpath(res_dir, "results.h5"), "r+") do file
        file["results"][i] = data  # Update the i-th element
    end
    nothing
end
save_data(res_dir, data)


function load_data(res_dir)
    return load(joinpath(res_dir, "results.jld2"))
end

d = load_data(res_dir)
d["results"]

"""