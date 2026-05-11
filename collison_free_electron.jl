#collison free electron should follow magnetic field line i.e. gyrocenter.

include("util.jl")
include("ode_boris_mover.jl")
include("setup.jl")
b_model = "dipole"
lim_pitch_deg = 20
pitch_angle_distribution = "single_angle"
E0 = 4e3 #ev
name = "r15_test_dipole_"
res_dir = joinpath("results", name * string(now()))
mkdir(res_dir)
save_commit_hash(res_dir)

include("espread.jl")


    Bin!(B, p) = 0
    if b_model == "dipole"
        Bin!(B, p) = dipole_field_earth!(B, p)
    elseif b_model == "vertical"
        Bin!(B, p) = convergent_vertical_field!(B, p)
    end

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
    hintervals = 1e2 #m
    densityf = make_densityf(hmin, hmax, hintervals, loc_geod)
    #stack(densityf_fast.(hmsis))' == atm
    #    > true

    #hmsis = 80e3:1e3: #km
    #densityf = atmospheric_model([[2020, 12, 12, 18, 0, 0]], hmsis, loc_geod[1], loc_geod[2])

    # Results directory    
    
    res_file = joinpath(res_dir, "res_$(E0)eV_$(lim_pitch_deg)deg.bin")
    io = open(res_file, "w")
        serialize(io,
        [E0,
        lim_pitch_deg,
        seed_value,
        hmin,
        hmax,
        hintervals])


    println("res_file = ", res_file)
    
    lim_pitch = lim_pitch_deg/180*pi

    n_e_sim = 1

    #partricle generation (1 for primaries, 2 for secondaries etc)
    generation = 1
    idx_scatter_rec = -1

    densityf = function f(x) return 0 end
    r0, v0 = initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c, b_model, nPerGyro, Bin!, densityf, pitch_angle_distribution)
    
    v = v0
    r = r0
    E = E_ev(norm(v))
    status = -1 #undef
    #status -1 will be retained if particle has low energy, such that the boris mover is not initiated

    #propagate electron until it runs out of energy to ionize
    
    n_mfp = rand(Exponential())
    status, r, v, t = ode_boris_mover_mfp(n_mfp, r, v, -c.qe, c.me, Bin!, cs_all_sum, densityf, trace = true, nPerGyro = nPerGyro)
    
    #propagate_electron(v0, r0, idx_scatter_rec, densityf, io, c, Bin!, nPerGyro, generation)
