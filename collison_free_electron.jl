#collison free electron should follow magnetic field line i.e. gyrocenter.

include("util.jl")
include("ode_boris_mover.jl")
include("setup.jl")
include("ionosphere.jl")

b_model = "dipole"
lim_pitch_deg = 20
pitch_angle_distribution = "single_angle"
E0 = 25.0 #ev
name = "r20_ee_loss_comparison_"
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

    ne_int = function null(x) return 0 end
    Te_int = function null(x) return 0 end

    r0, v0 = initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c, b_model, nPerGyro, Bin!, densityf, ne_int, Te_int, pitch_angle_distribution)
    #r0, v0 = initialize_primary_electron(E0, loc_gmag, alt0, lim_pitch, c, b_model, nPerGyro, Bin!, densityf, pitch_angle_distribution)
    
    # no neutral collisions: lam = 0, so the mover runs until the electron hits the ground or escapes
    densityf = function null(x) return 0 end
    null_f(x) = 0.0  # ne = 0 disables e-e energy loss
    ne_on, Te_on = make_ionosphere()

    st_on, r_on, v_on, t_on = ode_boris_mover_mfp(Inf, r0, v0, -c.qe, c.me, Bin!, cs_all_sum, densityf, ne_on, Te_on, trace = true, nPerGyro = nPerGyro)
    st_off, r_off, v_off, t_off = ode_boris_mover_mfp(Inf, r0, v0, -c.qe, c.me, Bin!, cs_all_sum, densityf, null_f, null_f, trace = true, nPerGyro = nPerGyro)

    alt_on = altitude.(collect.(eachcol(r_on)))
    alt_off = altitude.(collect.(eachcol(r_off)))
    E_on = E_ev.(norm.(eachcol(v_on)))
    E_off = E_ev.(norm.(eachcol(v_off)))

    for (lbl, st, alt, E, t) in (("with e-e loss", st_on, alt_on, E_on, t_on), ("without e-e loss", st_off, alt_off, E_off, t_off))
        println("--- ", lbl, " ---")
        println("  status          = ", st)
        println("  steps           = ", length(t))
        println("  duration [s]    = ", t[end])
        println("  min altitude [km] = ", minimum(alt)/1e3)
        println("  E start [eV]    = ", E[1])
        println("  E end [eV]      = ", E[end])
        println("  rel. dE         = ", (E[end] - E[1])/E[1])
    end
    println("E_end(off) - E_end(on) [eV] = ", E_off[end] - E_on[end])

    serialize(io, (; st_on, r_on, v_on, t_on, st_off, r_off, v_off, t_off))
    close(io)

    include("util.jl"); include("constants.jl")
    E0, lim_pitch_deg, seed_value, hmin, hmax, hintervals, tr..., alt_on, alt_off, E_on, E_off = load_ee_comparison(res_file)
    st_on, r_on, v_on, t_on, st_off, r_off, v_off, t_off = tr
    #lines(d.t_on, d.E_on)

    using WGLMakie
#    CairoMakie.activate!()
    fig = Figure(size = (900, 700))
    ax1 = Axis(fig[1, 1], xlabel = "time [s]", ylabel = "altitude [km]",
        title = "$(E0) eV, $(lim_pitch_deg)° pitch")
    lines!(ax1, t_on, alt_on ./ 1e3, label = "with e-e loss")
    lines!(ax1, t_off, alt_off ./ 1e3, label = "without e-e loss", linestyle = :dash)
    axislegend(ax1)
    ax2 = Axis(fig[2, 1], xlabel = "time [s]", ylabel = "energy [eV]")
    lines!(ax2, t_on, E_on, label = "with e-e loss")
    lines!(ax2, t_off, E_off, label = "without e-e loss", linestyle = :dash)
    axislegend(ax2, position = :lb)
    fig_file = joinpath(res_dir, "ee_loss_comparison.png")
    save(fig_file, fig)
    println("figure saved: ", fig_file)
