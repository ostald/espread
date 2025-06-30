include("phase_functions.jl")
include("constants.jl")
include("rot_around_v.jl")
using DataInterpolations: LinearInterpolation

#scattering angle theta

function scatter_angle_el_inel(E0, phase_fcn, scatter_mode, nsample::Int)
    dA = 0.1 #degree
    nA = Int(180/dA + 1)
    angles_mean_deg = LinRange(0, 180, nA)
    angles_lim_deg = [angles_mean_deg .- dA/2; 180]
    angles_lim_deg[1] = 0
    angles_mean = angles_mean_deg/180*pi
    angles_lim = angles_lim_deg/180*pi

    ph_el, ph_in = vec.(phase_fcn(angles_mean, E0))
    #ph_el, ph_in = vec.(phase_fcn_O(angles_mean, E0))

    if scatter_mode == "elastic"
        phfcn = ph_el
    elseif scatter_mode == "inelastic"
        phfcn = ph_in
    end
    pdf_discrete = convert_phase_fcn_to_3D(phfcn, angles_lim)
    """
    lines(angles_mean_deg, convert_phase_fcn_to_3D(ph_el, angles_lim))
    lines!(angles_mean_deg, convert_phase_fcn_to_3D(ph_in, angles_lim))
    """
    cdf_discrete= [0; cumsum(pdf_discrete)] 
    cdf_disc_norm = cdf_discrete ./ cdf_discrete[end]
    """
    lines(angles_lim_deg, cdf_disc_norm)
    lines!(angles_mean_deg, pdf_discrete/cdf_discrete[end])
    """
    inverse_cdf_itp = LinearInterpolation(limits, cdf_disc_norm)
    if nsample != 1
        random_theta = inverse_cdf_itp(rand(nsample))
    else
        random_theta = inverse_cdf_itp(rand())
    end
    return random_theta

    """
    nsample = Int(1e7)
    hist(scatter_angle_el_inel(E0, phase_fcn_O, "elastic", nsample)/pi*180,
        bins = angles_lim_deg
        )
    lines!(angles_mean_deg, pdf_discrete/cdf_discrete[end]*nsample)
    display(current_figure())
    """
end

function scatter(v_in, E_out, theta, phi)
    u1, u2, u_par = local_orthogonal_basis(v_in)
    # use u1 as reference for first rotation
    # rotate scattering angle theta first
    v_r1 = rot_around_v(u1, theta)*v_in
    # rotate resulting vector around v_in, using random phase phi
    v_r2 = rot_around_v(v_in, phi)*v_r1
    # scale vector according to outgoing energy
    abs_vout = v_abs(E_out)
    v_out = v_r2/norm(v_r2) * abs_vout
    return v_out
end


"""
v_in = [0, 0, 1]
theta = pi/2
phi = pi
E_out = 2 * c.me / c.qe
v_out = scatter(v_in, E_out, theta, phi)
fig = arrows([Point3f([0, 0, 0])], [Vec3f(v_in)], arrowcolor = :red, label = "v_in")
arrows!([Point3f([0, 0, 0])], [Vec3f(v_r1)], arrowcolor = :blue)
arrows!([Point3f([0, 0, 0])], [Vec3f(v_r2)], arrowcolor = :green)
arrows!([Point3f([0, 0, 0])], [Vec3f(v_out)], arrowcolor = :green)

v_in = [0, 0, 1] * v_abs(500)
theta = pi/4
phi = pi/2
E_out = E_ev(norm(v_in)) *2 
v_out = scatter(v_in, E_out, theta, phi)
"""

function scatter_inelastic_N2(v_in, E_exc)
    E_in = E_ev(norm(v_in))
    theta = scatter_angle_el_inel(E_in, phase_fcn_N2, "inelastic", 1)
    phi = rand()*2*pi
    E_out = E_in - E_exc
    @assert E_out > 0
    v_out = scatter(v_in, E_out, theta, phi)
    return v_out
end

function scatter_elastic_N2(v_in, E_exc = 0)
    E_in = E_ev(norm(v_in))
    theta = scatter_angle_el_inel(E_in, phase_fcn_N2, "elastic", 1)
    println(string(theta/pi*180))
    phi = rand()*2*pi
    E_out = E_in
    @assert E_out > 0
    v_out = scatter(v_in, E_out, theta, phi)
    return v_out
end

"""
v_in = [0, 0, 1] * v_abs(500)
v_out = scatter_elastic_N2(v_in)
fig = arrows([Point3f([0, 0, 0])], [Vec3f(v_in/v_abs(500))], arrowcolor = :red, label = "v_in")
arrows!([Point3f([0, 0, 0])], [Vec3f(v_out/v_abs(500))], arrowcolor = :green)
"""

function scatter_ion_N2(v_in, E_ionization)
    E_primary = E_ev(norm(v_in))
    Es_out = E_secondary_e_N2(E_primary, E_ionization)
    Ep_out = E_primary - E_sec - E_ionization
    if E_sec > E_primary
        theta = pi/4
        phi = rand()*2*pi
        vp_out = scatter(v_in, Ep_out, theta, phi)

        vs_out = scatter(v_in, Es_out, theta, pi-phi)
    else
        #scatter isotrop:
        vp_out = v_in/nprm(v_in) * v_abs(Ep_out)
        vs_out = ???
    end
    #record position, velocity vector e_secondary
    return vp_out, vs_out
end


function scatter_inelastic_O2()
    return 0
end
function scatter_elastic_O2()
    return 0
end
function scatter_ion_O2()
    return 0

end
function scatter_inelastic_O()
    return 0
end
function scatter_elastic_O()
    return 0
end
function scatter_ion_O()
    return 0
end





"""
unnormalized_pdf = DCSN2(angles_mean_rad, E0)
normalized_pdf = unnormalized_pdf/sum(unnormalized_pdf)
lines(angles_mean, normalized_pdf)
lines(angles_mean, cdf_discrete)
lines(cdf_discrete, angles_mean)
lines!(cdf_discrete[2:end] - cdf_discrete[1:end-1], angles_lim[2:end-1])
"""

"""
using DataInterpolations: LinearInterpolation
include("phase_functions.jl")

# resolution at small angles actually might have a big influence in the spreading!
# resolution of 0.1deg seems to be good, histogram flatten out at low angles.
dA = 0.1 #degree
nA = Int(180/dA + 1)

# forward scattering cone: limits (-0.5 ... 0.5 deg) with mean 0 deg
# but wo only want the scattering angle from 0 to 180 deg
# therefore limits = [0,    0.5,    1.5,    2.5, ....]
# and          means = [  0,      1,      2,     3, ...]
# such that the forward cone in also 1 degree wide

angles_mean_deg = LinRange(0, 180, nA)
angles_lim_deg = [angles_mean_deg .- dA/2; 180]
angles_lim_deg[1] = 0
angles_mean = angles_mean_deg/180*pi
angles_lim = angles_lim_deg/180*pi

E0 = 500

ph_el, ph_in = vec.(phase_fcn_O(angles_mean, E0))
ph_el_normalized = convert_phase_fcn_to_3D(ph_el, angles_lim)
ph_in_normalized = convert_phase_fcn_to_3D(ph_in, angles_lim)

cdf_discrete = cumsum(ph_el_normalized)#/sum(unnormalized_pdf)

lines(angles_mean, ph_el_normalized)
lines!(angles_mean, ph_in_normalized)

#produces descrete angles, as determined in angles_mean.
# instead, interpolate!
nsample = Int(1e6)
hist([angles_mean[findfirst(cdf_discrete .> r)] for r in rand(nsample)],
    bins = angles_lim
    )
lines!(angles_mean, ph_el_normalized*nsample)
display(current_figure())


# correct interpolation of cdf:
# 0. Define Intervals for cdf and mean values for pdf
limits = angles_lim
mean = angles_mean
# 1. create discrete pdf
pdf_disc = convert_phase_fcn_to_3D(ph_el, angles_lim)
# 2. creat discrete cdf, add 0, normalize
cdf_disc = [0; cumsum(pdf_disc)] 
cdf_disc = cdf_disc ./ cdf_disc[end]
lines(limits, cdf_disc)
lines!(mean, pdf_disc)
# 3. linearly interpolate
inverse_cdf_itp = LinearInterpolation(limits, cdf_disc)
# 4. uniformly distributed random number as argument in itpf produces random value sample from original pdf
r_uniform = rand()
r_pdf = inverse_cdf_itp(r_uniform)

rand_test = 0:0.0001:1
lines!(inverse_cdf_itp(rand_test), rand_test)


function sample_pdf(pdf_disc, limits, nsamples)
    # produces one sample from an arbitrary pdf
    # pdf has to be discrete, i.e. evaluated at n points
    # which are between n+1 limits (usually the mean values)
    cdf_disc = [0; cumsum(pdf_disc)] 
    cdf_disc = cdf_disc ./ cdf_disc[end]    
    inverse_cdf_itp = LinearInterpolation(limits, cdf_disc)
    return inverse_cdf_itp(rand(nsamples))
end

#apply sampling with interpolation to all random numbers!
#instead of returning a sampled value, return sampler, i.e. interpolated function
"""