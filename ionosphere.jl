using PCHIPInterpolation: Interpolator


function make_ionosphere()
    # returns interpolators, call ne_itp.(height_array) to get interpolated electron densities

    ne = [0,
        0.0,
        1.9539662964982204e9,
        5.686028841257099e9,
        2.991124126010653e10,
        1.6810931195392596e11,
        2.878097675250895e11,
        1.040554770249579e11,
        9.211827022450801e10,
        2.3801360801107352e11,
        0.0]

    h_ne  = [0,
         60000.0,
         77494.78897351192,
         81351.22576432955,
        104754.0075795146,
        141588.59851676473,
        185697.04248022253,
        282302.0212821382,
        463218.5961653356,
        544078.1414334535,
        800000.0]

    ne_itp = Interpolator(h_ne, ne, extrapolate = true)

    Te = [
         222.80432947543585,
         222.80432947543585,
         225.38518690389742,
         430.15397867483375,
        1129.8497935860632,
        2193.229997821504,
        2876.5584484149776,
        3276.0396689989534,
        3276.0396689989534]

    h_Te = [
         60000.0,
         77473.28920354624,
         81328.08973538037,
        113477.31619360304,
        161562.25281756447,
        198773.22237450292,
        301079.2038587333,
        568769.3511646576,
        800000.0]


    Te_itp = Interpolator(h_Te, Te, extrapolate = true)
    
    return ne_itp, Te_itp
end


# IRI-2026 via IRI/libiri.so, built with:
#   cd IRI && gfortran -shared -fPIC -O2 -std=legacy -w -o libiri.so irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for rocdrift.for
# IRI/ig_rz.dat and IRI/apf107.dat must be present (and current), see IRI/00readme.txt.
const IRI_DIR = joinpath(@__DIR__, "IRI")
const LIBIRI = joinpath(IRI_DIR, "libiri.so")
const IRI_INDICES_READ = Ref(false)

function iri_jf_default()
    # IRI-recommended switches (see iritest.for)
    jf = ones(Int32, 50)
    jf[[4, 5, 6, 23, 30, 33, 35, 39, 40, 47]] .= 0
    jf[34] = 0  # messages off
    return jf
end

"""
    iri(time, lat, lon, heights; jf = iri_jf_default())

Run IRI_SUB for one location and time.
time: [year, month, day, hour, minute, second] in UT
lat, lon: geographic latitude/longitude [rad]
heights: altitudes [m], equally spaced (a range)
Returns (ne [m^-3], Te [K], outf) with outf the full IRI OUTF array (20 x length(heights)).
IRI marks values it does not compute with -1.
"""
function iri(time, lat, lon, heights::AbstractRange; jf = iri_jf_default())
    year, month, day, hour, minute, second = time
    mmdd = Int32(100 * month + day)
    dhour = Float32(hour + minute / 60 + second / 3600 + 25)  # +25 means UT
    nh = length(heights)
    outf = zeros(Float32, 20, nh)
    oarr = zeros(Float32, 100)
    out_chunk = zeros(Float32, 20, 1000)

    # Fortran opens its coefficient and index files by relative path
    cd(IRI_DIR) do
        if !IRI_INDICES_READ[]
            ccall((:read_ig_rz_, LIBIRI), Cvoid, ())
            ccall((:readapf107_, LIBIRI), Cvoid, ())
            IRI_INDICES_READ[] = true
        end
        # IRI_SUB returns at most 1000 heights per call
        for i0 in 1:1000:nh
            i1 = min(i0 + 999, nh)
            oarr .= 0
            ccall((:iri_sub_, LIBIRI), Cvoid,
                (Ptr{Int32}, Ref{Int32}, Ref{Float32}, Ref{Float32}, Ref{Int32}, Ref{Int32}, Ref{Float32},
                 Ref{Float32}, Ref{Float32}, Ref{Float32}, Ptr{Float32}, Ptr{Float32}),
                jf, Int32(0), Float32(rad2deg(lat)), Float32(rad2deg(lon)), Int32(year), mmdd, dhour,
                Float32(heights[i0] / 1e3), Float32(heights[i1] / 1e3), Float32(step(heights) / 1e3),
                out_chunk, oarr)
            outf[:, i0:i1] .= out_chunk[:, 1:(i1 - i0 + 1)]
        end
    end
    return Float64.(outf[1, :]), Float64.(outf[4, :]), outf
end

"""
    make_ionosphere_iri(time, loc_geod; hmin = 60e3, hmax = 1000e3, hstep = 1e3)

Same interface as make_ionosphere(), but ne and Te profiles come from IRI.
Heights where IRI does not compute ne are set to 0.
"""
function make_ionosphere_iri(time, loc_geod; hmin = 60e3, hmax = 1000e3, hstep = 1e3)
    hh = hmin:hstep:hmax
    ne, Te, _ = iri(time, loc_geod[1], loc_geod[2], hh)
    ne = max.(ne, 0.0)
    ne_itp = Interpolator(collect(hh), ne, extrapolate = true)
    Te_itp = Interpolator(collect(hh), Te, extrapolate = true)
    return ne_itp, Te_itp
end

#=
using WGLMakie
ne_itp, Te_itp = make_ionosphere()
hh = (1:700e3)
ne = ne_itp.(hh)
Te = Te_itp.(hh)

fig = lines(ne, hh,
     axis = (#xscale = log10,
        limits = ((0, 1e2), nothing),),
    )
#display(fig)

lines(Te, hh)
=#


# Example usage of the IRI functions
#=
    using CairoMakie
    CairoMakie.activate!()

    loc_geod = deg2rad.([69.58, 19.23])  # Tromsø, as in setup.jl
    hh = 80e3:1e3:800e3
    hkm = hh ./ 1e3

    # 1) profiles from make_ionosphere_iri (interpolators) vs. the hardcoded make_ionosphere
    cases = [
        ("IRI 2020-12-12 18 UT (winter night)", [2020, 12, 12, 18, 0, 0]),
        ("IRI 2020-06-12 11 UT (summer noon)", [2020, 6, 12, 11, 0, 0]),
        ("IRI 2024-12-12 18 UT (winter night, solar max)", [2024, 12, 12, 18, 0, 0]),
    ]
    fig1 = Figure(size = (1000, 600))
    ax_ne = Axis(fig1[1, 1], xscale = log10, xlabel = "nₑ [m⁻³]", ylabel = "altitude [km]",
        title = "electron density, Tromsø")
    ax_Te = Axis(fig1[1, 2], xlabel = "Tₑ [K]", title = "electron temperature")
    linkyaxes!(ax_ne, ax_Te)
    for (lbl, t) in cases
        ne_itp, Te_itp = make_ionosphere_iri(t, loc_geod)
        lines!(ax_ne, max.(ne_itp.(hh), 1e6), hkm, label = lbl)
        lines!(ax_Te, Te_itp.(hh), hkm, label = lbl)
    end
    ne_itp, Te_itp = make_ionosphere()
    lines!(ax_ne, max.(ne_itp.(hh), 1e6), hkm, color = :black, linestyle = :dash, label = "make_ionosphere (hardcoded)")
    lines!(ax_Te, Te_itp.(hh), hkm, color = :black, linestyle = :dash, label = "make_ionosphere (hardcoded)")
    xlims!(ax_ne, 1e8, 1e12)
    Legend(fig1[2, 1:2], ax_ne, orientation = :horizontal, nbanks = 2, framevisible = false)
    save(joinpath(@__DIR__, "figures", "iri_profiles.png"), fig1)
    display(fig1)

    # 2) daily variation from iri() directly: ne(altitude, UT)
    uts = 0:0.5:24
    ne_day = stack([iri([2020, 12, 12, floor(Int, ut), round(Int, 60 * (ut % 1)), 0],
                        loc_geod[1], loc_geod[2], hh)[1] for ut in uts])  # (height, time)
    ne_day[ne_day .<= 0] .= NaN  # -1 = not computed by IRI
    fig2 = Figure(size = (900, 500))
    ax = Axis(fig2[1, 1], xlabel = "UT [h]", ylabel = "altitude [km]",
        title = "IRI log₁₀ nₑ [m⁻³], Tromsø, 2020-12-12", xticks = 0:3:24)
    hm = heatmap!(ax, uts, hkm, log10.(ne_day'), colormap = :viridis, colorrange = (8, 11.5))
    Colorbar(fig2[1, 2], hm)
    save(joinpath(@__DIR__, "figures", "iri_daily_ne.png"), fig2)
    display(fig2)
=#