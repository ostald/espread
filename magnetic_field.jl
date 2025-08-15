using LinearAlgebra
using StaticArrays
include("constants.jl")


function make_convergent_vertical_field_demo(c)
    z0 = (c.re + 80e3) #z0 at 80km above earth radius
    c1 = 8.22e15 #Tm3      c1 = mu_0 m / 4 pi ≈ m * 1e-7
    dBdz = -6*c1 / z0^4 * 2e3                               #factor 2e3
    function convergent_vertical_field!(B, p)
        B[1] = 1/2 * p[1] * dBdz
        B[2] = 1/2 * p[2] * dBdz
        B[3] = -c1*2 / p[3]^3 
        nothing
    end

    function convergent_vertical_field(p)
        B = [1/2 * p[1] * dBdz, 
             1/2 * p[2] * dBdz,
             -c1*2 / p[3]^3 ]
        return B
    end
    return convergent_vertical_field!, convergent_vertical_field
end


convergent_vertical_field_demo!, convergent_vertical_field_demo = make_convergent_vertical_field_demo(c)
"""
zz = [100, 200, 300, 400, 500]*1e3 .+ 500e3
yy = [-200, -100, 0, 100, 200.0] .*2
using GLMakie
p = [Point3([0, y, z]) for y in yy for z in zz]
v = Vec3.(convergent_vertical_field_demo.(p))
f = arrows3d(p, v, lengthscale = 1e6, axis=(type=Axis3, azimuth = 0, elevation = 0),)
save("figures/conicB.png", f)
"""

function make_convergent_vertical_field(c)
    z0 = (c.re + 80e3) #z0 at 80km above earth radius
    c1 = 8.22e15 #Tm3      c1 = mu_0 m / 4 pi ≈ m * 1e-7
    dBdz = -6*c1 / z0^4
    function convergent_vertical_field!(B, p)
        B[1] = 1/2 * p[1] * dBdz
        B[2] = 1/2 * p[2] * dBdz
        B[3] = -c1*2 / p[3]^3 
        nothing
    end

    function convergent_vertical_field(p)
        B = [1/2 * p[1] * dBdz, 
             1/2 * p[2] * dBdz,
             -c1*2 / p[3]^3 ]
        return B
    end
    return convergent_vertical_field!, convergent_vertical_field
end

convergent_vertical_field!, convergent_vertical_field = make_convergent_vertical_field(c)



function dipole_field_earth!(B, p)
    vm = @SVector [0, 0, -8.22e15] # Tm3
    dipole_field!(B, p, vm)
    nothing
end

function dipole_field!(B, p, vm)
    #so far the fastest!
    r = norm(p)
    B .= 1/r^5 * (3* p * dot(p, vm) - vm * r^2) 
    nothing
end

"""
B = @MArray zeros(Float64, 3)
p = @MArray rand(3)
@btime dipole_field_earth!(B, p)
"""



"""

function make_dipole_field_earth()
    vm = [0, 0, -8.22e15] # Tm3
    function dipole_field_earth!(B::Vector{Float64}, p::Vector{Float64})
        dipole_field!(B, p, vm)
        nothing
    end
    return dipole_field_earth!
end

"""

function dipole_field_earth(p::Vector{Float64})
    vm = [0, 0, -8.22e15] # Tm3
    return dipole_field(p, vm)
end


function dipole_field(p::Vector{Float64}, vm::Vector{Float64})
    #so far the fastest!
    r = norm(p)
    B = 1/r^5 * (3* p * dot(p, vm) - vm * r^2) 
    return B
end


"""

function dipole_field_(p, vm)
    # p - point to calculate the magnetic field
    # vm - moment vector for the magnetic field orientation
    #       8.22e15 Tm3 for earth
    #       vm = [0, 0, vm_z] with vm_z = mu_0 m / 4 pi ≈ m * 1e-7
    #           with m - magnetic dipole moment, m = 8.22e22 Am-2 for earth

    #returns B - magnetic field vector
    #               at point p
    #               of a magnetic moment vm
    #               assuming vm is at origin coordinate system

    r = norm(p)
    dipole = [
        3*p[1]^2 - r^2  3*p[1]*p[2]     3*p[1]*p[3];
        3*p[1]*p[2]     3*p[2]^2 - r^2  3*p[2]*p[3];
        3*p[1]*p[3]     3*p[2]*p[3]     3*p[3]^2 - r^2
        ]
    B = 1/r^5 * dipole * vm
    return B
end

#B = Vector{Float64}(undef, 3)

function dipole_field!(B::Vector{Float64}, p::Vector{Float64}, vm::Vector{Float64})
    #so far the fastest!
    r = norm(p)
    B .= 1/r^5 * (3* p * dot(p, vm) - vm * r^2) 
    nothing
end

using BenchmarkTools
# With normal vectors
B = zeros(3)
vm = [0, 0, -8.22e15] # Tm3
p = rand(3)
r = norm(p)
@btime dipole_field(p, vm)      # 219ns
@btime dipole_field(r, p, vm)   # 236ns
@btime dipole_field!(B, r, p, vm)   # 235ns
@btime dipole_field!(B, p, vm)      # 218ns

using StaticArrays
# With static arrays/vectors
B = @MArray zeros(3)
vm = @MArray [0, 0, -8.22e15] # Tm3
p = @MArray rand(3)
r = norm(p)
@btime dipole_field(p, vm)      # 40ns
@btime dipole_field(r, p, vm)   # 41ns
@btime dipole_field!(B, r, p, vm)   #36ns
@btime dipole_field!(B, p, vm)      # 34ns


function dipole_field4(p, vm)

    r = norm(p)
    B = 1/r^5 * (3* p * dot(p, vm)) - vm / r^3

    return B
end

function dipole_field2(p, vm)
    #half as fast :(
    r = norm(p)
    dipole = Hermitian([
        p[1]^2 - r^2  p[1]*p[2]     p[1]*p[3];
        0             p[2]^2 - r^2  p[2]*p[3];
        0             0             p[3]^2 - r^2
        ])
    B = 1/r^5 * 3 * dipole * vm
    return B
end
"""
