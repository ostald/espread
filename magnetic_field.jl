using LinearAlgebra


function dipole_field_earth(p::Vector{Float64})
    vm = [0, 0, -8.22e15] # Tm3
    return dipole_field(p, vm)
end


function dipole_field_earth!(B::Vector{Float64}, p::Vector{Float64})
    vm = [0, 0, -8.22e15] # Tm3
    dipole_field!(B, p, vm)
    nothing
end


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


B = Vector{Float64}(undef, 3)

function dipole_field(p::Vector{Float64}, vm::Vector{Float64})
    #so far the fastest!
    r = norm(p)
    B = 1/r^5 * (3* p * dot(p, vm) - vm * r^2) 
    return B
end




function dipole_field!(B::Vector{Float64}, p::Vector{Float64}, vm::Vector{Float64})
    #so far the fastest!
    r = norm(p)
    B = 1/r^5 * (3* p * dot(p, vm) - vm * r^2) 
    nothing
end


"""
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
