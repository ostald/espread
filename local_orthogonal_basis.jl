using LinearAlgebra

function local_orthogonal_basis(vector)
    # returns an orthogonal basis with the supplied vector along the last base vector
    xu = [1, 0, 0]
    if vector/norm(vector) == xu
        xu = [0, 1, 0]
    end
    n1 = cross(xu, vector)
    n1u = n1/norm(n1)
    n2 = cross(n1u, vector)
    n2u = n2/norm(n2)
    n3u = vector/norm(vector) # the parellel base vector
    return n1u, n2u, n3u
end