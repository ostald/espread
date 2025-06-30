function rot_around_v(e_rot, phi)
    # ROT_AROUND_V - matrix for rotation PHI radians around E_ROT
    #   
    # Calling:
    # R = rot_around_v(e_rot,phi)
    # 
    # Multiply from left: v_out = R * v_in
    #
    #   Copyright 2002 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
    #   This is free software, licensed under GNU GPL version 2 or later
    
    r = [cos(phi) -sin(phi) 0;              #rotation matrix around z
         sin(phi)  cos(phi) 0;
         0         0        1];
    e_Rot = e_rot ./ sum(e_rot.^2)^(1/2);   #unit vector in e_rot direction
    az = atan(e_Rot[1], e_Rot[2]);
    ze = acos(e_Rot[3]);

    T = [1      0         0;
         0      cos(ze)  -sin(ze);
         0      sin(ze)   cos(ze)]*         #rotation around x
        [cos(az)  -sin(az)  0;
         sin(az)   cos(az)  0;
         0         0        1];             #rotation around z
    
    R = T'*r*T;
    return R
    
end

"""
phi = 0:0.1:2*pi
e_rot = [0, 0, 1]
v0 = [1, 0, 1]
ring = reduce(hcat, [rot_around_v(e_rot, p)*v0 for p in phi])

using GLMakie
cm = GLMakie

fig = cm.scatter(ring[1, :], ring[2, :], ring[3, :])
cm.arrows!(v0...)
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(v0)])
cm.arrows!([Point3f([0, 0, 0])], [Vec3f(e_rot)])
"""