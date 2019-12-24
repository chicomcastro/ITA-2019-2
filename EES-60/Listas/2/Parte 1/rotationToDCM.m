function [DCM] = rotationToDCM(rotation)

rotation = rotation(:);
psi = rotation(1);
theta = rotation(2);
phi = rotation(3);

Dx = [
    1 0 0
    0 cos(psi) sin(psi)
    0 -sin(psi) cos(psi)
    ];

Dy = [
    cos(theta) 0 -sin(theta)
    0 1 0
    sin(theta) 0 cos(theta)
    ];

Dz = [
    cos(phi) sin(phi) 0
    -sin(phi) cos(phi) 0
    0 0 1
    ];

DCM = Dz*Dy*Dx;

end

