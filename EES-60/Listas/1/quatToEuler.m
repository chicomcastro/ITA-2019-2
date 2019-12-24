function [euler] = quatToEuler(q)

q = q(:);
qw = q(1);
qx = q(2);
qy = q(3);
qz = q(4);

% roll (x-axis rotation)
sinr_cosp = 2.0 * (qw * qx + qy * qz);
cosr_cosp = 1.0 - 2.0 * (qx * qx + qy * qy);
euler(3) = atan2(sinr_cosp, cosr_cosp);

% pitch (y-axis rotation)
sinp = 2.0 * (qw * qy - qz * qx);
if (abs(sinp) >= 1)
    euler(2) = (pi / 2) * sign(sinp); % use 90 degrees if out of range
else
    euler(2) = asin(sinp);
end

% yaw (z-axis rotation)
siny_cosp = +2.0 * (qw * qz + qx * qy);
cosy_cosp = +1.0 - 2.0 * (qy * qy + qz * qz);  
euler(1) = atan2(siny_cosp, cosy_cosp);

euler = euler(:);

end

