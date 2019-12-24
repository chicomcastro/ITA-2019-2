syms lambda%% = 23.20995*pi/180;   %rad

syms ax ay az wx wy wz

g = sqrt(mean(ax)^2+mean(ay)^2+mean(az)^2)
Omega = sqrt(mean(wx)^2+mean(wy)^2+mean(wz)^2)

Asp_b = [ax; ay; az];
Omega_b = [wx; wy; wz];
p1 = cross(Asp_b, Omega_b);

Asp_ned = [0; 0; -g];
Omega_ned = [Omega*cos(lambda); 0; -Omega*sin(lambda)]
p2 = cross(Asp_ned, Omega_ned);

B = [Asp_b Omega_b p1];
A = [Asp_ned Omega_ned p2];

Db_ned = B*inv(A)