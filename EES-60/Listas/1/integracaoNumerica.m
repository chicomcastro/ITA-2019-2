y = [
        q0_base3;
        [0,0,0]';
        lambda;
        Lambda;
        h;
        h
    ];

tf = 19.5*60;               % tempo final
i0 = 1;                     % index do início do movimento
i = i0;                     % index inicial
freqAquisicao = 100;        % Hz
p = 1/freqAquisicao;        % passo de integração
result = zeros(tf/p,length(y)+2);

h_m = h;            % constante e igual à altitude inicial
estabVert.B = 10;
estabVert.C = 10;
estabVert.T_h = 60;

modTerra.R_0 = R_0; % [m]
modTerra.g_0 = g_0; % [m/s^2]
modTerra.e = e;     % achatamento

for t=0:p:tf
    omega_B_i = w_b(i,:)';  % retirado dos dados a cada passo
    Asp_B_i = Asp_b(i,:)';  % idem
    
    k1 = f(t, y, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k2 = f(t + p/2, y + p/2*k1, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k3 = f(t + p/2, y + p/2*k2, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);
    k4 = f(t + p, y + p*k3, omega_B_i, Asp_B_i, modTerra, estabVert, h_m);

    y = y + p/6*(k1 + 2*k2 + 2*k3 + k4);
    result(i-i0+1,1) = t;       % tempo
    result(i-i0+1,2:12) = y;    % vetor de estados
    i = i+1;
end