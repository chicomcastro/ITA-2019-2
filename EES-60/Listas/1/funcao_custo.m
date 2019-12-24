B = x(1);
C = x(2);
T_h = x(3);

%% Integração numérica

y = [
    q0_base3;
    [0,0,0]';
    lambda;
    Lambda;
    h;
    h
    ];

%% Condições de parada
tf = 800;                   % tempo final
i0 = 36000;                 % index do início do movimento
i = i0;                     % index inicial
freqAquisicao = 100;        % Hz
p = 1/freqAquisicao;        % passo de integração
result = zeros(tf/p,length(y)+1);

%% Runge-kutta
for t=0:p:tf
    omega_B_i = w_b(i,:)';
    Asp_B_i = Asp_b(i,:)';
    
    k1 = f(t, y, omega_B_i, Asp_B_i, R_0, g_0, e, B, C, T_h);
    k2 = f(t + p/2, y + p/2*k1, omega_B_i, Asp_B_i, R_0, g_0, e, B, C, T_h);
    k3 = f(t + p/2, y + p/2*k2, omega_B_i, Asp_B_i, R_0, g_0, e, B, C, T_h);
    k4 = f(t + p, y + p*k3, omega_B_i, Asp_B_i, R_0, g_0, e, B, C, T_h);

    y = y + p/6*(k1 + 2*k2 + 2*k3 + k4);
    result(i-i0+1,1) = t;       % tempo
    result(i-i0+1,2:11) = y;    % vetor de estados
    i = i+1;
end