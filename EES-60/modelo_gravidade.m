%% Modelo gravidade

syms R_0 e g_0 h lambda
R_E = R_0*(1 + e*(sin(lambda))^2);              % raio leste-oeste
R_N = R_0*(1 - e*(2 - 3*(sin(lambda))^2));      % raio norte-sul
R_e = R_0*(1 - e*(sin(lambda))^2);              % raio terrestre

g = g_0*(1+.0053*(sin(lambda))^2)*(1-2*h/R_e);

%%
% Parâmetros do elipsóide de referência (WGS-84)
e = 1/298.25; % achatamento
R_0 = 6.378138e6; % m