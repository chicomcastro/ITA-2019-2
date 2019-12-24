function [y_ponto] = f(t, y, omega_B_i, Asp_B_i, modTerra, estabVert, h_m)

    %% Mapeamento de variáveis
    q = y(1:4);
    V_N = y(5);
    V_E = y(6);
    V_D = y(7);
    lambda = y(8);
    Lambda = y(9);
    h = y(10);
    h_aux = y(11);
    
    
    %% Parâmetros auxiliares e valores de apoio
    x = 1;
    y = 2;
    z = 3;
    N = 1;
    E = 2;
    D = 3;
    Omega = 7.292115e-5;

    
    %% Modelo da Terra
    R_0 = modTerra.R_0;
    e = modTerra.e;
    R_E = R_0*(1 + e*(sin(lambda))^2);            % raio leste-oeste
    R_N = R_0*(1 - e*(2 - 3*(sin(lambda))^2));    % raio norte-sul
    R_e = R_0*(1 - e*(sin(lambda))^2);

    
    %% Modelo da gravidade
    %g_0, R_0, e foram definidos anteriormente
    g_0 = modTerra.g_0;
    g = g_0*(1+.0053*(sin(lambda))^2)*(1-2*h/R_e);

    
    %% Navegação inercial
    lambda_c_ponto = V_N/(R_N+h);
    Lambda_c_ponto = V_E/((R_E+h)*cos(lambda));
    omega_NED_i = [
        (Omega+Lambda_c_ponto)*cos(lambda);
        -lambda_c_ponto;
        -(Omega+Lambda_c_ponto)*sin(lambda)
    ];


    %% Determinação de atitude
    q_NED_B_comp = q;
    q_NED_B_comp = q_NED_B_comp/norm(q_NED_B_comp);
    Asp_B_medido = Asp_B_i;
    Asp_B_medido_q = [0;Asp_B_medido];
    
    Asp_NED_medido_q = quatProd(quatProd(q_NED_B_comp,Asp_B_medido_q),quatConj(q_NED_B_comp));
    Asp_NED_medido = Asp_NED_medido_q(2:4);
    
    Omega_B_i = [
        0               -omega_B_i(x)    -omega_B_i(y)   -omega_B_i(z);
        omega_B_i(x)    0                omega_B_i(z)    -omega_B_i(y);
        omega_B_i(y)    -omega_B_i(z)    0               omega_B_i(x);
        omega_B_i(z)    omega_B_i(y)     -omega_B_i(x)   0
        ];

    Omega_NED_i = [
        0                  omega_NED_i(x)     omega_NED_i(y)    omega_NED_i(z);
        -omega_NED_i(x)    0                  omega_NED_i(z)    -omega_NED_i(y);
        -omega_NED_i(y)    -omega_NED_i(z)    0                 omega_NED_i(x);
        -omega_NED_i(z)    omega_NED_i(y)     -omega_NED_i(x)   0
        ];
    
    %% Estabilização vertical
    B = estabVert.B;
    C = estabVert.C;
    T_h = estabVert.T_h;
    
    %% Output
    % y' = f(t,y)
    y_ponto = [
        1/2*(Omega_B_i+Omega_NED_i)*q_NED_B_comp;
        Asp_NED_medido(N) + (V_N/(R_N + h))*V_D - (2*Omega*sin(lambda) + V_E/(R_E + h)*tan(lambda))*V_E;
        Asp_NED_medido(E) + (2*Omega*sin(lambda) + V_E/(R_E + h)*tan(lambda))*V_N + (2*Omega*cos(lambda) + V_E/(R_E + h))*V_D;
        Asp_NED_medido(D) - (2*Omega*cos(lambda) + V_E/(R_E + h))*V_E - (V_N/(R_N + h))*V_N + g + B*(h-h_aux);
        V_N/(R_N+h);
        V_E/((R_E+h)*cos(lambda));
        -V_D - C*(h-h_aux);
        (h_m-h_aux)/T_h
    ];
end

