function [y_new] = fs(t, y, omega_B_i, Asp_B_i, modTerra, estabVert, h_m, T)

    %% Preparação
    N = 1;
    E = 2;
    D = 3;

    R_0 = modTerra.R_0;
    e = modTerra.e;
    g_0 = modTerra.g_0;

    B = estabVert.B;
    C = estabVert.C;
    T_h = estabVert.T_h;

    Omega = 7.29e-5;

    %% Mapeamento de variáveis
    q_NEDold_bold = y(1:4);
    V_NED = y(5:7);
    lat = y(8);
    lon = y(9);
    alt = y(10);
    h_aux = y(11);
    
    %% Passo 0: quaternion de rotação do NEDold para NEDnew a cada 4 amostras de freqüência rápida dos sensores inerciais, computado com frequência lenta
    %% Modelo da Terra
    R_E = R_0*(1 + e*(sin(lat))^2);            % raio leste-oeste
    R_N = R_0*(1 - e*(2 - 3*(sin(lat))^2));    % raio norte-sul

    rho_NED = [...
        V_NED(E)/(R_E+alt), ...
        V_NED(N)/(R_N+alt), ...
        V_NED(E)/(R_E+alt)*tan(lat) ...
        ]';
    Omega_NED = [Omega*cos(lat),0,-Omega*sin(lat)]';

    omega_NEDi_NED = rho_NED + Omega_NED;  % taxa de transporte usa estimativas de velocidade terrestre e posição mais
     % recentes disponíveis 

    omega_versor = omega_NEDi_NED/norm(omega_NEDi_NED); % eixo de rotação instantânea unitário

    q_NEDold_NEDnew = [...
         cos(norm(omega_NEDi_NED*4*T)/2);...
         omega_versor(:)*sin(norm(omega_NEDi_NED*4*T)/2)
         ];

    %% Passo 1: incremento da velocidade de empuxo ?Uf ,b,k computada com frequência rápida. Índice k representa o instante inicial de um conjunto de 4 amostras; k+1 o instante inicial do próximo conjunto de 4 amostras, sem superposição com o anterior.
    alpha(:,1) = omega_B_i(:,1);
    alpha(:,2) = omega_B_i(:,2);
    alpha(:,3) = omega_B_i(:,3);
    alpha(:,4) = omega_B_i(:,4);

    delta_beta(:,1) = Asp_B_i(:,1);
    delta_beta(:,2) = Asp_B_i(:,2);
    delta_beta(:,3) = Asp_B_i(:,3);
    delta_beta(:,4) = Asp_B_i(:,4);

    W_k_0 = zeros(3,1);   % W_k(0) = (0,0,0)^T
    W_k_old = W_k_0;

    for m=1:4 % sculling correction
        W_k_new = delta_beta(:,m) - cross(alpha(:,m), W_k_old) + W_k_old;
        W_k_new = delta_beta(:,m) - cross(alpha(:,m), W_k_new) + W_k_old;
        W_k_old = W_k_new;
    end

    delta_U_f_b_k = W_k_new;

    %% Passo 2: Incremento angular computado com quatro amostras incrementais e quaternion de rotação bold qbnew do corpo na atitude anterior (bold) para o corpo na nova atitude (bnew). (new é a estampa de tempo k+1 e old é a estampa de tempo k) 
    % coning correction
    P1 = crossToMatrix(alpha(:,1));
    P2 = crossToMatrix(alpha(:,2));
    P3 = crossToMatrix(alpha(:,3));
    P4 = crossToMatrix(alpha(:,4));

    delta_phi = alpha(:,1)+alpha(:,2)+alpha(:,3)+alpha(:,4)+...
        2/3*(P1*alpha(:,2)+P3*alpha(:,4))+...
        1/2*(P1+P2)*(alpha(:,3)+alpha(:,4))+...
        1/30*(P1-P2)*(alpha(:,3)-alpha(:,4));

    delta_phi_versor = delta_phi/norm(delta_phi);

    q_bold_bnew = [...
         cos(norm(delta_phi)/2);...
         delta_phi_versor(:)*sin(norm(delta_phi)/2)
         ];

    % quaternion associado a ?Uf ,b,k
    delta_U_f_b_k_q = quat(delta_U_f_b_k);

    %% Passo 3: Computa em freqüência lenta o quaternion de rotação bnew qNEDnew mediante atualização de bold qNEDold devido à rotação de SNED para posterior transformação da força específica de Sb para SNED. (new é a estampa de tempo k+1 e old é a estampa de tempo k)
    q_NEDnew_bnew = quatInv(...
        quatProd(...
            quatProd(...
                quatInv(q_bold_bnew), ...
                quatInv(q_NEDold_bold)), ...
            q_NEDold_NEDnew)...
        );
    delta_U_f_NED_k_q = quatProd(...
        quatProd(...
            q_NEDnew_bnew, ...
            delta_U_f_b_k_q), ...
        quatInv(q_NEDnew_bnew)...
        );

    %% Passo 4: Atualiza em baixa freqüência (a cada período 4T) velocidade terrestre e posição. Usa posição ? k( ? ),1 k(h ? )1 e velocidade terrestre VNED,k?1 disponíveis para computar incremento de velocidade terrestre devido à gravidade.
    R_e = R_0*(1 - e*(sin(lat))^2);
    g = g_0*(1+.0053*(sin(lat))^2)*(1-2*alt/R_e);
    g_NED = [0,0,g]';

    delta_Vg_NED = (-cross((rho_NED + 2*Omega_NED), V_NED) + g_NED + [0,0,B*(alt-h_aux)]')*4*T;
    delta_V_NED = delta_Vg_NED + quatToVector(delta_U_f_NED_k_q);
    V_NED = V_NED + delta_V_NED;  %V_NED_new

    V_N = V_NED(N);
    V_E = V_NED(E);
    V_D = V_NED(D);

    % Atualizar lambda, Lambda e h para o instante k
    lat = lat + V_N/(R_N+alt);
    lon = lon + V_E/((R_E+alt)*cos(lat));
    alt = alt + -V_D - C*(alt-h_aux);
    h_aux = h_aux + (h_m-h_aux)/T_h;
    
    y_new = [...
        q_NEDnew_bnew; ...
        V_NED; ...
        lat; ...
        lon; ...
        alt; ...
        h_aux ...
    ];
end

