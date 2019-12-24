%% Análise de resultados para os dados não incrementais

format shortg

%% Verificação de resultados
% *Ortonormalidade das DCMs*
determinantesDCM = [det(D0_NED_B_base1), det(D0_NED_B_base2), det(D0_NED_B_base3)]

%%
% *Normalidade dos quaternion*
normasQuat = [norm(q0_base1), norm(q0_base2), norm(q0_base3)]
%%
% Ou seja, os quaternions, com exceção do último, não são de rotação.


%% Análises de normalização e ortonormalização iterativa
% *Ortonormalização iterativa*
D0_base1_orto = ortonormalizacaoIterativo(D0_NED_B_base1);
D0_base2_orto = ortonormalizacaoIterativo(D0_NED_B_base2);
D0_base3_orto = ortonormalizacaoIterativo(D0_NED_B_base3);
%%
desalinhamentoOrtoIt = [...
    desalinhamento(D0_base1_orto, D0_NED_B_base1),...
    desalinhamento(D0_base2_orto, D0_NED_B_base2),...
    desalinhamento(D0_base3_orto, D0_NED_B_base3) ...
    ]

%%
% *Normalização iterativa*
%%
q_rot1 = normalizacaoIterativa(q0_base1);	% Gera um quaternion de rotação
nova_DCM1 = quatToDCM(q_rot1);              % Calcula a DCM associada

%%
q_rot2 = normalizacaoIterativa(q0_base2);	% Gera um quaternion de rotação
nova_DCM2 = quatToDCM(q_rot2);              % Calcula a DCM associada

%%
q_rot3 = normalizacaoIterativa(q0_base3);	% Gera um quaternion de rotação
nova_DCM3 = quatToDCM(q_rot3);              % Calcula a DCM associada

%%
desalinhamentoNormIt = [...
    desalinhamento(nova_DCM1, D0_NED_B_base1),... % Calcula o desalinhamento
    desalinhamento(nova_DCM2, D0_NED_B_base2),... % Calcula o desalinhamento
    desalinhamento(nova_DCM3, D0_NED_B_base3) ... % Calcula o desalinhamento
    ]

%%
desalinhamentoCruzado = [...
    desalinhamento(D0_NED_B_base1, D0_NED_B_base3),...
    desalinhamento(D0_NED_B_base2, D0_NED_B_base3),...
    desalinhamento(D0_NED_B_base3, D0_NED_B_base3);...
    desalinhamento(D0_base1_orto, D0_NED_B_base3),...
    desalinhamento(D0_base2_orto, D0_NED_B_base3),...
    desalinhamento(D0_base3_orto, D0_NED_B_base3);...
    desalinhamento(nova_DCM1, D0_NED_B_base3),...
    desalinhamento(nova_DCM2, D0_NED_B_base3),...
    desalinhamento(nova_DCM3, D0_NED_B_base3)...   
    ]
    
