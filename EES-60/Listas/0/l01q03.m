%% a) Sequência de ângulos:
%  E  ->  N  ->  E
% pi/4   pi/4   pi/4
N = [0,1,0];
E = [1,0,0];

q1 = rotationToQuat(pi/4, E);
q2 = rotationToQuat(pi/4, N);
q3 = rotationToQuat(pi/4, E);

q_NEDo_NEDf = quatProd(quatProd(q3,q2),q1);

%% b) DCM

DCM_NEDo_NEDf = quatToDCM(q_NEDo_NEDf)


%% c) Vetor de Euler

phi_NEDo_NEDf = quatToRotation(q_NEDo_NEDf)