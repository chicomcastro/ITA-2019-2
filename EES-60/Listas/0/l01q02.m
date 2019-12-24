%% Questão 02

% S_I -(pi/2;3)-(pi/2;2-> S_V

q_IA = rotationToQuat(pi/2,[0,0,1]);
q_AV = rotationToQuat(pi/2,[0,1,0]);

q_IV = quatProd(q_IA,q_AV)  % qual a ordem??

phi_IV = quatToRotation(q_IV)

r_I_q = [0,10,-5,0];
r_V_q = quatProd(quatProd(quatConj(q_IV),r_I_q),q_IV);
r_v = r_V_q(2:4)