%% Resultados para as DCMs iniciais
%%
% Base 1
D0_NED_B_base1 = double(subs(subs(D_NED_B_1)));
q0_base1 = DCMtoQuaternion(D0_NED_B_base1);
euler0_base1 = quatToEuler(q0_base1);

%%
% Base 2
D0_NED_B_base2 = double(subs(subs(D_NED_B_2)));
q0_base2 = DCMtoQuaternion(D0_NED_B_base2);
euler0_base2 = quatToEuler(q0_base2);

%%
% Base 3
D0_NED_B_base3 = double(subs(subs(D_NED_B_3)));
q0_base3 = DCMtoQuaternion(D0_NED_B_base3);
euler0_base3 = quatToEuler(q0_base3);
