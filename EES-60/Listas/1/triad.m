

%% DCM inicial utilizando o método TRIAD
syms a_x a_y a_z w_x w_y w_z
Asp_b = [a_x ; a_y ; a_z];
Omega_b = [w_x; w_y; w_z];

syms g Omega lambda
Asp_NED = [0;0;-g];
Omega_NED = [Omega*cos(lambda); 0 ; -Omega*sin(lambda)];

%%
% Base 1: Asp, ? e Asp×?. (base não ortonormal)
A = [Asp_NED, Omega_NED, cross(Asp_NED,Omega_NED)];
B = [Asp_b, Omega_b, cross(Asp_b,Omega_b)];

D_NED_B_1 = B*inv(A);

%%
% Base 2:  Asp, Asp×? e Asp×(Asp×?). (base ortogonal, mas não normalizada)
A = [Asp_NED, cross(Asp_NED,Omega_NED), cross(Asp_NED,cross(Asp_NED,Omega_NED))];
B = [Asp_b, cross(Asp_b,Omega_b), cross(Asp_b,cross(Asp_b,Omega_b))];

D_NED_B_2 = B*inv(A);

%%
% Base 3: Asp/|Asp|, Asp×?/|Asp×?| e Asp×(Asp×?)/|Asp×(Asp×?)|. (base ortonormalizada)
A = [Asp_NED/norm(Asp_NED), cross(Asp_NED,Omega_NED)/norm(cross(Asp_NED,Omega_NED)), cross(Asp_NED,cross(Asp_NED,Omega_NED))/norm(cross(Asp_NED,cross(Asp_NED,Omega_NED)))];
B = [Asp_b/norm(Asp_b), cross(Asp_b,Omega_b)/norm(cross(Asp_b,Omega_b)), cross(Asp_b,cross(Asp_b,Omega_b))/norm(cross(Asp_b,cross(Asp_b,Omega_b)))];

D_NED_B_3 = B*inv(A);
