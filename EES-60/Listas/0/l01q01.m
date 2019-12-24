r_I = [
    1;
    0;
    0
    ];
r_I_q = quat(r_I);

q_I_V = rotationToQuat(pi/2,[0,0,1]);
q1 = quatConj(q_I_V);
q2 = q_I_V;

r_V_q = quatProd(quatProd(q1,r_I_q),q2)

r_V = r_V_q(2:4)