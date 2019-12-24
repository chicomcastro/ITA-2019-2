function [q3] = quatProd(q1,q2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

lambda1 = q1(1);
lambda2 = q2(1);

r1 = q1(2:4);
r2 = q2(2:4);

r3 = cross(r1(:),r2(:)) + lambda1*r2(:) + lambda2*r1(:);

q3 = [    lambda1*lambda2 - dot(r1,r2);    r3    ];

end

