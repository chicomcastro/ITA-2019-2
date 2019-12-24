figure
plot3(X0,Y0,Z0,'o')
hold on

t = 0:0.01:1;

for i = 1:4
    
plot3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'x')
hold on
plot3(X0 + r_e(1,i)*t,Y0 + r_e(2,i)*t,Z0 + r_e(3,i)*t)
hold on

end

grid on
t = t*1e4;
x_ned = D_eNED(LAT_receptor,LON_receptor)*[1 0 0]'*t;
hold on;
plot3(X0 + x_ned(1,i)*t,Y0 + x_ned(2,i)*t,Z0 + x_ned(3,i)*t);

x_ned = D_eNED(LAT_receptor,LON_receptor)*[0 1 0]'*t;
hold on;
plot3(X0 + x_ned(1,i)*t,Y0 + x_ned(2,i)*t,Z0 + x_ned(3,i)*t);

x_ned = D_eNED(LAT_receptor,LON_receptor)*[0 0 1]'*t;
hold on;
plot3(X0 + x_ned(1,i)*t,Y0 + x_ned(2,i)*t,Z0 + x_ned(3,i)*t);
