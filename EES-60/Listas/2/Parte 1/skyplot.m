function skyplot(azimuth,elevation)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

axis = polaraxes;

polarplot(axis, deg2rad(azimuth), elevation);

axis.ThetaDir = 'clockwise';
axis.ThetaZeroLocation = 'top';

axis.RLim = [0 90];
axis.RDir = 'reverse';
grid on

end

