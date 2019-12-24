function [x,y,z] = geodToCart(lat,lon,h,modTerra)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

e = modTerra.e;
f = 0.0818;
R_0 = modTerra.R_0;
x = (R_0./sqrt(1-(e*sin(lat)).^2)+h).*cos(lat).*cos(lon);
y = (R_0./sqrt(1-(e*sin(lat)).^2)+h).*cos(lat).*sin(lon);
z = (R_0*(1-f^2)./sqrt(1-(e*sin(lat)).^2)+h).*cos(lat);

wgs84 = wgs84Ellipsoid('meter');
[x,y,z] = geodetic2ecef(wgs84,lat,lon,h, 'radians');

wgs84 = wgs84Ellipsoid('meter');
[x,y,z] = geodetic2enu(lat,lon,h,lat(1263),lon(1263),h(1263),wgs84,'radians');

end

