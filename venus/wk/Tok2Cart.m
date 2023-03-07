%Takes TOKAMAK cylindrical coords (R,Z,phi) and maps them into cartesian
%coords (x,y,z)

function [x,y,z]=Tok2Cart(R,Z,phi);

x=R.*cos(phi);
y=R.*sin(phi);
z=Z;