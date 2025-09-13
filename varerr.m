%rtklibÖÐµÄnorm error model
function [var]=varerr(eli,bl,dt)

Clight=299792458.0;
sclkstab=5e-12;
c=0*bl/1e4;
d=Clight*sclkstab*dt;
fact=1;

a=fact*0.003;
b=fact*0.003;

var=2*(a*a+b*b/sin(eli)/sin(eli)+c*c)+d*d;
end