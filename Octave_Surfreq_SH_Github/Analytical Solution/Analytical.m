%This code returns Analytical Solution of SH wave

function [Ush_anly,R] = Analytical(int,vs,rho)

vel=vs(1);
density=rho(1);
[X,Z] = meshgrid(int.x,int.z);
R = sqrt((X-int.dx*round(int.nx/2)).^2 +(Z-int.dx*round(int.nz/2)).^2);
ricker=(2.0./sqrt(pi)).*(int.mf^2/int.fc.^3).*exp(-(int.mf/ int.fc).^2);
%spike=1;
w=2*pi*int.mf;
const =(w*density)./4; %-1i./(4.*(vel.^2).*density);
H01 = besselh(0,1,w*R./vel);
%H01(1)=eps;
Green = const*H01;
Ush_anly=Green*ricker;

end


