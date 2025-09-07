% This function returns Right Hand Side
function [ Rhs] = Rhs_Mul_Shot(int,f)
int.mf=f; %current frequncy
Rhs = sparse(double(int.nz*int.nx),double(int.nsrc));
src_index = int.srcx + (int.srcz-1)*int.nx;

%f0=int.fc;
f=int.mf;
  % t0 = -2/(pi*f0);
    %% Wavelets
 %ricker2= 1i*(2*pi*f).*sqrt(pi/f0).*exp(-.5*(f/f0).^2 + 1i*(2*pi*f)*t0);
 ricker1= sqrt(4/(pi*int.fc^2))*((f/int.fc).^2).*exp(-(f/int.fc).^2).*exp(-1i*2*pi*f/int.fc);
 %amp = (2.0./sqrt(pi)).*(f.^2/int.fc.^3).*exp(-(f/int.fc).^2);
Rhs(sub2ind([int.nx*int.nz,int.nsrc],src_index(1:int.nsrc),1:int.nsrc)) = ricker1./int.dx^2;

end




