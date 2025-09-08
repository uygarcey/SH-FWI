%This function calculates gardient and data residuals for each frequency.
%for this, parallel libary of OCTAVE is used.
function [Lb1]=func(int,f,vs, rho,S,vsint,rhoint)
% initialization of data residuals  (e) and gardient(g)
 g=zeros(int.nx*int.nz,1);
 e=0;
parfor i=1:length(f)
[Lb]=Grad(int,f(i),vs, rho,S,vsint,rhoint);

g=g+Lb.grad;
e=e+Lb.norm;
end
 Lb1.grad=g;
 Lb1.norm=e;





