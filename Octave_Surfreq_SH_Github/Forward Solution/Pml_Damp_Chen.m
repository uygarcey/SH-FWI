% PML.Reference paper
% Chen, Z., Cheng, D.,Feng,W.,and Wu,T. 2013.An optimal 9-point finite difference
% scheme for the Helmholtz equation with PML. International Journal Of Numerical
% Analysis and Modeling,10,389-410.
%In order to to create Pml_Damp_Chen function, The PML function(written in C language) in Germania repository was adapted to OCTAVE language.(No for loops are necessary in this version)
function[pml]= Pml_Damp_Chen( int, vs, rho,f )

freq = 2.0*pi*f;
lPML = int.pml_thc * int.dx;

ix=1:int.pml_thc-1;
jx=int.nx - int.pml_thc:int.nx;

damp_x=[freq .* int.a0_pml .* (((ix(end)+1)-ix).*int.dx./lPML).^2,...
    zeros(1,length(ix(end)+1:(jx(1))-1)),...
    freq .* int.a0_pml .* ((int.dx.*((jx(1)+1)-jx))./lPML).^2];

jz=int.nz - int.pml_thc:int.nz;
if (int.surface ~=1)
 damp_z=[freq .* int.a0_pml .* (((ix(end)+1)-ix).*int.dx./lPML).^2,...
    zeros(1,length(int.pml_thc:(jz(1))-1)),...
     freq .* int.a0_pml .* ((int.dx.*((jz(1)+1)-jz))./lPML).^2];
else
    damp_z=[zeros(1,length(ix)), zeros(1,length(int.pml_thc:(jz(1))-1)),...
     freq .* int.a0_pml .* ((int.dx.*((jz(1)+1)-jz))./lPML).^2];
end

sx=repmat(1-(damp_x(1:int.nx)*1i/freq),int.nz,1);
sz=repmat(1-(damp_z(1:int.nz)'*1i/freq),1,int.nx);
%
% Add q factor
%wref=int.fc;
 %q = 30;
%  vs = 1./((1./vs)+(1./(pi*vs.*q))*log(wref./freq) + 1i*0.5./(vs.*q)) ;

% Frequency-Domain Q-Compensated Reverse Time Migration Using a Stabilization Scheme
%  absor=(1./vs).*(1-(1./(pi*q))*log(wref./freq)).*(1-1i*0.5/q);
%  vs=1./absor;
% SH wave propagation in viscoelastic media
%P.F. Daley and E.S. Krebes
% dispersion=(1-(1./(pi*q))*log(freq./wref));
% vsd=vs.*(1./dispersion);
% vsa=(1./vsd).*(1-1i*0.5/q);
% vs=1./vsa;
%


pml.A = (1./rho).*sz./sx;
pml.B = (1./rho).*sx./sz;
pml.C = sx.*sz;

pml.Ax =pml.A;
pml.By = pml.B;

j=2:int.nz-1;
i=2:int.nx-1;

pml.Ax(j,i) =(pml.A(j,i) +pml.A(j,i+1)).*0.5;
pml.By(j,i) = (pml.B(j,i) + pml.B(j+1,i)).*0.5;


 mu = rho .* vs.*vs;


pml.b=0.7926;
pml.d=0.3768;
pml.e=-0.0064;
pml.freq2=freq^2;
pml.k2=pml.freq2*(1./mu);


if strcmp(int.parameter,'Velocity')
%  pml.dk2=2*pml.freq2.*(1./(rho.*vs.^3));
pml.dk2=(pml.freq2.*(2.*rho.*vs))./(mu.^2);
elseif strcmp(int.parameter,'Squared_Slowness')

 pml.dk2=pml.freq2.*(1./(rho));
end
pml.idh2=1.0/int.dx^2;

end





