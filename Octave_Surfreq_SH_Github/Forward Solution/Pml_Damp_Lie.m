% PML. Reference paper
%Lei,W.,Liu,Y.,Li,G.,Zhu,S.,Chen,G., and Li,G. 2023.2D frequency-domain finite
%difference acoustic wave modeling using optimized perfectlly matched layers.
%Geophysics, 88, F1-F13.
function [ pml ] = Pml_Damp_Lie( int, vs, rho,f )


int.omega = 2.0*pi*f;

sigma_max=(3*pi*(int.m2_pml+1).*min(min(vs))./(10*int.dx));

lPML = int.pml_thc * int.dx;

ix=1:int.pml_thc-1;
jx=int.nx - int.pml_thc:int.nx;



sigma_x=[sigma_max.* (((ix(end)+1)-ix).*int.dx./lPML).^int.m2_pml,...
    zeros(1,length(ix(end)+1:(jx(1)-1))),...
    sigma_max.* (((int.dx.*((jx(1)+1)-jx))./lPML).^int.m2_pml)];

beta_x=[int.b0_pml*int.fc.* (1-((ix(end)-ix).*int.dx./lPML).^int.n_pml),...
    zeros(1,length(ix(end)+1:(jx(1))-1)),...
    int.b0_pml*int.fc.*((lPML-((int.dx.*((jx(1)+1)-jx))))./lPML).^int.n_pml];

jz=int.nz - int.pml_thc:int.nz;

if (int.surface ~=1)
 sigma_z=[sigma_max.* ((ix(end)-ix).*int.dx./lPML).^int.m2_pml,...
    zeros(1,length(int.pml_thc:(jz(1))-1)),...
     sigma_max.* ((int.dx.*((jz(1)+1)-jz))./lPML).^int.m2_pml];

 beta_z=[int.b0_pml*int.fc.* (1-((ix(end)-ix).*int.dx./lPML).^int.n_pml),...
    zeros(1,length(int.pml_thc:(jz(1))-1)),...
    int.b0_pml*int.fc.*((lPML-((int.dx.*((jz(1)+1)-jz))))./lPML).^int.n_pml];
else
    sigma_z=[zeros(1,length(ix)), zeros(1,length(int.pml_thc:(jz(1))-1)),...
     sigma_max.* ((lPML-((int.dx.*((jz(1)+1)-jz)))./lPML)).^int.m2_pml];

  beta_z=[zeros(1,length(ix)),...
    zeros(1,length(int.pml_thc:(jz(1))-1)),...
    int.b0_pml*int.fc.*((lPML-((int.dx.*((jz(1)+1)-jz)))./lPML)).^int.m2_pml];
end

sx=repmat(1+(sigma_x(1:int.nx)./(beta_x(1:int.nx)+(1i*int.omega))),int.nz,1);
sz=repmat(1+(sigma_z(1:int.nz)'./(beta_z(1:int.nz)'+(1i*int.omega))),1,int.nx);
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
pml.freq2=int.omega^2;
pml.k2=pml.freq2*(1./mu);
%*(1./(mu.^2));

if strcmp(int.parameter,'Velocity')
%  pml.dk2=2*pml.freq2.*(1./(rho.*vs.^3));
pml.dk2=(pml.freq2.*(2.*rho.*vs))./(mu.^2);
elseif strcmp(int.parameter,'Squared_Slowness')

 pml.dk2=pml.freq2.*(1./(rho));
end
pml.idh2=1.0/int.dx^2;


end


