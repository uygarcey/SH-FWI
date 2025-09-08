%This function returns Inr structure which contains corner coefficients
 %of the coefficient matrix and,  row and column numbers of them
 
function[Inr]= Interior(int,pml)

% Non zero node numbers
%Interior Nodes=(nx-2)*(nz-2)*int.stencil
%Upper+Lower Nodes=(nz-2)*(int.stencil-int.sd)*2
%Right + Left Nodes=(nx-2)*(int.stencil-int.sd)*2
%Corner Nodes =(int.stencil-(2*int.sd-1))*4
%
% Corner and Side Names of Calculation Domain
%[NW,N,NE,W,C,E,SW,S,SE]

% Corner and Side Names of 9-point Stencil and Indexing them
    %UL=(zz-2)*nx+ (xx-1); / (zz-1,xx-1) Upper Left
    %  U=(zz-2)*nx + xx;    / ((zz-1,xx)Upper
    %  UR=(zz-2)*nx + (xx+1); / (zz-1,xx+1) Upper Right
    %  L=(zz-1)*nx + xx-1;    /(zz,xx-1)Left
    %  C=(zz-1)*nx + xx;    /(zz,xx)Center
    %  R=(zz-1)*nx + (xx+1);    /(zz,xx) Right
    %  BL=(zz)*nx + (xx-1);     /  (zz+1,xx-1) Lower Left
    % B=(zz)*nx+xx; / (zz+1,xx) Lower
    % BR=(zz)*nx+(xx+1); / (zz+1,xx+1) Lower Right

nx=int.nx;
nz=int.nz;
index_IN=(nx-2)*(nz-2)*int.stencil;


b=pml.b;
d=pml.d;
e=pml.e;
k2=pml.k2;
dk2=pml.dk2;
idh2=pml.idh2;



j=2:nz-1;
i=2:nx-1;
  [Inx,Inz]=meshgrid(i,j);

 xx=reshape(Inx,length(i)*length(j),1);
 zz=reshape(Inz,length(i)*length(j),1);

    J_UL=(zz-2)*nx+ (xx-1);
      J_U=(zz-2)*nx + xx;
         J_UR=(zz-2)*nx + (xx+1);
         J_L=(zz-1)*nx + xx-1;
         J_C=(zz-1)*nx + xx;
        J_R=(zz-1)*nx + (xx+1);
        J_BL=(zz)*nx + (xx-1);
       J_B=(zz)*nx+xx;
       J_BR=(zz)*nx+(xx+1);

V_UL= reshape(((1-b).*idh2./2.0) .* (pml.Ax(j-1,i-1) + pml.By(j-1,i-1))  + (e./4.0) .* pml.C(j-1,i-1) .* k2(j-1,i-1),length(j)*length(i),1);
  V_U=reshape(-((1-b).*idh2) .* pml.A(j-1,i) + (b.*idh2) .* pml.By(j-1,i) + (d./4.0) .* pml.C(j-1,i) .* k2(j-1,i),length(j)*length(i),1);
  V_UR=reshape(((1-b).*idh2./2.0) .* (pml.Ax(j-1,i) + pml.By(j-1,i+1))+ (e./4.0) .* pml.C(j-1,i+1) .* k2(j-1,i+1),length(j)*length(i),1);

 V_L=reshape((b.*idh2) .* pml.Ax(j,i-1) - ((1-b).*idh2) .* pml.B(j,i-1)+ (d./4.0) .* pml.C(j,i-1) .* k2(j,i-1),length(j)*length(i),1);

 V_C=reshape((1-d-e) .* pml.C(j,i) .* k2(j,i)  - (2.0.*b.*idh2) .* (pml.A(j,i)+pml.B(j,i)),length(j)*length(i),1);

 V_R=reshape((b.*idh2) .* pml.Ax(j,i) - ((1-b).*idh2) .* pml.B(j,i+1) + (d./4.0) .* pml.C(j,i+1) .* k2(j,i+1),length(j)*length(i),1);

V_BL=reshape(((1-b).*idh2./2.0) .* (pml.Ax(j+1,i-1) + pml.By(j,i-1)) + (e./4.0) .* pml.C(j+1,i-1) .* k2(j+1,i-1),length(j)*length(i),1);

 V_B=reshape(-((1-b).*idh2) .* pml.A(j+1,i) + (b.*idh2) .* pml.By(j,i)  + (d./4.0) .* pml.C(j+1,i) .* k2(j+1,i),length(j)*length(i),1);

V_BR=reshape(((1-b).*idh2./2.0) .* (pml.Ax(j+1,i) + pml.By(j,i+1))  + (e./4.0) .* pml.C(j+1,i+1) .* k2(j+1,i+1),length(j)*length(i),1);
    % Derivative Matrix
  DV_UL= reshape((e./4.0) .* pml.C(j-1,i-1) .* dk2(j-1,i-1),length(j)*length(i),1);
  DV_U=reshape((d./4.0) .* pml.C(j-1,i) .* dk2(j-1,i),length(j)*length(i),1);
  DV_UR=reshape( (e./4.0) .* pml.C(j-1,i+1) .* dk2(j-1,i+1),length(j)*length(i),1);
 DV_L=reshape( (d./4.0) .* pml.C(j,i-1) .* dk2(j,i-1),length(j)*length(i),1);
 DV_C=reshape((1-d-e) .* pml.C(j,i) .* dk2(j,i) ,length(j)*length(i),1);
 DV_R=reshape( (d./4.0) .* pml.C(j,i+1) .* dk2(j,i+1),length(j)*length(i),1);
DV_BL=reshape( + (e./4.0) .* pml.C(j+1,i-1) .* dk2(j+1,i-1),length(j)*length(i),1);
 DV_B=reshape(- (d./4.0) .* pml.C(j+1,i) .* dk2(j+1,i),length(j)*length(i),1);
DV_BR=reshape( (e./4.0) .* pml.C(j+1,i+1) .* dk2(j+1,i+1),length(j)*length(i),1);

Inr.I=reshape(repmat(J_C,1,int.stencil)',index_IN,1);
Inr.J=reshape([J_UL,J_U,J_UR,J_L,J_C,J_R,J_BL,J_B,J_BR]',index_IN,1);
Inr.V=double(reshape([V_UL,V_U,V_UR,V_L,V_C,V_R,V_BL,V_B,V_BR].',index_IN,1));
Inr.DV=double(reshape([DV_UL,DV_U,DV_UR,DV_L,DV_C,DV_R,DV_BL,DV_B,DV_BR].',index_IN,1));
Inr.I=reshape(repmat(J_C,1,int.stencil)',index_IN,1);
Inr.J=reshape([J_UL,J_U,J_UR,J_L,J_C,J_R,J_BL,J_B,J_BR]',index_IN,1);
Inr.V=reshape([V_UL,V_U,V_UR,V_L,V_C,V_R,V_BL,V_B,V_BR].',index_IN,1);
end




