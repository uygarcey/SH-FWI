%This function returns Crs structure which contains corner parts of the coefficient matrix with row and column numbers of them.
%Developed by Uygar Ceyhan
function[Crs]= Corners(int,pml)
int.sd=3; %Square root of 9
b=pml.b;
d=pml.d;
e=pml.e;
k2=pml.k2;
dk2=pml.dk2;
idh2=pml.idh2;

nx=int.nx;
nz=int.nz;


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

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corner Nodes
% I indicates rows, J indicates columns, V indicates values

index_CR= (int.stencil-(2*int.sd-1));
% NW (1,1)

J_NWC= 1;
J_NWR=2;
J_NWB=nx+1;
J_NWBR=nx+2;


V_NWC=(1-d-e) .* pml.C(1,1) .* k2(1,1)- (2.0.*b.*idh2) .* (pml.A(1,1)+pml.B(1,1));
V_NWR=(b.*idh2) .* pml.Ax(1,1) - ((1-b).*idh2) .* pml.B(1,2) + (d./4.0) .* pml.C(1,2) .* k2(1,2);
V_NWB=-((1-b).*idh2) .* pml.A(2,1) + (b.*idh2) .* pml.By(1,1) + (d./4.0) .* pml.C(2,1) .* k2(2,1);
V_NWBR=((1-b).*idh2./2.0) .* (pml.Ax(2,1) + pml.By(1,2)) + (e./4.0) .* pml.C(2,2) .* k2(2,2);

% Derivative matrix (It is used for gradient calculation)
DV_NWC=(1-d-e) .* pml.C(1,1) .* dk2(1,1);
DV_NWR= (d./4.0) .* pml.C(1,2) .* dk2(1,2);
DV_NWB= (d./4.0) .* pml.C(2,1) .*dk2(2,1);
DV_NWBR= (e./4.0) .* pml.C(2,2) .* dk2(2,2);
%
I_NW=repmat(J_NWC,index_CR,1);
J_NW=[J_NWC;J_NWR;J_NWB;J_NWBR];
V_NW=[V_NWC;V_NWR;V_NWB;V_NWBR];
DV_NW=[DV_NWC;DV_NWR;DV_NWB;DV_NWBR];

%  NE (1,nx)
  J_NEL=nx-1;
  J_NEC=nx;
  J_NEBL=2*nx-1;
  J_NEB=2*nx;

  V_NEL=(b.*idh2) .* pml.Ax(1,nx-1) - ((1-b).*idh2) .* pml.B(1,nx-1)  + (d./4.0) .* pml.C(1,nx-1) .* k2(1,nx-1);
  V_NEC=(1-d-e) .* pml.C(1,nx) .* k2(1,nx)- (2.0.*b.*idh2) .* (pml.A(1,nx)+pml.B(1,nx));
  V_NEBL=((1-b).*idh2./2.0) .* (pml.Ax(2,nx-1) + pml.By(1,nx-1))+ (e./4.0) .* pml.C(2,nx-1) .* k2(2,nx-1);
  V_NEB=-((1-b).*idh2) .* pml.A(2,nx) + (b.*idh2) .* pml.By(1,nx) + (d./4.0) .* pml.C(2,nx) .* k2(2,nx);
  % Derivative matrix
   DV_NEL=(d./4.0) .* pml.C(1,nx-1) .* dk2(1,nx-1);
  DV_NEC=(1-d-e) .* pml.C(1,nx) .* dk2(1,nx);
  DV_NEBL= (e./4.0) .* pml.C(2,nx-1) .* dk2(2,nx-1);
  DV_NEB= (d./4.0) .* pml.C(2,nx) .* dk2(2,nx);
    %
I_NE=repmat(J_NEC,index_CR,1) ;
J_NE=[J_NEL;J_NEC; J_NEBL;J_NEB];
V_NE=[ V_NEL;V_NEC; V_NEBL;V_NEB];

DV_NE=[ DV_NEL;DV_NEC; DV_NEBL;DV_NEB];
% SW (nz,1)

 J_SWU=(nz-2)*nx + 1 ;
 J_SWUR=(nz-2)*nx +2;
 J_SWC=(nz-1)*nx + 1;
 J_SWR=(nz-1)*nx + 2;

  V_SWU=-((1-b).*idh2) .* pml.A(nz-1,1) + (b.*idh2) .* pml.By(nz-1,1) + (d./4.0) .* pml.C(nz-1,1) .* k2(nz-1,1);
  V_SWUR=((1-b).*idh2./2.0) .* (pml.Ax(nz-1,1) + pml.By(nz-1,2))  + (e./4.0) .* pml.C(nz-1,2) .* k2(nz-1,2);
  V_SWC=(1-d-e) .* pml.C(nz,1) .* k2(nz,1)- (2.0.*b.*idh2) .* (pml.A(nz,1)+pml.B(nz,1));
  V_SWR=(b.*idh2) .* pml.Ax(nz,1) - ((1-b).*idh2) .* pml.B(nz,2) + (d./4.0) .* pml.C(nz,2) .* k2(nz,2);
  % Derivative matrix

  DV_SWU= (d./4.0) .* pml.C(nz-1,1) .* dk2(nz-1,1);
  DV_SWUR= (e./4.0) .* pml.C(nz-1,2) .* dk2(nz-1,2);
  DV_SWC=(1-d-e) .* pml.C(nz,1) .* dk2(nz,1);
  DV_SWR= (d./4.0) .* pml.C(nz,2) .* dk2(nz,2);
  %
I_SW=repmat(J_SWC,index_CR,1) ;
J_SW=[ J_SWU; J_SWUR;J_SWC; J_SWR];
V_SW=[V_SWU; V_SWUR; V_SWC;  V_SWR];
DV_SW=[DV_SWU; DV_SWUR; DV_SWC; DV_SWR];
 % SE (nz,nx)

J_SEUL=(nz-2)*nx+ (nx-1);
J_SEU=(nz-2)*nx + nx;
J_SEL=(nz-1)*nx + nx-1;
J_SEC=(nz-1)*nx + nx;

V_SEUL=((1-b).*idh2./2.0) .* (pml.Ax(nz-1,nx-1) + pml.By(nz-1,nx-1))  + (e./4.0) .* pml.C(nz-1,nx-1) .* k2(nz-1,nx-1);
V_SEU=-((1-b).*idh2) .* pml.A(nz-1,nx) + (b.*idh2) .* pml.By(nz-1,nx) + (d./4.0) .* pml.C(nz-1,nx) .* k2(nz-1,nx);
V_SEL=(b.*idh2) .* pml.Ax(nz,nx-1) - ((1-b).*idh2) .* pml.B(nz,nx-1)  + (d./4.0) .* pml.C(nz,nx-1) .* k2(nz,nx-1);
V_SEC=(1-d-e) .* pml.C(nz,nx) .* k2(nz,nx)- (2.0.*b.*idh2) .* (pml.A(nz,nx)+pml.B(nz,nx));
 % Derivative matrix
 DV_SEUL=(e./4.0) .* pml.C(nz-1,nx-1) .* dk2(nz-1,nx-1);
DV_SEU= (d./4.0) .* pml.C(nz-1,nx) .* dk2(nz-1,nx);
DV_SEL= (d./4.0) .* pml.C(nz,nx-1) .* dk2(nz,nx-1);
DV_SEC=(1-d-e) .* pml.C(nz,nx) .* dk2(nz,nx);
%
I_SE=repmat(J_SEC,index_CR,1) ;
J_SE=[ J_SEUL; J_SEU;J_SEL; J_SEC];
V_SE=[V_SEUL; V_SEU;V_SEL; V_SEC] ;
DV_SE=[DV_SEUL; DV_SEU;DV_SEL; DV_SEC] ;
% Combine all corners
Crs.I=[I_NW;I_NE;I_SW;I_SE];
Crs.J=[J_NW;J_NE;J_SW;J_SE];
Crs.V=double([V_NW;V_NE;V_SW;V_SE]);
Crs.DV=double([DV_NW;DV_NE;DV_SW;DV_SE]);
 Crs.I=[I_NW;I_NE;I_SW;I_SE];
 Crs.J=[J_NW;J_NE;J_SW;J_SE];
 Crs.V=[V_NW;V_NE;V_SW;V_SE];

end






