%This function returns Sds structure which contains side parts of the coefficient matrix with row and column numbers of them.
%
function[Sds]= Sides(int,pml)
int.sd=3;
% Non-zero node numbers
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Corner Nodes
% I indicates rows, J indicates columns, V indicates values


b=pml.b;
d=pml.d;
e=pml.e;
k2=pml.k2;
idh2=pml.idh2;
dk2=pml.dk2;
nx=int.nx;
nz=int.nz;

% W
index_LR=(nz-2)*(int.stencil-int.sd);
index_TB=(nx-2)*(int.stencil-int.sd);
zz=(2:nz-1)'; xx=1;

    J_WU=(zz-2)*nx + 1;
    J_WUR=(zz-2)*nx + 2;
    J_WC=(zz-1)*nx + 1;
    J_WR=(zz-1)*nx + 2;
    J_WB=(zz)*nx+1;
    J_WBR=(zz)*nx+2;


   V_WU=-((1-b).*idh2) .* pml.A(zz-1,xx) + (b.*idh2) .* pml.By(zz-1,xx) + (d./4.0) .* pml.C(zz-1,xx) .* k2(zz-1,xx);
   V_WUR=((1-b).*idh2./2.0) .* (pml.Ax(zz-1,xx) + pml.By(zz-1,xx+1))  + (e./4.0) .* pml.C(zz-1,xx+1) .* k2(zz-1,xx+1);
   V_WC=(1-d-e) .* pml.C(zz,xx) .* k2(zz,xx)- (2.0.*b.*idh2) .* (pml.A(zz,xx)+pml.B(zz,xx));
   V_WR=(b.*idh2) .* pml.Ax(zz,xx) - ((1-b).*idh2) .* pml.B(zz,xx+1) + (d./4.0) .* pml.C(zz,xx+1) .* k2(zz,xx+1);
   V_WB=-((1-b).*idh2) .* pml.A(zz+1,xx) + (b.*idh2) .* pml.By(zz,xx) + (d./4.0) .* pml.C(zz+1,xx) .* k2(zz+1,xx);
   V_WBR=((1-b).*idh2./2.0) .* (pml.Ax(zz+1,xx) + pml.By(zz,xx+1)) + (e./4.0) .* pml.C(zz+1,xx+1) .* k2(zz+1,xx+1);
 % Derivative matrix
  DV_WU= (d./4.0) .* pml.C(zz-1,xx) .* dk2(zz-1,xx);
   DV_WUR= (e./4.0) .* pml.C(zz-1,xx+1) .* dk2(zz-1,xx+1);
   DV_WC=(1-d-e) .* pml.C(zz,xx) .* dk2(zz,xx);
   DV_WR= (d./4.0) .* pml.C(zz,xx+1) .* dk2(zz,xx+1);
   DV_WB= (d./4.0) .* pml.C(zz+1,xx) .* dk2(zz+1,xx);
   DV_WBR= (e./4.0) .* pml.C(zz+1,xx+1) .* dk2(zz+1,xx+1);
   %
 I_W=reshape(repmat(J_WC,1,int.stencil-int.sd)',index_LR,1);
 J_W=reshape([J_WU, J_WUR,J_WC,J_WR,J_WB,J_WBR]',index_LR,1);
V_W=reshape([V_WU, V_WUR,V_WC,V_WR,V_WB,V_WBR].',index_LR,1) ;
DV_W=reshape([DV_WU, DV_WUR,DV_WC,DV_WR,DV_WB,DV_WBR].',index_LR,1) ;
% E
 xx=nx; %zz=(2:nz-1)';
J_EUL=(zz-2)*nx+ (nx-1);
  J_EL=(zz-1)*nx + nx-1;
  J_EU=(zz-2)*nx + nx;
  J_EC=(zz-1)*nx + nx;
J_EBL=(zz)*nx + (nx-1);
J_EB=(zz)*nx+nx;

V_EUL=((1-b).*idh2./2.0) .* (pml.Ax(zz-1,xx-1) + pml.By(zz-1,xx-1))  + (e./4.0) .* pml.C(zz-1,xx-1) .* k2(zz-1,xx-1);
V_EL=(b.*idh2) .* pml.Ax(zz,xx-1) - ((1-b).*idh2) .* pml.B(zz,xx-1)  + (d./4.0) .* pml.C(zz,xx-1) .* k2(zz,xx-1);
V_EU= -((1-b).*idh2) .* pml.A(zz-1,xx) + (b.*idh2) .* pml.By(zz-1,xx) + (d./4.0) .* pml.C(zz-1,xx) .* k2(zz-1,xx);
V_EC=(1-d-e) .* pml.C(zz,xx) .* k2(zz,xx)- (2.0.*b.*idh2) .* (pml.A(zz,xx)+pml.B(zz,xx));
V_EBL=((1-b).*idh2./2.0) .* (pml.Ax(zz+1,xx-1) + pml.By(zz,xx-1))+ (e./4.0) .* pml.C(zz+1,xx-1) .* k2(zz+1,xx-1);
V_EB=-((1-b).*idh2) .* pml.A(zz+1,xx) + (b.*idh2) .* pml.By(zz,xx) + (d./4.0) .* pml.C(zz+1,xx) .* k2(zz+1,xx);
 % Derivative matrix
 DV_EUL= (e./4.0) .* pml.C(zz-1,xx-1) .* dk2(zz-1,xx-1);
DV_EL= (d./4.0) .* pml.C(zz,xx-1) .* dk2(zz,xx-1);
DV_EU= (d./4.0) .* pml.C(zz-1,xx) .* dk2(zz-1,xx);
DV_EC=(1-d-e) .* pml.C(zz,xx) .* dk2(zz,xx);
DV_EBL= (e./4.0) .* pml.C(zz+1,xx-1) .* dk2(zz+1,xx-1);
DV_EB=(d./4.0) .* pml.C(zz+1,xx) .* dk2(zz+1,xx);
%
I_E=reshape(repmat(J_EC,1,int.stencil-int.sd)',index_LR,1);
J_E=reshape([J_EUL,J_EL,J_EU,J_EC,J_EBL,J_EB]',index_LR,1);
V_E=reshape([V_EUL,V_EL,V_EU,V_EC,V_EBL,V_EB].',index_LR,1);
%
DV_E=reshape([DV_EUL,DV_EL,DV_EU,DV_EC,DV_EBL,DV_EB].',index_LR,1);


% N

xx=(2:nx-1)'; zz=1;
      J_NL= xx-1;
       J_NC= xx;
        J_NR= xx+1;
       J_NB=nx + xx;
       J_NBL=nx + (xx-1);
       J_NBR=nx+(xx+1);

V_NL=(b.*idh2) .* pml.Ax(zz,xx-1) - ((1-b).*idh2) .* pml.B(zz,xx-1)  + (d./4.0) .* pml.C(zz,xx-1) .* k2(zz,xx-1);
 V_NC=(1-d-e) .* pml.C(zz,xx) .* k2(zz,xx)- (2.0.*b.*idh2) .* (pml.A(zz,xx)+pml.B(zz,xx));
V_NR=(b.*idh2) .* pml.Ax(zz,xx) - ((1-b).*idh2) .* pml.B(zz,xx+1) + (d./4.0) .* pml.C(zz,xx+1) .* k2(zz,xx+1);
V_NB=-((1-b).*idh2) .* pml.A(zz+1,xx) + (b.*idh2) .* pml.By(zz,xx) + (d./4.0) .* pml.C(zz+1,xx) .* k2(zz+1,xx);
V_NBL=((1-b).*idh2./2.0) .* (pml.Ax(zz+1,xx-1) + pml.By(zz,xx-1))+ (e./4.0) .* pml.C(zz+1,xx-1) .* k2(zz+1,xx-1);
V_NBR=((1-b).*idh2./2.0) .* (pml.Ax(zz+1,xx) + pml.By(zz,xx+1)) + (e./4.0) .* pml.C(zz+1,xx+1) .* k2(zz+1,xx+1);
% Derivative matrix
DV_NL= (d./4.0) .* pml.C(zz,xx-1) .* dk2(zz,xx-1);
 DV_NC=(1-d-e) .* pml.C(zz,xx) .* dk2(zz,xx);
DV_NR= (d./4.0) .* pml.C(zz,xx+1) .* dk2(zz,xx+1);
DV_NB= (d./4.0) .* pml.C(zz+1,xx) .* dk2(zz+1,xx);
DV_NBL= (e./4.0) .* pml.C(zz+1,xx-1) .* dk2(zz+1,xx-1);
DV_NBR= (e./4.0) .* pml.C(zz+1,xx+1) .* dk2(zz+1,xx+1);
%
I_N=reshape(repmat( J_NC,1,int.stencil-int.sd)',index_TB,1);
J_N=reshape([J_NL,J_NC,J_NR,J_NB,J_NBL,J_NBR]',index_TB,1);
V_N=reshape([V_NL',V_NC',V_NR',V_NB',V_NBL',V_NBR']',index_TB,1);
%
DV_N=reshape([DV_NL',DV_NC',DV_NR',DV_NB',DV_NBL',DV_NBR']',index_TB,1);
% S

zz=nz;
%xx=(2:nx-1)';
      J_SL=(zz-1)*nx + xx-1;
       J_SU=(zz-2)*nx + xx;
       J_SC=(zz-1)*nx + xx;
       J_SR=(zz-1)*nx + xx+1;
       J_SUL=(zz-2)*nx+ (xx-1);
       J_SUR=(zz-2)*nx + (xx+1);

V_SL=(b.*idh2) .* pml.Ax(zz,xx-1) - ((1-b).*idh2) .* pml.B(zz,xx-1)  + (d./4.0) .* pml.C(zz,xx-1) .* k2(zz,xx-1);
V_SU=-((1-b).*idh2) .* pml.A(zz-1,xx) + (b.*idh2) .* pml.By(zz-1,xx) + (d./4.0) .* pml.C(zz-1,xx) .* k2(zz-1,xx);
V_SC=(1-d-e) .* pml.C(zz,xx) .* k2(zz,xx)- (2.0.*b.*idh2) .* (pml.A(zz,xx)+pml.B(zz,xx));
V_SR= (b.*idh2) .* pml.Ax(zz,xx) - ((1-b).*idh2) .* pml.B(zz,xx+1) + (d./4.0) .* pml.C(zz,xx+1) .* k2(zz,xx+1);
V_SUL=((1-b).*idh2./2.0) .* (pml.Ax(zz-1,xx-1) + pml.By(zz-1,xx-1))  + (e./4.0) .* pml.C(zz-1,xx-1) .* k2(zz-1,xx-1);
V_SUR=((1-b).*idh2./2.0) .* (pml.Ax(zz-1,xx) + pml.By(zz-1,xx+1))  + (e./4.0) .* pml.C(zz-1,xx+1) .* k2(zz-1,xx+1);
  % Derivative matrix
DV_SL= (d./4.0) .* pml.C(zz,xx-1) .* dk2(zz,xx-1);
DV_SU= (d./4.0) .* pml.C(zz-1,xx) .* dk2(zz-1,xx);
DV_SC=(1-d-e) .* pml.C(zz,xx) .* dk2(zz,xx);
DV_SR=  (d./4.0) .* pml.C(zz,xx+1) .* dk2(zz,xx+1);
DV_SUL= (e./4.0) .* pml.C(zz-1,xx-1) .* dk2(zz-1,xx-1);
DV_SUR=(e./4.0) .* pml.C(zz-1,xx+1) .* dk2(zz-1,xx+1);
 %
I_S=reshape(repmat( J_SC,1,int.stencil-int.sd)',index_TB,1);
J_S=reshape([J_SL,J_SU,J_SC,J_SR,J_SUL,J_SUR]',index_TB,1);
V_S=reshape([V_SL',V_SU',V_SC',V_SR',V_SUL',V_SUR']',index_TB,1);
%
DV_S=reshape([DV_SL',DV_SU',DV_SC',DV_SR',DV_SUL',DV_SUR']',index_TB,1);
% Combine all sides
Sds.I=[I_W;I_N;I_S;I_E];
Sds.J=[J_W;J_N;J_S;J_E];
Sds.V=double([V_W;V_N;V_S;V_E]);
Sds.DV=double([DV_W;DV_N;DV_S;DV_E]);
Sds.I=[I_W;I_N;I_S;I_E];
Sds.J=[J_W;J_N;J_S;J_E];
Sds.V=[V_W;V_N;V_S;V_E];
end








