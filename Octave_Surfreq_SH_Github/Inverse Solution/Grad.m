% This function calculates gradient for each frquency
function [Lb]=Grad(int,f,vs, rho,S,vsint,rhoint)

[pml]= Pml_Damp_Lie( int, vs, rho,f );

% Build coefficient matrix
[Crs]= Corners(int,pml);
[Sds]= Sides(int,pml);
[Inr]= Interior(int,pml);

I=[Crs.I;Sds.I;Inr.I];
J=[Crs.J;Sds.J;Inr.J];
V=[Crs.V;Sds.V;Inr.V];

%M=sparse2(I,J,V);
M=sparse(I,J,V);
% Solve Linear System and Calculate displacement field for true model
[Rhs] = Rhs_Mul_Shot(int,f);
% [Rhs] = Rhs_Shot(int,f);
%if(int.nsrc < 8)
U_true=mldivide(M,Rhs);
%else
%opts.nthreads = 12;
%opts.grain = 32;
 %U_true = spqr_solve(M,Rhs,opts);
%end
 D_true=S'*U_true;

% Build coefficient matrix
[pml2]= Pml_Damp_Lie( int,vsint,rhoint,f );
[Crs]= Corners(int,pml2);
[Sds]= Sides(int,pml2);
[Inr]= Interior(int,pml2);
I=[Crs.I;Sds.I;Inr.I];
J=[Crs.J;Sds.J;Inr.J];
V=[Crs.V;Sds.V;Inr.V];

%M=sparse2(I,J,V);
M=sparse(I,J,V);

% Derivative of coefficient matrix
DV=[Crs.DV;Sds.DV;Inr.DV];
%DM=sparse2(I,J,DV);
DM=sparse(I,J,DV);

% Solve Linear System and Calculate displacement field for current model
[Rhs] = Rhs_Mul_Shot(int,f);

% [Rhs] = Rhs_Shot(int,f);
%if(int.nsrc < 8)
U_initial=mldivide(M,Rhs);
%else
%opts.nthreads = 12;
%opts.grain = 32;
%U_initial= spqr_solve(M,Rhs,opts);
%end

% Adjoint

 D_initial=S'*U_initial;

 e=(D_true-D_initial);

 Adj_Rhs=S*(conj(e));


 %Calculate adjoint wave field
 %if(int.nsrc < 8)
  V_adjoint=mldivide(M,Adj_Rhs);
%else
%opts.nthreads = 12;
%opts.grain = 32;
%V_adjoint = spqr_solve(M,Adj_Rhs,opts);
% end
 % Gradient

% Zero lag cross correlation of  forward and adjoint wave field
if strcmp(int.parameter,'Velocity')

 G=real(sum((DM'*U_initial).*(V_adjoint),2));

elseif strcmp(int.parameter,'Squared_Slowness')

G=real(sum((DM*U_initial).*(V_adjoint),2));

 end

 Lb.grad=-G;

 Lb.norm = .5*(norm( (D_initial-D_true),'fro')^2);

end




