clear all;close all;clc;
tic

% True Model
%examplevs = matfile('gridvs');
%examplerho = matfile('gridrho');
examplevs = load('Data\gridvs');
examplerho = load('Data\gridrho');
vs_true= examplevs.vs;
rho_true=examplerho.rho;


% Build inital model
[y,x]=size(vs_true);
vs_initial = ones(y,x).*500;
rho_initial=rho_true;

% Number of nodes per wavelnght
G=10;

tic
 f=45;

fmx=80;
[int] = initialization(fmx,G,f,vs_true);

int_initial=int;
if strcmp(int.parameter,'Velocity')
   model.true = vs_true;
model.initial = vs_initial;
elseif strcmp(int.parameter,'Squared_Slowness')
 model.true = 1./vs_true.^2;
model.initial = 1./vs_initial.^2;
end


[model.true,rho_true,int] =  Model_Extension( int,model.true,rho_true);

[model.initial,rho_initial,int] =  Model_Extension( int_initial, model.initial, rho_initial );

[Smp] = Extract(int);

m0=model.initial;
rho0=rho_initial;

[a0,b0]=size(m0);
% Frequencies were selected according to Sirgue ve Pratt (2004)
   freq=[8 15.0000   20.6897   28.5375   39.3621   54.2926   74.8864 80 ];

 tol=1e-3;
c=0;

% Convert initial model matrix to a vector
x0=reshape(m0,[a0*b0,1]);
% Function handel
fh=@(m0)func(int,freq,model.true, rho_true,Smp,m0,rho_initial);
m=50; % stored vector pairs
iter=1;
maxit=125; % maximum iteration
n = length(x0);
sk = zeros(n,m);
yk = zeros(n,m);





 while (iter < maxit) %(normest(g) > tol)&&  %convergence criteria will be change

   if iter==1
      [Lb0]=fh(m0);
       alpha0=1/2;
    [alpha]= Wolfe(int,fh,x0,Lb0,alpha0,Lb0.grad);
        x1 = x0 + (Lb0.grad*alpha);
        %convert vector x1 to a matrix
  m1 = (reshape(x1,int.nx,int.nz))';
 [Lb1]=fh(m1);

   else
       alpha0=1;
   end



sk0 = x1-x0;
    yk0 = (Lb1.grad- Lb0.grad);

    scaling = sk0'*yk0/(yk0'*yk0);

      p = zeros(length(Lb0.grad),1);
% Calculate Hessian by two-loop recursion. This part is adapted from optLBFGS-master
   if (iter<=m)

        sk(:,iter) = sk0;
        yk(:,iter) = yk0;
     p = -H(Lb1.grad,sk(:,1:iter),yk(:,1:iter),scaling);

      elseif (iter>m)
        sk(:,1:(m-1))=sk(:,2:m);
        yk(:,1:(m-1))=yk(:,2:m);
        sk(:,m) = sk0;
        yk(:,m) = yk0;
        p = -H(Lb1.grad,sk,yk,scaling);

    end


 % Line search

    [alpha_new]= Wolfe(int,fh,x1,Lb1,alpha0,p);

      x0 = x1;
      Lb0 =Lb1;
   % Update model
      x1 = x1 +(p*alpha_new);

 m1 = (reshape(x1,int.nx,int.nz))';
 % Calculate new gradient and norm
[Lb1]=fh(m1);
disp('Non-normalized data residual')
norm=Lb1.norm
 disp('Iteration')
   iter = iter + 1

end
% Smooth the model

    H = fspecial('disk',4);

    smooth =imfilter(m1 ,H,'replicate');
   m0= smooth;
 smooth =smooth(1+int.pml_thc:int.nz-int.pml_thc,1+int.pml_thc:int.nx-int.pml_thc);
 srcdx=int.srcdx*int.dx;

figure(1)

imagesc( int.x.*int.dx, int.z.*int.dx, smooth )
colormap  (jet)
caxis ([400 600])
xlabel('Shots (m)'); ylabel('Depths(m)');
   x=1:srcdx:49;
set(gca,'xtick',x, 'xaxisLocation','top')
c = colorbar;
title(c,'Vs (m/s)')
figure(2)
 imagesc( int.x.*int.dx, int.z.*int.dx,vs_true )
colormap  (jet)
xlabel('Shots (m)'); ylabel('Depths(m)');
caxis([400 600])
     x=1:srcdx:49;
set(gca,'xtick',x,'xaxisLocation','top')
c = colorbar;
title(c,'Vs (m/s)')
%
% axis equal;
% axis tight;
%end


