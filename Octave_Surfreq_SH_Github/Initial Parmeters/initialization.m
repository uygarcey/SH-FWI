% This function returns "int" structure which contains all parameters needed for Forward and Inverse solutions
%The idea of putting initial parameters into a structure was got from Germania repository

function[ int] = initialization(fmax,G,f,vs)

   int.stencil=9; %9-point stencil used for finite difference
     int.sd=3;

    [int.nz,int.nx]=size(vs);
     int.mf =f; %Current frequency. (All calculations are performed  frequency by frequency)
    int.fc=45; %Dominant frequency(... for calculation of Ricker wavelet).
% Numerical Dispersion Condition
    int.dx=min(min(vs))/(G*fmax);

    int.nxz = int.nx * int.nz;
    int.omega = 2.0*pi*int.mf;
%
% Display
     int.x = (1:int.dx:int.nx);
    int.z =(1:int.dx:int.nz);

% Parameterization
%int.parameter='Squared_Slowness';
 int.parameter='Velocity';

% Pml parameters
    int.pml_thc =55;
    int.a0_pml = 1.79;
    int.b0_pml=1.885;
    int.n_pml=2;

    int.m2_pml=2;
    int.m3_pml=3;

  % Receviers
     int.rcvdx=1;%interval between each receiver
     int.first_rcv=1; %First receiver point
     int.last_rcv=int.nx; %Last receiver point
 int.rcv1=int.pml_thc+int.first_rcv; %First receiver point for extended model domain (After adding PML layers to initial model domain)
int.rcvN=int.pml_thc+int.last_rcv; % Last receiver pont first reciver point for extended model domain
   int.rcvx=int.rcv1:int.rcvdx:int.rcvN;% Vector contains all receiver points
    int.nrcv=length(int.rcvx);%Number of receivers
     int.rcvz= int.pml_thc+1;%All receivers are on the surfaces(...just begining of the upper PML layer)
% Sources
int.srcdx=5; %interval between each source point
int.first_shot=2;%First source point
int.last_shot=98;%Last source point
int.shot1=int.pml_thc+int.first_shot;%First source point for extended model domain (After adding PML layers to initial model domain)
int.shotN=int.pml_thc+int.last_shot;%Last source point for extended model domain
int.srcx=int.shot1:int.srcdx:int.shotN;% Vector contains all source points
int.nsrc=length(int.srcx);%Number of sources

int.srcz= int.pml_thc+1;%All sources are  on the surface(...just begining of the upper PML layer)

% Free Surface Condition
    int.surface = 0; %No damping on upper bound of calculation domain
    %int.surface =1; % Apply damping


end





