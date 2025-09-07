% This function adds PML layer to all directions of initial model domain.
function [ vs, rho, int ] = Model_Extension( int, vs, rho )



 % New vs and rho
   vs=  [vs(1,1)*ones(int.pml_thc) , repmat(vs(1,:),int.pml_thc,1) , vs(1,end)*ones(int.pml_thc) ;...
         repmat(vs(:,1),1,int.pml_thc) , vs , repmat(vs(:,end),1,int.pml_thc) ; ...
        vs(end,1)*ones(int.pml_thc) , repmat(vs(end,:),int.pml_thc,1) , vs(end,end)*ones(int.pml_thc)];

  rho=  [rho(1,1)*ones(int.pml_thc) , repmat(rho(1,:),int.pml_thc,1) , rho(1,end)*ones(int.pml_thc) ;...
         repmat(rho(:,1),1,int.pml_thc) , rho , repmat(rho(:,end),1,int.pml_thc) ; ...
         rho(end,1)*ones(int.pml_thc) , repmat(rho(end,:),int.pml_thc,1) , rho(end,end)*ones(int.pml_thc)];

     % New nx,nz and nxz
      int.nx = int.nx + 2 * int.pml_thc;
      int.nz = int.nz + 2 * int.pml_thc;
      int.nxz = int.nx * int.nz;
%
end

