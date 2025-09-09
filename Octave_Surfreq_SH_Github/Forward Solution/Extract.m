%This function is used to extract wave fields in receiver points. It is a sparse matrix which contains only receiver points.
%Developed by Uygar Ceyhan.
function [S] = Extract(int)
S = sparse(double(int.nz*int.nx),double(int.nrcv));
rcv_index = int.rcvx + (int.rcvz-1)*int.nx; %Take receiver indexes
% row=rcv_index(1:int.nrcv);  % rows
% col=1:int.nrcv;%colums
S(sub2ind([int.nx*int.nz,int.nrcv],rcv_index(1:int.nrcv),1:int.nrcv)) = 1; % It equals 1 in receiver points ;
end



