function [alpha] = Wolfe(int,fh,x0,Lb,alpha0,p)
% Simple Wolfe linesearch, adapted from
% Lecture notes:High resolution geophysicsal imaging using full waveform
% modeling and inversion
%L.Metivier,2022
%

c1     = 1e-4;
c2     = 0.9;
g0=Lb.grad;
f0=Lb.norm;
counter=0;
finish= 1;
alpha_r  = 0;
alpha_l  = 0;

alpha=alpha0;
while finish

    if counter < 5
        mt         = x0 + (p*alpha);
         mt = (reshape(mt,int.nx,int.nz))';

        [Lb] = fh(mt);
        ft=Lb.norm;
        gt=Lb.grad;

       counter     = counter + 1;
    else
        alpha = alpha0;
        break;
    end

%
   if ( ft) >= (f0 + c1*alpha*g0'*p)

    alpha_r=alpha;

       alpha=0.5*(alpha_l +alpha_r);

    elseif (gt'*p) <= (c2*g0'*p)
       alpha_l =alpha;
        if alpha_r == 0
            alpha_r=10*alpha;
        else
           alpha=0.5*(alpha_l +alpha_r);
        end

    else
           finish = 0;
    end


end

