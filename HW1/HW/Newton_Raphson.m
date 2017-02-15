%%CEE576 - Nonlinear Finite Element
%Assignment#1
%Name: Ahmed Ghareeb
%NetID: ghareeb2
%Newton-Raphson Method
clear; close all; clc;
%Input Parameters:
syms N d real;
N=(0.19*d^3-2*d^2+6*d)*exp(0.02*d);
dN=diff(N,d);
%Tolerance Parameter
eps_s=10^-12;
%Load Value and Required no of Steps
F_ext_P=5.5;
del_F=0.5;
n_step=F_ext_P/del_F;
%Initial Values
dn=0; F_int=0;
%Load Steps
for Step=1:n_step
    F_ext(Step)=Step*F_ext_P/n_step;
    R_o=abs(F_ext(Step)-F_int); R=R_o;       %Residual (0)
    sub_step=0;
    %Loop till Convergence
    while R > eps_s*R_o
            sub_step=sub_step+1;                             %Counter
            K=double(subs(dN,d,dn));                       %Calculate K
            del_d=R/K;                            %Disp. diff.
            dn=dn+del_d;                                    %Update Disp.
            F_int=double(subs(N,dn));
            R=abs(F_ext(Step)-F_int);                    %Calculate Residual
            Residual(sub_step,Step)=R;                     %Store esidual
            if sub_step == 1                               %Update R(0) for first iteration
                R_o = R;
            end
    end
    Nsub(Step)=sub_step;
    D(Step)=dn;  %Store Displacement
    FINT(Step)=F_int;
end   
%Plot
FINT=[0 FINT]; D=[0 D];
plot(D,FINT,'-K'); xlabel('Displacement'); 
ylabel('Force'); title('Force- Displacement Curve');