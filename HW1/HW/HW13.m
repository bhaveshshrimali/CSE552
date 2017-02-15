%Nonlinear Finite Element Method : Fall 2016

% HW Assignment #1
% Problem #3
%Date: 09/02/2016

clear all; close all; clc;
disp ('%-------------------------------------------------------------------------%')
disp ('%         CEE 576: Nonlinear Finite Element Method                        %')
disp ('%      HW Assignment #1  (Modified Newton Raphson Code ), Problem #3      %')
disp ('%-------------------------------------------------------------------------%')

syms x;
N = (0.19*x^3 -2*x^2 + 6*x)*exp(0.02*x);
slope = diff(N,x);
Fint = double(subs(N,0));
Fext = 0.5:0.5:5.5;
count = 1;
epsilon = 1e-12;
dn = 0;
y = linspace(0,7,100); 
displacement = zeros (length(Fext),1);
Internal_force = zeros (length(Fext),1);

for i=1:length(Fext)
    Res(i) = abs(double(Fext(i)-Fint));
    j =1;
    Res_int = Res(i);
    Resi = Res_int;
               i
    while Resi  >  (epsilon * Res_int) 
%         Residual(j,i) = Resi;
        s = double(subs(slope,dn));
        del_d = double(Resi/s);
        dn = dn+del_d;        
        F_int(j,1) = double(subs(N,dn));
        Resi = Fext(i) - F_int(j,1);
        j
        j=j+1;
    end
    displacement(i,1) = dn;
    Internal_force(i,1) = F_int(j-1,1);
    Fint = double(subs(N,dn));
end
N_plot = double(subs(N,y));
plot(y,N_plot,displacement,Internal_force,'o');
grid on;
xlabel('Displacement (d)');
ylabel('Internal Force (F_{int})');
legend('Actual Curve','Obtained Curve');
title('F_{int} vs d (Newton Raphson)');
