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
In_slope = subs(slope,0);
Fint = double(subs(N,0));
Fext = 0.5:0.5:5.5;
count = 1;
epsilon = 1e-12;
d1 = 0;
 
for i=1:length(Fext)
    Res(i) = (double(Fext(i)-Fint));
    j =1;
    Res_int = Res(i);
    Resi = Res_int;
    
    if i ==1
        d_in = d1;
    end
    s = double(subs(slope,d_in));       

    while Resi  >  epsilon * Res_int
        if j==1
            dn = d_in;
        end
        d(count,1) = dn;
        del_d = double(Resi/s);
        dn = dn+del_d;        
        d(count+1,1) = dn+del_d;
        F_int(count,1) = double(subs(N,dn));
        Resi = Fext(i) - F_int(count,1);
        j=j+1;
        count = count+1;
    end
    Resi
    j
    d_in = dn;
    Fint = double(subs(N,d_in));
end
d1 = 0:0.05:7;
N_plot = double(subs(N,d1));
plot(d1,N_plot,d,F_int,'o');
grid on;
xlabel('Displacement (d)');
ylabel('Internal Force (F_{int})');
legend('Actual Curve','Obtained Curve');
title('F_{int} vs d (Newton Raphson)');
