%Nonlinear Finite Element Method : Fall 2016

% HW Assignment #1
% Problem #3
%Date: 09/02/2016

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
    Res_temp = abs(double(Fext(i)-Fint));
    del_d_temp = Res_temp/double(subs(slope,dn));
    dn = dn+del_d_temp;
    Fint = double(subs(N,dn));
    
    Res(i) = abs(double(Fext(i)-Fint));
    j =1;
    Res_int = Res(i)
    Resi = Res_int
    while Resi  >  (epsilon)*Res_int 
        j
        Residual(j,i) = Resi;
        s = double(subs(slope,dn));
        del_d = double(Resi/s);
        dn = dn+del_d;        
        F_int(j,1) = double(subs(N,dn));
        Resi = double(Fext(i) - F_int(j,1))
        Residual(j+1,i) = Resi;
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
filename = 'NR.xlsx';
xlswrite(filename,Residual);
Fall = [displacement Internal_force];
filename = 'NRFint.xlsx';
xlswrite(filename,Fall);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


