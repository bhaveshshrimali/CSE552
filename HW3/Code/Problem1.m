%Nonlinear Finite Element Method : Fall 2016

%HW Assignment 1
%Date: 09/02/2016

clear all; close all; clc;
disp ('%-----------------------------------------------------------------------------------------%')
disp ('%         CEE 576: Nonlinear Finite Element Method                                        %')
disp ('%      HW Assignment #1  (Modified Newton Raphson With Arc-Length and without Line Search %')
disp ('%-----------------------------------------------------------------------------------------%')

N = @(d) (3.5*sin(d)+(1/3)*d*cos(d)+2*d)*exp(0.02*d);
b=0.5;
del_a = 0.25; fext =1;
tol = 1e-12;
max_iter = 1000;
max_step = 200;
dN = @(d) 0.02*(3.5*sin(d)+1/3*d*cos(d)+2*d)*exp(0.02*d)+...
    (23/6*cos(d)-1/3*d*sin(d)+2)*exp(0.02*d);
n = 200;

% Initial Value of the Tangent Stiffness
K0 = dN(0);
q_til = fext/K0;
c = (1-b)*(K0);
epsilon_p = 1e-12;
epsilon_f = 1e-12;

% Initializing the parameters required in the Computation
displacement = zeros(n,1);
F_ext = zeros(n,1);
K = zeros(n,1);
R = zeros(n,1);
rf = zeros(n,1);
Ni = zeros(n,1);

for n = 1:max_step
    K(n) =dN(displacement(n));
    q = fext/K(n);
    f_trial = (c*q*K0*q + b)^0.5;
    if n ==1
        sign_lambda = sign(K(n));
    elseif n == 2
        if sign(K(n))~=sign(K(n-1))
            sign_lambda = -sign_lambda;
        end
    else
        if (sign(K(n))~=sign(K(n-1)) && abs(K(n-1)) < abs(K(n-2)))
            sign_lambda = -sign_lambda;
        end
    end
    d_lambda = sign_lambda*del_a/f_trial;
    del_d = d_lambda*q;
    
    F_ext(n+1) = F_ext(n) + d_lambda*fext;
    displacement(n+1) = displacement(n)+del_d;
    
    R(1) = F_ext(n+1)-N(displacement(n));
    f(1) = sqrt(c*del_d^2*K0 + b*d_lambda^2);
    r(1) = del_a - f(1);
    
    for i=1:max_iter
        df_ddel_d = c/f(i)*K0*del_d;
        D_d_bar = R(i)/K(n);
        q_bar = fext/K(n);
        D_lambda = (r(i) - df_ddel_d*D_d_bar)/(df_ddel_d*q_bar +  b/f(i)*d_lambda);
        D_d = D_d_bar + D_lambda*q_bar; 
        
        displacement(n+1) = displacement(n+1)+D_d;
        F_ext(n+1) = F_ext(n+1) + D_lambda*(fext);
        
        d_lambda = (F_ext(n+1) - F_ext(n))/fext;
        del_d = displacement(n+1) - displacement(n);
        
        R(i+1) = F_ext(n+1) - N(displacement(n+1));
        f(i+1) = (c*del_d'*diag(K0)*del_d + b*d_lambda^2)^0.5;
        r(i+1) = del_a - f(i+1);
        
        R(n+1) = R(i+1);
        rf(n+1) = r(i+1);
        
        % convergence
        if (abs(R(i+1)) < abs(epsilon_p*R(1))) && ((r(i+1)) < (epsilon_f*del_a))
            Ni(n+1) = i;
            break
        end
    end
    
end
figure(1)
ezplot(N,[0,6.0]);
hold on;
scatter(displacement,F_ext,'MarkerFaceColor',[0 0 0],'LineWidth',2);
xlabel(' Displacement(d)','FontSize',16,'FontWeight','bold');
ylabel(' Force (F)','FontSize',16,'FontWeight','bold');
title('Modified Newton Raphson with No Line Search','FontSize',16,'FontWeight','bold');
legend({'Original Function','Modified NR: No Line Search'},'Location','northwest','FontSize',14,'FontWeight','bold');
grid on;
ylim([0 14]);
