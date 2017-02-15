% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

%Move Constrained DOF to RHS
% fdtilda: Fd - the force due to the prescribed displacements 
% Basically doing static condensation 
% And then solve.....


Fdtilda = zeros(neq,1);
rhs= zeros(neq,1);
for i = 1:neq
    
    rhs = zeros(neq,1);
    for j = 1:nieq
       rhs(i) = rhs(i) + Kdf(i,j)*ModelDc(j);
    end
    Fdtilda(i) = Fd(i) - rhs(i);  
    
end

% Solving the Linear System : 

Mstar = Mdd/beta/delt^2+(1+alpha)*Kdd  ;
R = (1+alpha)*(FEXT(:,n+1)-F_bar_int) - alpha*(FEXT(:,n)-IntF_store(:,n))- Mdd * acc(:,storej);

ModelDx = mldivide(Mstar,R);

% Updating the displacement iterate (corrector): 

dis(:,storej+1)  = dis(:,storej) + ModelDx;

% Updating the velocity and acceleration iterate (corrector) :
acc(:,storej+1) = 1/beta/delt^2*(dis(:,storej+1)-dn(:,n) - delt * vn(:,n))-(1/(2*beta)-1)*an(:,n);
vel(:,storej+1)  = vn(:,n) + delt*((1-gamma)*an(:,n) + gamma*acc(:,storej+1));

% Updating the Residual for the iterate : 
res = norm(abs(R));