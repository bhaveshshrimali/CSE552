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
for i = 1:neq
    rhs = 0;
    for j = 1:nieq
        rhs = rhs + Kdf(i,j)*ModelDc(j);
    end
    Fdtilda(i) = Fd(i) - rhs;
end

Fnp (:,count) = Fdtilda;
FINT (:,count) = F_bar_int;
En_til = 0.5*v'*(Mdd*v) + Un;

if count == 1
    ap = Mdd\(Fdtilda - F_bar_int);
else
    % Energy Conserving Algorithm Setting up the matrix system to solve for...
    % del_ModelDx and del_lambda
    
    A11 = (1+lambda)*((4/del_t^2)*Mdd) + (1+alpha+lambda)*(Kdd);
    A12 = 4*Mdd/delt^2*(ModelDx-dn)+F_bar_int-2*Mdd/delt*vn -0.5*(Fdtilda+Fnp(:,count-1));
    A21 = A12';
    b1 = (1+alpha+0.5*lambda)*Fdtilda - (alpha-0.5*lambda)*Fnp(:,count-1)-alpha*(F_bar_int-FINT(:,count-1))+4*Mdd/delt*vn*(1+0.5*lambda)+Mdd*an-(1+lambda)*(F_bar_int+4*Mdd/delt^2*(ModelDx-dn));
%      b1 = -(1+lambda)*((4/del_t^2)*Mdd*(ModelDx - dn) + F_bar_int) - ...
%             alpha*(F_bar_int - FINT(:,count - 1)) + (1+alpha+lambda/2)*Fdtilda ...
%             -(alpha - lambda/2)*Fnp(:,count - 1) + Mdd*an + 4/del_t*(1+lambda/2)*vn
%     
    b2 = En-En_til +0.5*(ModelDx-dn)'*(Fdtilda+Fnp(:,count-1));
    
    zi = A11\A12 ;
    yi = A11\b1;
    
    del_lambda = (A21*yi - b2)/(A21*zi);
    del_D = yi - del_lambda * zi;
    
    ModelDx = ModelDx +  del_D;
    lambda = lambda + del_lambda;
    
    % Correctors
    ap = 4/del_t^2 * (ModelDx - dn) - an - 4/del_t*vn;
    v = vn + gamma*delt*(ap+an);

    Res = norm(([b1;b2]));
end
  
