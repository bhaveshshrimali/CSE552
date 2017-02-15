% Solve Partitioned Finite Element Matrix System
%
% Copyright (C) Arif Masud and Tim Truster
% 7/2009
% UIUC

%Move Constrained DOF to RHS

Fdtilda = zeros(neq,1);

for i = 1:neq
    
    rhs = 0;
    for j = 1:nieq
       rhs = rhs + Kdf(i,j)*ModelDc(j);
    end
    Fdtilda(i) = Fd(i) - rhs;
    
end

%Solve Kd = F

ModelDx = Kdd\Fdtilda; %#ok<NASGU>