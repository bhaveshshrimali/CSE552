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


%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%

%Solve the linearized system
R = Fdtilda - F_bar_int;
del_ModelDx = Kdd\R; 
ModelDx = ModelDx + del_ModelDx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
