function [strain,stress] = CompStrainStress_Elem_Cee570(xl,ul,mateprop,nel,ndf,PSPS)
%
% Subroutine to compute strain and stress for linear
% 2-dimensional elasticity element. Element currently supports bilinear
% quadrilateral elements with the following node and shape function
% labelling scheme:
%
%  (-1, 1)  4 -------------- 3 ( 1, 1)
%           |       s        |
%           |       ^        |
%           |       |        |
%           |       .-> r    |
%           |                |
%           |                |
%  (-1,-1)  1 -------------- 2 ( 1,-1)
%
% Element local coordinates (r,s) are defined by a coordinate axis with the
% origin at the center of the element; the corners of the element have
% local coordinate values as shown in the figure.
%
% Definitions for input:
%
%   xl:              = local array containing (x,y) coordinates of nodes
%                      forming the element; format is as follows:
%                          Nodes    |        n1  n2  n3  n4
%                          x-coord  |  xl = [x1  x2  x3  x4
%                          y-coord  |        y1  y2  y3  y4];
%
%   mateprop:        = vector of material properties:
%                          mateprop = [E v t]; 
%                                   = [(Young's Modulus) (Poisson's Ratio)
%                                      (thickness)];
%
%   nel:             = number of nodes on current element (4)
%
%   ndf:             = max number of DOF per node (2)
%
%   ndm:             = space dimension of mesh (2)
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
% Definitions for output:
%
%   strain:          = strain array containing strain components
%                      at integration points:
%                                    xx   yy   xy
%                      int1  strain[ .    .    .   
%                      int2          .    .    .  
%                      int3          .    .    .  
%                      int4          .    .    .  ];
%
%   stress:          = stress array containing stress components
%                      at integration points:
%                                    xx   yy   xy
%                      int1  stress[ .    .    .   
%                      int2          .    .    .  
%                      int3          .    .    .  
%                      int4          .    .    .  ];                    
%
% Definitions of local constants:
%
%   nst:             = size of element arrays (ndf*nel)
%
%

% Set Material Properties
a = mateprop(1);
b = mateprop(2);
thick = mateprop(3);

% Load Guass Integration Points

if nel == 3
    lint = 4;
else
    lint = 4;
end


% Initialize Matrix and Vector

nst = nel*ndf;
ul_elem = reshape(ul,ndf*nel,1);
strain = zeros(lint,3);
stress = zeros(lint,3);
strain_temp = zeros(3,1);
stress_temp = zeros(3,1);

% Loop over integration points
for l = 1:lint

        if nel == 3
            [Wgt,r,s] =  intpntt(l,lint,0);
        else
            [Wgt,r,s] =  intpntq(l,lint,0);
        end

        % Evaluate local basis functions at integration point
        shp = shpl_2d(r,s,nel);

        % Evaluate first derivatives of basis functions at int. point
        [Qxy, Jdet] = shpg_2d(shp,xl,nel);

        % Form B matrix
       if nel == 3
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3)];
        else
        Bmat = [Qxy(1,1) 0        Qxy(1,2) 0        Qxy(1,3) 0        Qxy(1,4) 0 
                0        Qxy(2,1) 0        Qxy(2,2) 0        Qxy(2,3) 0        Qxy(2,4)
                Qxy(2,1) Qxy(1,1) Qxy(2,2) Qxy(1,2) Qxy(2,3) Qxy(1,3) Qxy(2,4) Qxy(1,4)];      
            
        end
        
        % Compute strain
        strain_temp = Bmat*ul_elem;
        
        strain(l,1) = strain_temp(1,1);
        strain(l,2) = strain_temp(2,1);
        strain(l,3) = strain_temp(3,1);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%

        
        % Computing the D Matrix at the updated strain values 
            tr_epsilon = strain(l,1)+strain(l,2);
            tr_epssq = strain(l,1)^2 + strain(l,2)^2 + 0.5*strain(l,3)^2;
            
            kappa = (2+tr_epsilon)/(1+tr_epsilon)^2;
            Dmat(1,1) = a*kappa + 3*b;
            Dmat(1,2) = a*kappa;
            Dmat(1,3) = 0;
            Dmat(2,2) = a*kappa+3*b;
            Dmat(2,3) = Dmat(1,3);
            Dmat(3,3) = 1.5*b;
            Dmat(2,1) = Dmat(1,2);
            Dmat(3,2) = Dmat(2,3);
            Dmat(3,1) = Dmat(1,3); 
        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Compute stress
        stress_temp = (a*log(1+tr_epsilon)+a*(tr_epsilon)/(1+tr_epsilon))*[1 1 0]'+3*b*strain_temp;
        stress(l,1) = stress_temp(1,1);
        stress(l,2) = stress_temp(2,1);
        stress(l,3) = 0.5*stress_temp(3,1);
end
end