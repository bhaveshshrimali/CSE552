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
mu1= mateprop(1);
mu2 = mateprop(2);
kappa = mateprop(3);
thick = mateprop(4);      % Load Guass Integration Points

if nel == 3
    lint = 4;
else
    lint = 4;
end


% Initialize Matrix and Vector

nst = nel*ndf;
ul_elem = reshape(ul,ndf*nel,1);
ul_elem2=[ul_elem(1:2:length(ul_elem)) ul_elem(2:2:length(ul_elem))];

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
        
        
% Compute Cauchy Stress  ----- Mofication beyond the lienar code 
 
        Bmat2=[Qxy(1,1) Qxy(1,2) Qxy(1,3) Qxy(1,4);         %Modification
               Qxy(2,1) Qxy(2,2) Qxy(2,3) Qxy(2,4)];
           
        grad_U=Bmat2*ul_elem2;
             
        F=(grad_U)'+eye(length(grad_U)) ;           
        J=det(F);                                           
        N=[0;1]  ;                                         
        
        FN=F*N;
        I4=dot(FN,FN);
        
        % Cauchy Stress : 
        sigma=mu1/J*(F*F') +2*mu2/J*(I4-1)*(FN*(FN)')+((kappa+mu1)*(J-1)-mu1)*eye(2);  
        
        stress(l,1)=sigma(1,1);
        stress(l,2)=sigma(2,2);
        stress(l,3)=sigma(1,2);
       
        % Compute Green-Lagrange Strain
        
        E=0.5*((F'*F)-eye(2));
        strain(l,1) = E(1,1);
        strain(l,2) = E(2,2);
        strain(l,3) = E(1,2);
end
end