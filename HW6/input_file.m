% Input File: One Quadrilateral Element Under Axial Load
%
% Copyright (C) Bhavesh Shrimali
%
% This input file should be run prior to executing the FEA_Program routine.
%
% Format of required input:
%
%   numnp:           = number of nodes in the mesh (length(NodeTable))
%
%   numel:           = number of elements in the mesh
%
%   nen:             = maximum number of nodes per element (4)
%
%   PSPS:            = flag for plane stress ('s') or plane strain ('n')
%
%   NodeTable:       = table of mesh nodal coordinates defining the
%                      geometry of the mesh; format of the table is as
%                      follows:
%                          Nodes  |             x-coord  y-coord
%                          n1     |  NodeTable = [x1     y1
%                          n2     |               x2     y2
%                          ...    |               ..     ..
%                          nnumnp |               xnumnp ynumnp];
%
%   ix:              = table of mesh connectivity information, specifying
%                      how nodes are attached to elements and how materials
%                      are assigned to elements; entries in the first nen
%                      columns correspond to the rows of NodeTable
%                      representing the nodes attached to element e;
%                      entries in the last nen+1 column are rows from MateT
%                      signifying the material properties assigned to
%                      element e; format of the table is as follows:
%                          Elements  |         n1    n2    n3    n4   mat
%                          e1        |  ix = [e1n1  e1n2  e1n3  e1n4 e1mat
%                          e2        |        e2n1  e2n2  e2n3  e2n4 e2mat
%                          ...       |         ..    ..    ..    ..   ..
%                          enumel    |        values for element numel   ];
%
%   MateT:           = table of mesh material properties for each distinct
%                      set of material properties; these sets are
%                      referenced by element e by setting the value of
%                      ix(e,nen+1) to the row number of the desired
%                      material set; format of the table is as follows:
%                          Materials  |           E   v   t
%                          mat1       |  MateT = [E1  v1  t1
%                          mat2       |           E2  v2  t2
%                          ...        |           ..  ..  ..];
%
%   BCLIndex:        = list of the number of boundary conditions and loads
%                      applied to the mesh; first entry is the number of
%                      prescribed displacements at nodes; second entry is
%                      the number of nodal forces
%
%   NodeBC:          = table of prescribed nodal displacement boundary
%                      conditions; it contains lists of nodes, the
%                      direction of the displacement prescribed (x=1, y=2),
%                      and the value of the displacement (set 0 for fixed
%                      boundary); the length of the table must match the
%                      entry in BCLIndex(1), otherwise an error will result
%                      if too few conditions are given or extra BCs will be
%                      ignored in the model input module;  format of the
%                      table is as follows:
%                          BCs  |            nodes direction value
%                          bc1  |  NodeBC = [bc1n   bc1dir   bc1u
%                          bc2  |            bc2n   bc2dir   bc2u
%                          ...  |             ..     ..       .. ];
%
%   NodeLoad:        = table of prescribed nodal forces; it contains lists
%                      of nodes, the direction of the force prescribed
%                      (x=1, y=2), and the value of the force; the length
%                      of the table must match the entry in BCLIndex(2),
%                      otherwise an error will result if too few conditions
%                      are given or extra loads will be ignored in the
%                      model input module; format of the table is as
%                      follows:
%                          Loads  |              nodes direction value
%                          P1     |  NodeLoad = [ P1n    P1dir    P1P
%                          P2     |               P2n    P2dir    P2P
%                          ...    |               ..     ..       .. ];
%
%
%           4 -------------- 3       2
%           |                |       | \
%           |                |       |  \
%           |                |       |   \
%           |                |       |    \
%           |                |       |     \
%           |                |       |      \
%           1 -------------- 2       3-------1
%

clear all; close all; clc;

% Mesh Nodal Coordinates
El.type='1x1';        El.loading = 'No Loading';
El.bf = 'False';

rho = 1.0;
endtime = 15;
beta = 0.25; gamma = 0.5;

alpha =-1/3;
do = 0.05;delt = 0.001;del_t = delt;

time = linspace(0,endtime,endtime/delt);
maxiter = 50;

switch El.loading
    case 'No Loading'
        P =0*linspace(0.1,5.94,10);
    case '4(a)Tension'
        P = linspace(0.1,1.46,10);   %FOR 4(A) TENSION
    case '4(b)Shear'
        P = linspace(0.1,1.25,10);   %FOR 4(B) SHEAR
    case '3(b)Comrpession'
        P = linspace(0.1,5.94,10);   %FOR 3(A) TENSION
    case '3(a)Tension'
        P = linspace(0.1,5.84,10);   %FOR 3(B) COMRPESSION
end

switch El.type
    case '1x1'
        %%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        NodeTable = [0 0
            1 0
            0 1
            1 1];    % TO BE RUN ONLY FOR 1x1 Mesh
        
        numnp = length(NodeTable);
        ix = [1 2 4 3 1];     % FOR 1x1 Mesh
        nen = 4;
        numel = length(ix(:,1));
        meshx = 1; meshy = 1;
        
        
        % Mesh Boundary Conditions and Loads
        
        BCLIndex = [meshy+2 0]';
        
        yvect = (1:meshy)'*(meshx+1) + 1;
        xvect = (1:meshy+1)'*(meshx+1);
        
        NodeBC = [1 1 0
            1 2 0
            yvect ones(length(yvect),1) zeros(length(yvect),1)];
        
        % Mesh Material Properties: Instead of Young's Modulus and Poisson's
        % ratio we specify the given parameters 'a' and 'b'
        
        a = 40;
        b = 10;
        thick = 1;
        PSPS = 's';
        Load_incr = 10; Res_0 = 0;
        MateT = [a b thick rho]; tol=1e-10;
        
        switch El.bf
            case 'True'
                ModelDx = [0 0 0 0 0]';       % Initial Displacement
                v_0 = [0 0 0 0 0]';
                v = v_0;
            case 'False'
                ModelDx = [do 0 0 do 0]';       % Initial Displacement
                v_0 = [0 0 0 0 0]';
                v = v_0;
        end
end
Stress_count = zeros(1,3);
Strain_count = zeros(1,3);

for count = 1:length(time)
      
    switch El.bf
        case 'True'
            Fb = [0.1*sin(pi*time(count)) 0]';
        case 'False'
            Fb = [0 0]';
    end
    
    %%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%   
    if count == 1
        lambda = 0;
        iter =1;
        FEA_Program
        
        dn = ModelDx;
        vn = v;
        an=ap;
        En = En_til;
    else
        lambda = 0; iter = 1;
        
        % Predictor Phase:
        ModelDx = dn + del_t*vn+0.25*(del_t)^2*an;
        v = vn+del_t*0.5*an;
       
        % NEWTON RAPHSON
        
        while iter>=1 
            if iter==1
                FEA_Program
                Res_o = Res;
                iter = iter+1;
            else
                FEA_Program
                if Res <= tol*Res_o
                    break;
                end
                iter=iter+1;
                
            end
        end
        dn = ModelDx;
        vn = v;
        an=ap;
        En = 0.5*vn'*(Mdd*vn) + Un;
    end
    
    Energ(count,1) = En;
    
    dN2x(count,1) = Node_U_V(2,1);
    dN4x(count,1) = Node_U_V(4,1);
    dN2y(count,1) = Node_U_V(2,2);
    dN4y(count,1) = Node_U_V(4,2);
    
    % Stress at the desired integration point
    Sigma(count,:) = stress(2,:);
    Epsilon(count,:) = strain(2,:);
end

Output