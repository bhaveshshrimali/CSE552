% Input File: One Quadrilateral Element Under Axial Load
%
% Copyright (C) Arif Masud and Tim Truster
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
L = 1; H = 1; meshx = 4; meshy = 4;

%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

NodeTable = [0.25*(0:meshx)' zeros(length(0:meshx),1)
            0.25*(0:meshx)' 0.25*ones(length(0:meshx),1)
            0.25*(0:meshx)' 0.5*ones(length(0:meshx),1)
            0.25*(0:meshx)' 0.75*ones(length(0:meshx),1)
            0.25*(0:meshx)' ones(length(0:meshx),1)
             ];         % TO BE RUN FOR 4x4 MESH
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%

% NodeTable = [0 0
%              1 0
%              0 1
%              1 1];    % TO BE RUN ONLY FOR 1x1 Mesh 

numnp = length(NodeTable);

% Mesh Element Connectivities
ix = [(1:19)' (2:20)' (7:25)' (6:24)' ones(length(6:24),1)];
ixl = zeros(1,5);
for i = 1:4
    ixl = [ixl;ix(5*(i-1)+1:5*(i-1)+4,:)];
end
ixl = ixl(2:length(ixl(:,1)),:);  % For 4x4 Mesh
ix = ixl;

% ix = [1 2 4 3 1];     % FOR 1x1 Mesh

nen = 4;
numel = length(ix(:,1));
% meshx = 1; meshy = 1;


% Mesh Boundary Conditions and Loads

BCLIndex = [meshy+2 meshy+1]';

yvect = (1:meshy)'*(meshx+1) + 1;
xvect = (1:meshy+1)'*(meshx+1);

NodeBC = [1 1 0
          1 2 0
          yvect ones(length(yvect),1) zeros(length(yvect),1)];
      
% Mesh Material Properties: Instead of Young's Modulus and Poisson's
% ratio we specify the given parameters 'a' and 'b'

a = 40;
b = 60;
thick = 1;
PSPS = 's';
Load_incr = 10; Res_0 = 0;
MateT = [a b thick]; tol=1e-8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%

%Initialization of the Displacements 
ModelDx = zeros(2*numnp-BCLIndex(1),1);
ModelDc = zeros(BCLIndex(1),1); 

P = linspace(0.1,1.46,10);   %FOR 4(A) TENSION
% P = linspace(0.1,1.25,10);  FOR 4(B) SHEAR
% P = linspace(0.1,5.94,10);  FOR 3(B) COMRPESSION
% P = linspace(0.1,5.84,10);  FOR 3(A) TENSION

Stress_count = zeros(1,3);
Strain_count = zeros(1,3);

% fint=zeros(8,1);
% F_bar_int = zeros(8,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%
% NEWTON RAPHSON

for count=1:Load_incr
    NodeLoad = [xvect ones(length(xvect),1) [P(count) 2*P(count) 2*P(count) 2*P(count) P(count)]'];   % For 4x4 Mesh
%     NodeLoad = [2 1 -1*P(count)
%                 4 1 -1*P(count)];  % For 1x1 Mesh
%    NodeLoad = [xvect 2*ones(length(xvect),1) [P(count) 2*P(count)
%    2*P(count) 2*P(count) P(count)]'];   % For 4x4 Mesh Shear Loading
    
iter = 1;
    
    while iter >= 1
    if iter ==1
        FEA_Program;
        Res_0 = R;
        iter = iter+1;
%     if abs(strain(3,1)) >= 0.05
%         break
%     end
    else
        FEA_Program
        if norm(R) < tol*norm(Res_0)
            break
        end
        iter = iter+1;

%         if abs(strain(3,1)) >= 0.05 
%         break
%         end
    end
    end
%     if strain(4,1) >= 0.05
%         break
%     end
    Niter(count) = iter;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Counting the total number of load steps
total_count = 0.25*(length(Stress_count)-1);
Stress_count = Stress_count(1:length(Stress_count(:,3)),:);
Strain_count = Strain_count(1:length(Stress_count(:,3)),:);

% Plotting the corresponding stress and strains. First collect the strain
% and stress on the closest integration point (which is 3 in this case), as
% we follow the counter-clockwise numbering of the integration points. 

for i = 1:total_count
    Sigma(i) = Stress_count(4*i-1,3);
    Epsilon(i) = 0.5*Strain_count(4*i-1,3);
end

% Plotting the Stress-Strain Results
figure(1)
plot(Epsilon,Sigma,'LineWidth',2);
hx = xlabel('\bf Strain at Closest Integration Point (Elem 16) ($\varepsilon_{xy}(\frac{1}{\sqrt{3}},\frac{1}{\sqrt{3}})$)','FontWeight','bold','FontSize',18);
set(hx,'Interpreter','latex');
hy = ylabel('\bf Stress at Closest Integration Point (Elem 16) ($\sigma_{xy}(\frac{1}{\sqrt{3}},\frac{1}{\sqrt{3}})$)','FontWeight','bold','FontSize',18);
set(hy,'Interpreter','latex');
leg = legend({'\bf Stress ($\sigma_{xy}$) - Strain($2 \varepsilon_{xy}$) Curve'},'FontWeight','bold','FontSize',18,'Location','northwest');
set(leg,'Interpreter','latex');
t = title('\bf ($\sigma_{xy}$) vs ($2 \varepsilon_{xy}$), in Shear (4b), (16 Elements)','FontWeight','bold','FontSize',18);
set(t,'Interpreter','latex');

%%%%%%%%%%%%%%%% MODIFICATION %%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: This piece of code is used to write data to excel file for record


% Writing the Stress and Strain at Integration Point to the Excel File
% filename = 'Stress4a.xlsx';
% xlswrite(filename,Sigma');
% filename = 'Strain4a.xlsx';
% xlswrite(filename,Epsilon');
% filename='Displacements4a .xlsx';
% xlswrite(filename,Node_U_V);
grid on;

% Node_U_V ;
% Epsilon;
% Sigma;


% Plotting the results
% grid on;
% figure(2)
% plotModelCont(NodeTable, Node_U_V(:,1), ix, numel, nen, 1, 1, 2, ' Displacement Contour - u_x');
% plotModelCont(NodeTable, Node_U_V(:,2), ix, numel, nen, 1, 2, 2, 'Displacement Contour - u_y')
% figure(3)
% plotModel(NodeTable2, ix, numel, nen, 2, 1, 1, 'Deformed Configuration', 'y', 'y')
% plotModel(NodeTable, ix, numel, nen, 2, 1, 1, 'Deformed Configuration', 'n', 'n')