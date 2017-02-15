% Input File: Patch Test For Shear Loading
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
% The following numbering convention is used for 4-node quadrilateral
% elements:
%
%           4 -------------- 3
%           |                |
%           |                |
%           |                |
%           |                |
%           |                |
%           |                |
%           1 -------------- 2
%

clc; clear all; close all;
meshx = [10 20 40];
meshy = [2 4 8];
H1_all = [0 0 0];
L2_all = [0 0 0];
for i=1:length(meshx)
    tic
    % Arbitrary data for assistance in defining the mesh
    UINsum=3;
    L=10+UINsum;
    C=L/4;
    H=2*C;
    P=5*UINsum;
    
    % Mesh Nodal Coordinates
    NodeTable = [1 0 -C
        2 L -C
        3 L C
        4 0 C];
    type = 'cart';
    rinc = meshx(i);
    sinc = meshy(i);
    node1 =1;
    elmt1 = 1;
    mat = 1;
    rskip = 0;
    btype = 0;
    nen = 4;
    [NodeTable,ix,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,NodeTable,nen);
    
    % Compute Force Vector for a Given Problem (iprob)
    iprob = 1;
    
    nf = 1:(sinc+1);
    %location of y
    locy = NodeTable(nf*(rinc+1),2);
    yel = 2*C/sinc;
    
    %Force Vector
    syms yy;
    sigma_xy = 3*P*(C^2-yy^2)/4/C^3;
    sigma_xx = 3*P*(0-L)*yy/2/C^3;
    for i =1:sinc
        f1y(i) = int(sigma_xy*(locy(i+1)-yy)/yel,yy,locy(i),locy(i+1));
        f2y(i) = int(sigma_xy*(yy-locy(i))/yel,yy,locy(i),locy(i+1));
        f1x(i) = int(sigma_xx*(locy(i+1)-yy)/yel,yy,locy(i),locy(i+1));
        f2x(i) = int(sigma_xx*(yy-locy(i))/yel,yy,locy(i),locy(i+1));
    end
    
    fy(1) = f1y(1);
    fy(2:sinc)=f1y(2:sinc)+f2y(1:(sinc-1));
    fy(sinc+1) = f2y(sinc);
    fx(1) = f1x(1);
    fx(2:sinc)=f1x(2:sinc)+f2x(1:(sinc-1));
    fx(sinc+1) = f2x(sinc);
    
    %% Mesh Boundary Conditions and Loads
    BCLIndex = [4 3*(sinc+1)]';
    NodeBC = [1 1 0
        12 1 0
        12 2 0
        23 1 0];
    
    NodeLoad = [nf'*(rinc+1) 2*ones(1,length(nf))' -fy';
        ((nf'-1)*(rinc+1))+1 ones(1,length(nf))' -fx';
        ((nf'-1)*(rinc+1))+1 2*ones(1,length(nf))' fy';];
    
    %%
    
    % Mesh Material Properties
    young = 10e5;
    pois = .25;
    thick = 1;
    PSPS = 's';
    density = 0;
    Gy = 0;
    MateT = [young pois thick density Gy];
    
    % Error Estimation
    hx=L/rinc; hy=H/sinc;
    h = hx;
    
    FEA_Program
    Error_Estimation
    
    %Store Values of H1 and L2
    H1_all(i) = H1;
    L2_all(i) = L2;
    toc
end