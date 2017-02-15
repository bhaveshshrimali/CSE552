% Compute strain and stress at integration points
%
% Soonpil Kang
% nen: Maximum Number of nodes per element
% nel: Number of nodes on the current element (3 or 4)
% ndm: Space dimension of the mesh (2)
% ndf: Maximum Number of dof per node (2)
% node:
% ElemFlag: 
% xl:       = local array containing (x,y) coordinates of nodes
%                  forming the element; format is as follows:
%                      Nodes    |        n1  n2  n3  n4
%                      x-coord  |  xl = [x1  x2  x3  x4
%                      y-coord  |        y1  y2  y3  y4];

% 3/2016
% UIUC

for elem = 1:numel
    
    %Determine element parameters
    if nen == 3
        nel = 3;
    elseif nen == 4
        if ix(elem,nen) == 0
            nel = 3;
        else
            nel = 4;
        end
    elseif nen == 6
        nel = 6;
    else
        if ix(elem,nen) == 0
            nel = 6;
        else
            nel = 9;
        end
    end
    nst = nel*ndf;
    
    %Extract element nodal coordinates
    xl = zeros(ndm, nel);
    ElemFlag = zeros(nel, 1);
    for k = 1:nel
        node = ix(elem,k);
        ElemFlag(k) = node;
        for l = 1:ndm
            xl(l,k) = NodeTable(node,l);
        end
    end
    
    %Extract element material properties
    ma = ix(elem,nen+1);
    mateprop = MateT(ma,:);

    %Extract element nodal displacements
    EDOFT = LocToGlobDOF(ElemFlag, NDOFT, nel, ndf);
    
    ul = zeros(ndm,nel);
    for i = 1:nel*ndf
        ndof_index = EDOFT(i);
        if(ndof_index<=neq)
            ul(i) = ModelDx(ndof_index);
        else
            ul(i) = ModelDc(ndof_index-neq);
        end
    end
    
    [strain,stress] = CompStrainStress_Elem_Cee570(xl,ul,mateprop,nel,ndf,PSPS);

    % Plot strain and stress arrays
    fprintf('\nElement\t%1.0f\n\n',elem);
    filename = 'Stress.xlsx';
    xlswrite(filename,stress);
    filename = 'Strain.xlsx';
    xlswrite(filename,strain);
end