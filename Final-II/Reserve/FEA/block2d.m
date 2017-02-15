function [x,ix,numnp,numel] = block2d(type,rinc,sinc,node1,elmt1,mat,rskip,btype,xl,nen,x,ix)

% Copyright (C) Arif Masud and Tim Truster

if strcmp(type,'cart') == 1

    if nargin < 10 || nargin == 11
        error('At least 10 arguments are required')
    elseif nargin == 10
        x = zeros(0,0);
        ix = zeros(0,0);
    elseif nargin == 12
        %x znd ix are given
    end
    
    xll = zeros(2,9);
    ixl = zeros(9,1);
    for i = 1:length(xl)
        ind = xl(i,1);
        ixl(ind) = ind;
        xll(1,ind) = xl(i,2);
        xll(2,ind) = xl(i,3);
    end
    nr = rinc;
    ns = sinc;

    ni = node1;
    ne = elmt1;
    ma = mat;
    nodinc = rskip;
    ntyp = btype;
    ndm = 2;
    nm = 9;
    nen1 = nen+1;

    dr = 2/nr;
    ds = 2/ns;

    nr = nr + 1;
    ns = ns + 1;

    x = sblkn(nr,ns,xll,ixl,dr,ds,ni,ndm,nodinc,nm,x);
    numnp = length(x);
    [ix,numel] = sblke6(nr,ns,ni,ne,ndm,nen1,nodinc,ntyp,nm,ma,ix);
    x = x';
    ix = ix';
end