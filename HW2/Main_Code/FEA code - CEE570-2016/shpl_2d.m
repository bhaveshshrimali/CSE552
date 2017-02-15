function shp = shpl_2d(r,s,nel)
%
% Copyright (C) Arif Masud and Tim Truster
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

shp = zeros(3,nel);
 
if nel == 3
    
    shp(1,1)= 1;
    shp(2,1)= 0;
    shp(3,1)= r;
 
    shp(1,2)= 0;
    shp(2,2)= 1;
    shp(3,2)= s;
 
    shp(1,3)=-1;
    shp(2,3)=-1;
    shp(3,3)= 1-r-s;
 
 
elseif nel == 4
   
    % shp(i,j). j = denotes the node number,
    % i = 1, derivative w.r.t r
    % i = 2, derivative w.r.t s
    % i = 3, shape funtion
    clear shp
    syms xc eta real
    shp(3,1)=0.25*(1-xc)*(1-eta);
    shp(1,1) = diff(shp(3,1),xc);
    shp(2,1) = diff(shp(3,1),eta);
    
    shp(3,2) = 0.25*(1+xc)*(1-eta);
    shp(1,2) = diff(shp(3,2),xc);
    shp(2,2) = diff(shp(3,2),eta);
    
    shp(3,3) = 0.25*(1+xc)*(1+eta);
    shp(1,3) = diff(shp(3,3),xc);
    shp(2,3) = diff(shp(3,3),eta);
    
    shp(3,4) = 0.25*(1-xc)*(1+eta);
    shp(1,4) = diff(shp(3,4),xc);
    shp(2,4) = diff(shp(3,4),eta);
    
    xc=r;
    eta=s;
    shp=double(subs(shp));
end