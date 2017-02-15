clear all; close all; clc;


syms xi

N1 = xi*(xi-1)/2;
N2 = (1+xi)*(1-xi);
N3 = xi*(xi+1)/2;


n1 = diff(N1,xi);n2 = diff(N2,xi);n3 = diff(N3,xi);

A = [n1;n2;n3]*[n1 n2 n3];
A = simplify(A)

A1 = double(subs(A,xi,-1/(3)^0.5))
A2 = double(subs(A,xi,1/(3)^0.5))

A1+A2


B = [n1;n2;n3]*[N1 N2 N3]

B1 = double(subs(B,xi,-1/(3)^0.5))
B2 = double(subs(B,xi,1/(3)^0.5))

B1+B2


C = [N1;N2;N3]*[N1 N2 N3]

C1 = double(subs(C,xi,-1/(3)^0.5))
C2 = double(subs(C,xi,1/(3)^0.5))

C1+C2

D = [N1;N2;N3]*[n1 n2 n3]

D1 = double(subs(D,xi,-1/(3)^0.5))
D2 = double(subs(D,xi,1/(3)^0.5))

D1+D2

