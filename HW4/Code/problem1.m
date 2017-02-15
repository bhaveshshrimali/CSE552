clear all; close all; clc;

syms z1 z2 z3  lam t real
F = [1+z2 z1 0;
    3*z2 1+3*z1 0;
    0 0 1];
C = simplify(F'*F);
E = 0.5*(C-eye(3))
z1 = 1;z2 = 1;z3 = 0;
C = double(subs(C))
E = double(subs(E))

[A,B]=eig(C)




% F = [1 t; 0.5*t 1]
% a = simplify(inv(F))
% 
% answer = simplify([0 1; 0.5 0]*a)
% 
% D = 0.5*(answer+answer')
% D = double(subs(D,t,0.5))
% Omega = 0.5*(answer-answer')
% Omega = double(subs(Omega,t,0.5))
% D+Omega
% 
