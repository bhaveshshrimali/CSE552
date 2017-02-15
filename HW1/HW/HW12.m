%Nonlinear Finite Element Method : Fall 2016
clear all; close all; clc;
syms xe1 he xe2 xe3 xc de1 de2 de3 real;

Jacb = 0.5*xe1*(2*xc-1) - 2*xc*xe2 + 0.5*xe3*(2*xc+1);
N1 = xc*(xc-1)/2;
N2 = (1-xc)*(1+xc);
N3 = xc*(xc+1)/2;

N = [N1*N1 N1*N2 N1*N3;
     N1*N2 N2*N2 N2*N3;
     N1*N3 N2*N3 N3*N3];
 
 N = N*Jacb;
 Fext_matrix = simplify(subs(N,xc,-(3)^(-0.5)) + subs(N,xc,(3)^(-0.5))) 
 % The above is true, in general, for any given values of xe1, xe2, xe3
 % However numerical evaluation is presented corresponding to the following
 % values
 
 J_1 = 1/Jacb;
 Na = [diff(N1,xc);diff(N2,xc);diff(N3,xc)];
 Na =simplify(Na);
 Nb = [N1 N2 N3];
 
 d = Na' * [de1;de2;de3];
 d = simplify(d);
 K_matrix1 = J_1 * Na * Nb * d;
 K_matrix1 = simplify(K_matrix1);
 
 K_matrix2 = J_1 * (Na) * (Na)';
 
 xc = -(3)^(-0.5);
 K1_1 = simplify(subs(K_matrix1));
 K2_1 = simplify(subs(K_matrix2));
 
 xc = (3)^(-0.5);
 K1_2 = simplify(subs(K_matrix1));
 K2_2 = simplify(subs(K_matrix2));
 
 % Sample Calculations for some specific values of xe_i 
 xe1 = 0;
 xe2 = 1/2; 
 xe3 = 1;
 
 Fext_matrix = simplify(subs(Fext_matrix));
K1_1 = simplify(subs(K1_1))
K1_2 = simplify(subs(K1_2))
K2_1 = simplify(subs(K2_1))
K2_2 = simplify(subs(K2_2))
 
 
 
 