clear all; close all; clc;

syms x real

N1 = 0.5*(2*x-1);
N2 = -2*x;
N3 = 0.5*(2*x + 1);

N = [N1 N2 N3];
Shape = 2*N'*N;

x1 = -1*(3)^(-0.5);
x2 = 1*(3)^(-0.5);


K1 = double(subs(Shape,x,x1))
K2 = double(subs(Shape,x,x2))

ans1 = 2*K1
ans2 = 2*K2