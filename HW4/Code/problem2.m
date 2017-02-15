clear all ; close all; clc; 

syms sigma1 sigma2 A1 A2 xi real

eval = ( A1*(1-xi) + A2*(xi)) * (sigma1*(1-xi)+sigma2*(xi))
eval = simplify(int(eval,xi,0,1))

N = [1-xi;xi]*( A1*(1-xi) + A2*(xi))
f = int(N,xi,0,1)