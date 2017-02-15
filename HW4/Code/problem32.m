clear all; close all; clc;

base = 2; height = 1;

x1 = [0 0];
x2 = @(t) [base 0.5*base*t];
x3 = @(t) [height*t height];

%Plotting the triangle at t = 0s
A2 = x2(0); A3 = x3(0);
X1 = [x1(1) A2(1) A3(1) 0];
Y1 = [x1(2) A2(2) A3(2) 0];
%Plotting the triangle at t = 1s
B2 = x2(1); B3 = x3(1);
X2 = [x1(1) B2(1) B3(1) 0];
Y2 = [x1(2) B2(2) B3(2) 0];
%Plotting the triangle at t = 0.5s
C2 = x2(0.5); C3 = x3(0.5);
X3 = [x1(1) C2(1) C3(1) 0];
Y3 = [x1(2) C2(2) C3(2) 0];
%Plotting the triangle at t = sqrt(2)s
D2 = x2((2)^0.5); D3 = x3((2)^0.5);
X4 = [x1(1) D2(1) D3(1) 0];
Y4 = [x1(2) D2(2) D3(2) 0];

plot(X1,Y1,X4,Y4,'Linewidth',2);
grid on;
xlabel('x','FontWeight','bold','FontSize',18);
ylabel('y','FontWeight','bold','FontSize',18);
h = legend({'Element at t = 0','Element at t = $\sqrt{2}$'},'FontWeight','bold','FontSize',18);
set(h,'Interpreter','latex')
title('Element snapshots at various time instants','FontWeight','bold','FontSize',18);
ylim([-0.05 1]);