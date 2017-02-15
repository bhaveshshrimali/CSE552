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
plot(X1,Y1,X2,Y2,X3,Y3,'Linewidth',2);
grid on;
xlabel('x','FontWeight','bold','FontSize',18);
ylabel('y','FontWeight','bold','FontSize',18);
legend({'Element at t = 0','Element at t = 1','Element at t = 0.5'},'FontWeight','bold','FontSize',18);
title('Element snapshots at various time instants','FontWeight','bold','FontSize',18);
ylim([-0.05 1.5]);