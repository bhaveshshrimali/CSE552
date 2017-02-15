% 1-D linear shape functions- coupled lienar thermoelastic problem

clear all; close all; clc;

dt  = 5;
% Material Parameters:
E = 1e2;
c = 1.0;
k=1.0;
h=0.25;
m = 0.3;
d0 = [0.1 0.25 0.1 0.25 0.1 0.25 0.1 0.25]';
d = d0;

% For Isothermal Split:

A1_is = [E/h 0 -E/h 0 0 0 0 0
    -m/(2*dt) c*h/(2*dt)+k/(2*h) m/(2*dt) -k/(2*h) 0 0 0 0
    -E/h 0 (2*E)/h 0 -E/h 0 0 0
    -m/(2*dt) -k/(2*h) 0 (c*h/(1*dt)+k/(1*h)) m/(2*dt) -k/(2*h) 0 0
    0 0 -E/h 0 2*E/h 0 -E/h 0
    0 0 -m/(2*dt) -k/(2*h) 0 (c*h/(1*dt)+k/(1*h)) m/(2*dt) -k/(2*h)
    0 0 0 0 -E/h 0 2*E/h 0
    0 0 0 0 -m/(2*dt) -k/(2*h) 0 (c*h/(1*dt)+k/(1*h))];

A0_is = [0 -m/2 0 -m/2 0 0 0 0
    -m/(2*dt) c*h/(2*dt)-k/(2*h) m/(2*dt) k/(2*h) 0 0 0 0
    0 m/2 0 0 0 -m/2 0 0
    -m/(2*dt) k/(2*h) 0 2*(c*h/(2*dt)-k/(2*h)) m/(2*dt) k/(2*h) 0 0
    0 0 0 m/2 0 0 0 -m/2
    0 0 -m/(2*dt) k/(2*h) 0 2*(c*h/(2*dt)-k/(2*h)) m/(2*dt) k/(2*h)
    0 0 0 0 0 m/2 0 0
    0 0 0 0 -m/(2*dt) k/(2*h) 0 2*(c*h/(2*dt)-k/(2*h))];

d_is(:,1) = d0;
for n = 2:6
    d_is(:,n) = A1_is\(A0_is*d);
    d=d_is(:,n);
end

% For Adiabatic Split:

A1_ad = [E/h m/4 -E/h m/4 0 0 0 0
    -m/(2*dt) c*h/(2*dt)+k/(h) m/(2*dt) -k/(h) 0 0 0 0
    -E/h -m/4 (2*E)/h 0 -E/h m/4 0 0
    -m/(2*dt) -k/(h) 0 2*(c*h/(2*dt)+k/(h)) m/(2*dt) -k/(h) 0 0
    0 0 -E/h -m/4 2*E/h 0 -E/h m/4
    0 0 -m/(2*dt) -k/(h) 0 2*(c*h/(2*dt)+k/(h)) m/(2*dt) -k/(h)
    0 0 0 0 -E/h -m/4 2*E/h 0
    0 0 0 0 -m/(2*dt) -k/(h) 0 2*(c*h/(2*dt)+k/(h))];

A0_ad = [0 -m/4 0 -m/4 0 0 0 0;
    -m/(2*dt) c*h/(2*dt) m/(2*dt) 0 0 0 0 0;
    0 m/4 0 0 0 -m/4 0 0;
    -m/(2*dt) 0 0 2*(c*h/(2*dt)) m/(2*dt) 0 0 0;
    0 0 0 m/4 0 0 0 -m/4;
    0 0 -m/(2*dt) 0 0 2*(c*h/(2*dt)) m/(2*dt) 0;
    0 0 0 0 0 m/4 0 0;
    0 0 0 0 -m/(2*dt) 0 0 2*(c*h/(2*dt))];

d=d0;
d_ad(:,1) = d0;
for n = 2:6
    d_ad(:,n) = A1_ad\(A0_ad*d);
    d=d_ad(:,n);
end
t = 1:6;

% Adiabatic-Plots: 

figure(1)
subplot(2,1,1)
plot(t,d_ad(1,:),t,d_ad(3,:),t,d_ad(5,:),t,d_ad(7,:),'LineWidth',1.5)
xlabel('Time-Step','FontWeight','bold','FontSize',12)
ylabel('Displacement','FontWeight','bold','FontSize',12)
legend({'Node-1','Node-2','Node-3','Node-4'},'LineWidth',2)
title('Adiabatic: \Delta t = 0.1','FontWeight','bold','FontSize',12)
grid on

subplot(2,1,2)
plot(t,d_ad(2,:),t,d_ad(4,:),t,d_ad(6,:),t,d_ad(8,:),'LineWidth',2)
xlabel('Time-Step','FontWeight','bold','FontSize',12')
ylabel('Temperature','FontWeight','bold','FontSize',12)
legend({'Node-1','Node-2','Node-3','Node-4'},'FontSize',12)
grid on

% Isothermal Plots: 
figure(2)
% subplot(2,1,1)
% plot(t,d_is(1,:),t,d_is(3,:),t,d_is(5,:),t,d_is(7,:),'LineWidth',1.5)
% xlabel('Time-Step','FontWeight','bold','FontSize',12)
% ylabel('Displacement','FontWeight','bold','FontSize',12)
% legend({'Node-1','Node-2','Node-3','Node-4'},'LineWidth',2)
% title('Isothermal: \Delta t = 0.1','FontWeight','bold','FontSize',12)
% grid on

subplot(2,1,1)
plot(t,d_is(2,:),t,d_is(4,:),t,d_is(6,:),t,d_is(8,:),'LineWidth',2)
xlabel('Time-Step','FontWeight','bold','FontSize',12')
ylabel('Temperature','FontWeight','bold','FontSize',12)
legend({'Node-1','Node-2','Node-3','Node-4'},'FontSize',12)
grid on

