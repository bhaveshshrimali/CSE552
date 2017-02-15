clear all; close all; clc; 

xa = 0;
xb = 1;
ya = 0;
yb = 1;
xinc = (xb-xa)/4;
yinc = (yb-ya)/4;

xvector = [xa xa+xinc xa+2*xinc xa+3*xinc xb xb xb xb xb xb-xinc xb-2*xinc xb-3*xinc xa xa xa xa xa];
yvector = [ya ya ya ya ya ya+yinc  ya+2*yinc ya+3*yinc yb yb yb yb yb yb-yinc yb-2*yinc yb-3*yinc ya];
plot(xvector,yvector,'LineWidth',2);
hx = xlabel('$\bf x$','FontWeight','bold','FontSize',18);
set(hx,'Interpreter','latex');
hy = ylabel('$\bf y$','FontWeight','bold','FontSize',18);
set(hy,'Interpreter','latex');
leg = legend({'\bf 1x1 Mesh element'},'FontWeight','bold','FontSize',18,'Location','northwest');
set(leg,'Interpreter','latex');
t = title('\bf Problem 3(a)','FontWeight','bold','FontSize',18);
xlim([-0.5 1.5]);
ylim([-0.5 1.5]);
set(t,'Interpreter','latex');
grid on;
set(gca,'xtick',[-0.5:0.25:1.5]);
set(gca,'ytick',[-0.5:0.25:1.5]);
