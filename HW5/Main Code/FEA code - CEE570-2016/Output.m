switch El.bf
    case 'False'
        
        figure(1)
        subplot(3,1,1)
        plot(time,dN2x,time,dN4x,'LineWidth',2);
        hy = ylabel('\bf $\Delta^x_{n_2}$ and $\Delta^x_{n_4}$(Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Node 2','\bf Node 4'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf X-Displacement Time History','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,2)
        plot(time,dN2y,time,dN4y,'LineWidth',2);
        hy = ylabel('\bf $\Delta^y_{n_2}$ and $\Delta^y_{n_4}$(Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Node 2','\bf Node 4'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Y-Displacement Time History','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,3)
        plot(time,Energ,'LineWidth',2);
        % ylim([0,1]);
        hy = ylabel('\bf Energy of the System $\mathcal{E}(u)$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Energy'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t2 = title('\bf Energy Time History','FontWeight','bold','FontSize',12);
        set(t2,'Interpreter','latex');
        ylim([0.1,0.2]);
        grid on
        
        figure(2)
        subplot(3,1,1)
        plot(time,Sigma(:,1),time,Epsilon(:,1),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{xx}$ and $\varepsilon_{xx}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{xx}(t)$','\bf $\varepsilon_{xx}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\sigma_{xx}$ \bf and $\varepsilon_{xx}(t)$','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,2)
        plot(time,Sigma(:,2),time,Epsilon(:,2),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{yy}$ and $\varepsilon_{yy}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{yy}(t)$','$\bf\varepsilon_{yy}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\bf\sigma_{yy}$ \bf and $\bf\varepsilon_{yy}$','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,3)
        plot(time,Sigma(:,3),time,Epsilon(:,3),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{xy}$ and $\varepsilon_{xy}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{xy}(t)$','$\bf\varepsilon_{xy}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\bf\sigma_{xy}$ \bf and $\bf\varepsilon_{xy}$','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
            
    case 'True'
        
        
        figure(1)
        subplot(3,1,1)
        plot(time,dN2x,time,dN4x,'LineWidth',2);
        hy = ylabel('\bf $\Delta^x_{n_2}$ and $\Delta^x_{n_4}$(Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Node 2','\bf Node 4'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf X-Displacement Time History \bf(Body Force)','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,2)
        plot(time,dN2y,time,dN4y,'LineWidth',2);
        hy = ylabel('\bf $\Delta^y_{n_2}$ and $\Delta^y_{n_4}$(Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Node 2','\bf Node 4'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Y-Displacement Time History \bf(Body Force)','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,3)
        plot(time,Energ,'LineWidth',2);
        % ylim([0,1]);
        hy = ylabel('\bf Energy of the System $\mathcal{E}(u)$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf Energy'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t2 = title('\bf Energy Time History \bf(Body Force)','FontWeight','bold','FontSize',12);
        set(t2,'Interpreter','latex');
%         ylim([0.1,0.2]);
        grid on
        
        figure(2)
        subplot(3,1,1)
        plot(time,Sigma(:,1),time,Epsilon(:,1),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{xx}$ and $\varepsilon_{xx}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{xx}(t)$','\bf $\varepsilon_{xx}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\sigma_{xx}$\bf \bf and $\varepsilon_{xx}(t)$ (Body Force)','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,2)
        plot(time,Sigma(:,2),time,Epsilon(:,2),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{yy}$ and $\varepsilon_{yy}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{yy}(t)$','$\bf\varepsilon_{yy}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\bf\sigma_{yy}$ \bf and $\bf\varepsilon_{yy}$ \bf(Body Force)','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
        
        subplot(3,1,3)
        plot(time,Sigma(:,3),time,Epsilon(:,3),'LineWidth',2);
        hy = ylabel('\bf $\sigma_{xy}$ and $\varepsilon_{xy}$ (Elem=1)','FontWeight','bold','FontSize',12);
        set(hy,'Interpreter','latex');
        hx = xlabel('\bf Time(t)','FontWeight','bold','FontSize',12);
        set(hx,'Interpreter','latex');
        leg = legend({'\bf $\sigma_{xy}(t)$','$\bf\varepsilon_{xy}(t)$'},'FontWeight','bold','FontSize',12,'Location','northeast');
        set(leg,'Interpreter','latex');
        t1 = title('\bf Time History $\bf\sigma_{xy}$ \bf and $\bf\varepsilon_{xy}$ \bf(Body Force)','FontWeight','bold','FontSize',12);
        set(t1,'Interpreter','latex');
        grid on
end