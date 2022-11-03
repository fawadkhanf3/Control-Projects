function plot_all(final,func,n)

plot_m = 10;

figure(n);
hold all; box on;

for i = 1:size(final,1)
    Au = final{i}{1};
    hu = real(final{i}{2});
    Ab = final{i}{3};
    hb = real(final{i}{4});
    subBound = final{i}{5};
    xp = subBound(1,1);
    xq = subBound(1,2);
    
    xplot = linspace(xp,xq,plot_m);
    X = xplot;%meshgrid(xplot);
    Y = func(X);
    
    Zb = Ab(1)*X+hb;
    Zu = Au(1)*X+hu;
    
    h1 = plot(X,Y,'color','k','linewidth',1.5);
    h2 = plot(X,Zb-eps,'color','r','linewidth',1.5);
    h3 = plot(X,Zu+eps,'color','b','linewidth',1.5);
    
%     
%     mesh(X,Y,'edgecolor','k','linewidth',1.5);
%     mesh(X,Zb-eps,'edgecolor','r','linewidth',1.5);
%     mesh(X,Zu+eps,'edgecolor','b','linewidth',1.5);
    
end
title('Piecewise Affine Abstraction','FontSize',18);
xlabel('x','FontSize',18);ylabel('f(x) = sin(x)','FontSize',18);
legend([h1 h2 h3],'Nonlinear Function','Approximated Upper Hyperplane',...
    'Approximated Lower Hyperplane','FontSize',12);
hold off;grid on
end