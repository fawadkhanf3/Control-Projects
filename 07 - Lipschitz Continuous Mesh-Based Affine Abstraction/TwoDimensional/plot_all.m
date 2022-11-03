function h = plot_all(final,param)

func   = param.func;
plot_m = param.m;

h = figure(1);clf;hold all;grid on;box on;

for i = 1:size(final,1)
    
    Au = final{i,1};
    hu = real(final{i,2});
    Ab = final{i,3};
    hb = real(final{i,4});
    subBound = final{i,6};
    xp = subBound(1,1);
    xq = subBound(1,2);

    yp = subBound(2,1);
    yq = subBound(2,2);
    
    xplot = linspace(xp,xq,plot_m);
    yplot = linspace(yp,yq,plot_m);    
    [X,Y] = meshgrid(xplot,yplot);
    Z = func(X,Y);
    
    Zb = Ab(1)*X+Ab(2)*Y+hb;
    Zu = Au(1)*X+Au(2)*Y+hu;
    mesh(X,Y,Z,'edgecolor','k','linewidth',1.5);
    mesh(X,Y,Zb-eps,'edgecolor','r','linewidth',1.5);
    mesh(X,Y,Zu+eps,'edgecolor','b','linewidth',1.5);
    
end
hold off; view(3);
end