function plot_all(hyperplanes,func,n)

plot_m = 50;

figure(n);
hold all; box on;

for i = 1:size(hyperplanes,1)
    Au = hyperplanes{i}{1};
    hu = real(hyperplanes{i}{2});
    Ab = hyperplanes{i}{3};
    hb = real(hyperplanes{i}{4});
    subBound = hyperplanes{i}{5};
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
    
    DATA(i).Z = Z;
    DATA(i).Zb = Zb;
    DATA(i).Zu = Zu;

    mesh(X,Y,Z,'edgecolor','k','linewidth',1.5);
    mesh(X,Y,Zb-eps,'edgecolor','r','linewidth',1.5);
    mesh(X,Y,Zu+eps,'edgecolor','b','linewidth',1.5);
    
end
hold off; view(3);
save Zs.mat DATA
end