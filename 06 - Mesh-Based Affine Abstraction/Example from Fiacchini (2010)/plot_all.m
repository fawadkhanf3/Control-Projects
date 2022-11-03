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
    
    zp = subBound(3,1);
    zq = subBound(3,2);
    
    xplot = linspace(xp,xq,plot_m);
    yplot = linspace(yp,yq,plot_m);    
    zplot = linspace(zp,zq,plot_m);
    
    [X,Y,Z] = meshgrid(xplot,yplot,zplot);
    
%     pointsize = 15;  %adjust as needed
%     [X, Y, Z] = ndgrid(xplot, yplot, zplot);
    F = func(X,Y,Z);
    
    isosurface(X, Y, Z, F,-3);
    
%     subplot(1,2,1)
    scatter3(X(:),Y(:),Z(:),20,F(:));
%     subplot(1,2,2);
%     F2 = F;
%     for level = 0.2:0.2:0.8
%         isosurface(X, Y, Z, F2, level);
%     end
%     scatter3(X(:), Y(:), Z(:), 1, F(:));
    
    
    Zb = Ab(1)*X+Ab(2)*Y+Ab(3)*Z+hb;
    Zu = Au(1)*X+Au(2)*Y+Au(3)*Z+hu;
%     hold on
%     isosurface(X, Y, Z, Zb,-3);
    
%     isosurface(X, Y, Z, Zu,-3);
%
%     surf(X,Y,Z,F)
%     mesh(X,Y,F,'edgecolor','k','linewidth',1.5);
%     mesh(X,Y,Zb-eps,'edgecolor','r','linewidth',1.5);
%     mesh(X,Y,Zu+eps,'edgecolor','b','linewidth',1.5);
    
end
hold off; view(3);
end