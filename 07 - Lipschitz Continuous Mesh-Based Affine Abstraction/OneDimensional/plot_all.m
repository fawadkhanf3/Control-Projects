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
    
    X = linspace(xp,xq,plot_m);
    Y = func(X);
    
    Zb = Ab(1)*X+hb;
    Zu = Au(1)*X+hu;
    
    h1 = plot(X,Y,'color','k','linewidth',1.5);
    h2 = plot(X,Zb-eps,'color','r','linewidth',1.5);
    h3 = plot(X,Zu+eps,'color','b','linewidth',1.5);
        
end

title('Piecewise Affine Abstraction','FontSize',16);
xlabel('x','FontSize',12,'FontWeight','Bold');
ylabel('f(x)','FontSize',12,'FontWeight','Bold');
legend([h1 h2 h3],{'Nonlinear Function','Upper Hyperplane',...
    'Lower Hyperplane'},'FontSize',10','Location','Best');

end