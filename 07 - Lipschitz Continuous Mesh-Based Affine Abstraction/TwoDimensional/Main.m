clearvars;clc;close all;
format short g;warning off;
set(0,'DefaultLineLineWidth',2);

%%

example = 3;

switch example
    case 1
        func = @(x,y) x.*cos(y);
    case 2
        func = @(x,y) x.*sin(y);
    case 3
        func = @(x,y) exp(x).*y;
end
    
% bound = [-2 2;0 2*pi];

bound = [0 pi;0 5];

syms x y;  

deriv{1,1}=matlabFunction(diff(func,x),'Vars',[x y]);
deriv{1,2}=matlabFunction(diff(func,y),'Vars',[x y]);

hess{1,1} = matlabFunction(diff(deriv{1,1},x),'Vars',[x y]);
hess{1,2} = matlabFunction(diff(deriv{1,1},y),'Vars',[x y]);
hess{2,1} = matlabFunction(diff(deriv{1,2},x),'Vars',[x y]);
hess{2,2} = matlabFunction(diff(deriv{1,2},y),'Vars',[x y]);

lambda{1,1} = value(lipchitz(deriv,bound));
lambda{2,1} = lambda{1,1};
lambda{3,1} = case4(deriv,bound);
lambda{4,1} = case5(hess,bound);

%%

param.algo   = 4;
param.m      = 10;
param.lambda = lambda;
 
param.func  = func;
param.bound = bound;

param.des_error   = 0.5;
param.num_regions = 256; % m*n (if num_regions = 64 -> 8x8 grid)

if ~((mod(sqrt(param.num_regions),1)==0) == 1)
    error('Number of Regions is not a Perfect Square');
end
    
flag = 0; %{0,1}:{desired n-regions, desired precision}

err = 1e3;
if (flag == 0)
    final = lipschitz_continuous_abstraction(param);
else
    nm = 2;
    param.num_regions = nm^2;
    while (err > param.des_error)
        final = lipschitz_continuous_abstraction(param);
        nm = nm+1; param.num_regions = nm^2;
        err = max(cell2mat(final(:,5)));
        disp(size(final));
        disp(err);
    end
end

%%

fprintf('\n     Average Error = %10.4f',mean(cell2mat(final(:,5))));
fprintf('\n     Maximum Error = %10.4f',max(cell2mat(final(:,5))));
fprintf('\n Number of Regions = %d\n',size(final,1));

h = plot_all(final,param);