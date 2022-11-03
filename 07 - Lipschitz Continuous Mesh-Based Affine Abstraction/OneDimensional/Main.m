clearvars;clc;close all;
format short g;warning off;
set(0,'DefaultLineLineWidth',2);

%%

func  = @(x) -5.*x.*(1+sin(x).^2);
bound = [-pi,pi];

syms x;
deriv{1,1} = matlabFunction(diff(func,x),'Vars',x);
hess{1,1}  = matlabFunction(diff(deriv{1,1},x),'Vars',x);

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

param.des_error   = 0.1;
param.num_regions = 24;

flag = 0; %{0,1}:{desired n-regions, desired precision}

err = 1e3;
if (flag == 0)
    final = lipschitz_continuous_abstraction(param);
else
    param.num_regions = 4;
    while (err > param.des_error)
        final = lipschitz_continuous_abstraction(param);
        param.num_regions = param.num_regions*2;
        err = max(cell2mat(final(:,5)));
    end
end

%%

fprintf('\n     Average Error = %10.4f',mean(cell2mat(final(:,5))));
fprintf('\n     Maximum Error = %10.4f',max(cell2mat(final(:,5))));
fprintf('\n Number of Regions = %d\n',size(final,1));

h = plot_all(final,param);