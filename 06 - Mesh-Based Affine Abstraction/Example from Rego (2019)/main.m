%% Abstraction (Adaptive Cruise Control Example - Nilsson 2016 Paper)

clc; clear; close all

des_error = 2;%0.2;%0.8; % precision

% Bounds
x1_min = -2.5;
x1_max = 0.5;

x2_min = -1.5;
x2_max = 3.5;

bound = [x1_min x1_max; x2_min x2_max]; % entire domain

m = 10; % number of grids along one dimension per subdomain

%%  Different continuity assumptions based on Proposition 1
%1 -> C0
%2 -> Lipchitz
%3 -> C1
%4 -> C2

algo = 4;

%% Non-linear Equations
example = 1;

if example == 1
    func = @(x1,x2) 3.*x1 - x1.^2/7 - 4.*x1.*x2./(4+x1);
elseif example == 2
    func = @(x1,x2) -2.*x2 + 3.*x1.*x2./(4+x1);
end

%% Constants for Proposition 1 (can be easily substituted with other algorithms)
syms x1 x2;

deriv{1,1} = matlabFunction(diff(func,x1) ,'Vars',[x1 x2]);
deriv{1,2} = matlabFunction(diff(func,x2),'Vars',[x1 x2]);

hess{1,1} = matlabFunction(diff(deriv{1,1},x1) ,'Vars',[x1 x2]);
hess{1,2} = matlabFunction(diff(deriv{1,1},x2),'Vars',[x1 x2]);
hess{2,1} = matlabFunction(diff(deriv{1,2},x1) ,'Vars',[x1 x2]);
hess{2,2} = matlabFunction(diff(deriv{1,2},x2),'Vars',[x1 x2]);

lambda  = value(lipchitz(deriv,bound));
lambda4 = case4(deriv,bound);
lambda5 = case5(hess,bound);

%% Main Code
tic

param.des_error = des_error;
param.func      = func;
param.m         = m;

param.algo      = algo;
param.lambda    = lambda;
param.lambda4   = lambda4;
param.lambda5   = lambda5;

final     = epsCover(bound,param);

toc

%%
num_grids = size(final,1);
figure_number = 1;
plot_all(final,func,figure_number);

display(num_grids) % Number of subregions