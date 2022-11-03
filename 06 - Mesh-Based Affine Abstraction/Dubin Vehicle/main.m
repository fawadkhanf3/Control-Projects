% 1D Dubin Vehicle with Constant v

clc; clear; close all

des_error = 0.5; % precision

bound = [0 2*pi]; % entire domain [x_min x_max,; y_min y_max]
m = 10; % number of grids along one dimension per subdomain
v = 1;

%%  Different continuity assumptions based on Proposition 1
    %1 -> C0
    %2 -> Lipchitz
    %3 -> C1
    %4 -> C2
    
algo = 4; 

%% Different functions 
% Examples for v*cos(y) and v.*sin(y)

example = 2;

switch example
    case 1
        func = @(y) v.*cos(y);
    case 2
        func = @(y) v.*sin(y);
end
    
%% Constants for Proposition 1 (can be easily substituted with other algorithms)
syms y;  

deriv{1,1} = matlabFunction(diff(func,y),'Vars',y);
hess{1,1} = matlabFunction(diff(deriv{1,1},y),'Vars',y);

lambda = value(lipchitz(deriv,bound));
lambda4 = case4(deriv,bound);
lambda5 = case5(hess,bound);

%% Main Code
tic

param = cell(7,1);

param{1} = des_error;
param{2} = lambda;
param{3} = func;
param{4} = m;
param{5} = algo;
param{6} = lambda4;
param{7} = lambda5;

final = recursive_error(bound,param);

toc

%% 
num_grids = size(final,1);
figure_number = example;
plot_all(final,func,figure_number);

display(num_grids) % Number of subregions