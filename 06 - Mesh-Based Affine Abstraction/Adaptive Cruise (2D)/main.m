%% Abstraction (Adaptive Cruise Control Example - Nilsson 2016 Paper)

clc; clear; close all

des_error = 0.005;%0.2;%0.02; % precision

% Constants
mass = 1370;
g0 = 9.82;

f0    = 3.8e-3*mass*g0;
f1    = 2.6e-5*mass*g0;
f2    = 0.4161;

% Bounds
v_min = 0;
v_max = 35;

Fw_min = -0.3*mass*g0;
Fw_max = +0.2*mass*g0;

bound = [v_min v_max; Fw_min Fw_max]; % entire domain

m = 50; % number of grids along one dimension per subdomain

%%  Different continuity assumptions based on Proposition 1
    %1 -> C0
    %2 -> Lipchitz
    %3 -> C1
    %4 -> C2
    
algo = 4; 

%% Non-linear Equation (vdot)

func = @(v,Fw) (1/mass).*(Fw-f0-f1.*v-f2.*v.^2);
    
%% Constants for Proposition 1 (can be easily substituted with other algorithms)
syms v Fw;  

deriv{1,1} = matlabFunction(diff(func,v) ,'Vars',[v Fw]);
deriv{1,2} = matlabFunction(diff(func,Fw),'Vars',[v Fw]);

hess{1,1} = matlabFunction(diff(deriv{1,1},v) ,'Vars',[v Fw]);
hess{1,2} = matlabFunction(diff(deriv{1,1},Fw),'Vars',[v Fw]);
hess{2,1} = matlabFunction(diff(deriv{1,2},v) ,'Vars',[v Fw]);
hess{2,2} = matlabFunction(diff(deriv{1,2},Fw),'Vars',[v Fw]);

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

hyperplanes     = epsCover(bound,param);

toc

%% 
num_grids = size(hyperplanes,1);
figure_number = 1;
plot_all(hyperplanes,func,figure_number);

display(num_grids) % Number of subregions