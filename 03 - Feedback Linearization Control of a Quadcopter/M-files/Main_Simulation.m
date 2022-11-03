% Fault Detection of a Quadrotor with Nonlinaer State Estimation using 
% Constrained Zonotopes

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

rng(8); % for reproducibility

%%
fp = struct;
fp = constants(fp);
fp = navigation_parameters(fp);
fp = guidance_parameters(fp);
fp = controller_parameters(fp);

% Initial Conditions
t0    = -fp.T;
x0    = [0;0;0;pi/3;0;0;0;0;0;0;0;0];
zeta0 = 0.1;
xi0   = 0.1;
tf    = 30;

% Simulation

[t,y] = simulate_quadrotor(fp,[t0,tf],x0,zeta0,xi0);

figure(1);hold on;grid on;box on
plot3(fp.xref(t),fp.yref(t),fp.zref(t));
plot3(y(:,1),y(:,2),y(:,3));
view(3);