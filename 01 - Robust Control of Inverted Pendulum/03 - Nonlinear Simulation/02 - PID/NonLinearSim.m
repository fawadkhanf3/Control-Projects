%% Nonlinear Simulation for Inverted Pendulum

clear;close all;
format long g;
warning off;

%%

load ControllerPID

%% Step Response

x0  = [0;0;0;0];
tic;sim('InvPendSimPID');toc;

%%

figure(1);hold on;grid on;
plot(position,'.-');shg
xlabel('t [sec]');ylabel('x [m]');
title('Position vs Time (Step Response)')

figure(2);hold on;grid on;
plot(angle,'.-');shg
xlabel('t [sec]');ylabel('\theta [deg]');
title('Angle vs Time (Step Response)')
