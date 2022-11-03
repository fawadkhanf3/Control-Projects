%% Quadcopter Simulation

clear;close all;clc
format long g
set(0,'DefaultLineLineWidth',2);

%%

% Initial Conditions

x0    = [0;0;0;pi/3;0;0;0;0;0;0;0;0];
zeta0 = 0.1;
xi0   = 0.1;

%Simulation

out = sim('Quadcopter_01',[0 30]);

figure(1);hold on;grid on;box on
t = 1:30;
plot3(0.5*cos(t/2),1/2*sin(t/2),1+t/10);
plot3(out.simout.data(:,1),out.simout.data(:,2),out.simout.data(:,3));
view(3);%legend('Reference','Tracking');