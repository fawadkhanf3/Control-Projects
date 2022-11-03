%% Test Piecewise Controller for Inverted Pendulum (without Cart)

clc;clear;close all
format short g
set(0,'DefaultLineLineWidth',2);

load LinearController K
p.K = K;

%
p.g = 9.81;  % gravitational acceleration
p.l = p.g;     % length of pendulum
p.m = 1/p.l^2; % mass of pendulum

x0 = [deg2rad(60);0];
ts = [0 10];

options = odeset('MaxStep',0.1);

[t,y] = ode45(@(t,x) dots(t,x,p),ts,x0,options);

figure(1);hold on;grid on;box on;
plot(t,rad2deg(y(:,1)),'k.-');shg
xlabel('t [sec]','FontSize',18,'FontWeight','Bold');
ylabel('\theta [deg]','FontSize',18,'FontWeight','Bold');

%%
function dx = dots(~,x,p)

u = p.K*x;

dx = [x(2); (p.g/p.l)*sin(x(1)) + 1/(p.m*p.l^2)*u];

end