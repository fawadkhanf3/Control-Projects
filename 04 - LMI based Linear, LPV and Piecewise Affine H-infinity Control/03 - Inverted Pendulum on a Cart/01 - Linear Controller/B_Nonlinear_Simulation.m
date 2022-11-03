%% Test Linear Controller for Inverted Pendulum on Cart

clc;clear;close all
format short g
set(0,'DefaultLineLineWidth',2);

load LinearController K

p.g0 = 9.80665; % Acceleration due to gravity [m/s^2]
p.Mc = 1.0;     % Mass of Cart [kg]
p.Mr = 0.25;    % Mass of rod [kg]
p.l  = 0.5;     % Length of rod [m]
p.b  = 0.05;    % Viscous Force Coefficient

p.K = K;

x0 = [0;60*pi/180;0;0];
ts = [0 10];

[t,y] = ode45(@(t,x) dots(t,x,p),ts,x0);

figure(1);
subplot(211);hold on;grid on;box on;
plot(t,rad2deg(y(:,2)),'b.-');shg;
xlabel('t','FontSize',18,'FontWeight','Bold');
ylabel('\theta','FontSize',18,'FontWeight','Bold');
subplot(212);hold on;grid on;box on;
plot(t,y(:,1),'b.-');shg;
xlabel('t','FontSize',18,'FontWeight','Bold');
ylabel('z','FontSize',18,'FontWeight','Bold');

%%
function dx = dots(~,x,p)

u = p.K*x;

% yy     = x(1);
theta  = x(2);
dy     = x(3);
dtheta = x(4);

M = [p.Mr+p.Mc               p.Mr*p.l/2*cos(theta);
     p.Mr*p.l/2*cos(theta)   p.Mr*p.l^2/4];
 
c = [p.Mr*p.l/2*dtheta^2*sin(theta) + u - p.b*dy;
      p.Mr*p.l/2*p.g0*sin(theta)];
    
tmp = M\c;

dx = [dy;dtheta;tmp(:)];

end