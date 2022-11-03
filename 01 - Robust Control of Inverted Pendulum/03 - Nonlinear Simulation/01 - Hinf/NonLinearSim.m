%% Nonlinear Simulation for Inverted Pendulum

clear;close all;
format long g;
warning off;

%%

Controller = 'LS';%{'LS','MS'};

switch(Controller)
    case 'LS'
        load ControllerLS
    case 'MS'
        load ControllerMS
end

%% Step Response

ref  = 0; ref2 = 0;
x0  = [0;0;0;0];
tic;sim('InvPendSim');toc;

%%

figure(1);hold on;grid on;
plot(position,'.-');shg
xlabel('t [sec]');ylabel('x [m]');
title('Position vs Time (Step Response)')

figure(2);hold on;grid on;
plot(angle,'.-');shg
xlabel('t [sec]');ylabel('\theta [deg]');
title('Angle vs Time (Step Response)')

%% Impulse Response

ref  = 0; ref2 = -1;
x0  = [0;0;0;0];
tic;sim('InvPendSim');toc;

%%

figure(3);hold on;grid on;
plot(position,'.-');shg
xlabel('t [sec]');ylabel('x [m]');
title('Position vs Time (Impulse Response)')

figure(4);hold on;grid on;
plot(angle,'.-');shg
xlabel('t [sec]');ylabel('\theta [deg]');
title('Angle vs Time (Impulse Response)')

%% Disturbance Rejection

ref = -1;
x0  = [0;0;pi/6;0];
tic;sim('InvPendSim');toc;

%%

figure(5);hold on;grid on;
plot(position,'.-');shg
xlabel('t [sec]');ylabel('x [m]');
title('Position vs Time (Disturbance Rejection)')

figure(6);hold on;grid on;
plot(angle,'.-');shg
xlabel('t [sec]');ylabel('\theta [deg]');
title('Angle vs Time (Disturbance Rejection)')
