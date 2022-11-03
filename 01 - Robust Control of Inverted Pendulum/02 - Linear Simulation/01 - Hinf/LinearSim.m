%% Linear Simulation for Inverted Pendulum

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

%% State Space Inverted Pendulum

sys = InvPendulum(); % State Space

aug_sys = feedback(Kth*sys,1,1,2);
aug_sys = feedback(Kx*aug_sys,1,1,1);

G   = tf(aug_sys);   % Transfer Functions

position = G(1,:);   % Transfer Function for Position
theta    = G(2,:);   % Transfer Function for Angle

%%

t = 0:0.01:3;
u = t<0.1;
y = lsim(aug_sys,u,t);

%%

figure(1);hold on;
impulse(position,t,'b.-');grid on;shg
xlabel('t');ylabel('x [m]');
title('Position vs Time (Impulse Response)')

figure(2);hold on
impulse(theta,t,'b.-');grid on;shg
xlabel('t');ylabel('\theta [deg]');
title('Angle vs Time (Impulse Response)')

%%

figure(3);hold on;
step(position,t,'b.-');grid on;shg
xlabel('t');ylabel('x [m]');
title('Position vs Time (Step Response)')

figure(4);hold on
step(theta,t,'b.-');grid on;shg
xlabel('t');ylabel('\theta [deg]');
title('Angle vs Time (Step Response)')

%% Impulse

figure(5);hold on;
plot(t,y(:,1),'b.-');

figure(6);hold on
plot(t,y(:,2)*180/pi,'b.-');

%%
ref = 0;
tic;sim('LinearInvPendSim');toc;

figure(5);hold on;
plot(position,'r--');grid on;shg;
xlabel('t');ylabel('x [m]');
title('Position vs Time (Impulse Response)')
legend('MATLAB','SIMULINK');

figure(6);hold on
plot(theta,'r--');grid on;shg;
xlabel('t');ylabel('\theta [deg]');
title('Angle vs Time (Impulse Response)')
legend('MATLAB','SIMULINK');

%%

u = t>0.1;
y = lsim(aug_sys,u,t);

%% Step

figure(7);hold on;
plot(t,y(:,1),'b.-');

figure(8);hold on
plot(t,y(:,2)*180/pi,'b.-');

%%
ref = -1;
tic;sim('LinearInvPendSim');toc;

figure(7);hold on;
plot(position,'r--');grid on;shg;
xlabel('t');ylabel('x [m]');
title('Position vs Time (Step Response)')
legend('MATLAB','SIMULINK');

figure(8);hold on
plot(theta,'r--');grid on;shg;
xlabel('t');ylabel('\theta [deg]');
title('Angle vs Time (Step Response)')
legend('MATLAB','SIMULINK');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sys = InvPendulum()

g  = 9.80665;
Mc = 1.0;
Mr = 0.25;
l  = 0.5;
b  = 0.05;

A = [0 0 1 0;
    0 0 0 1;
    0 -(g*Mr/Mc) -b/Mc 0;
    0 (2*g*(Mc+Mr))/(Mc*l) 2*b/(Mc*l) 0];

B = [0;
    0;
    1/Mc;
    -2/(Mc*l)];

C = [1 0 0 0;
    0 1 0 0];

D = [0;0];

sys = ss(A,B,C,D);
sys = minreal(sys);

end