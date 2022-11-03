%% Linear Simulation for Inverted Pendulum

clear;close all;
format long g;
warning off;

%% State Space Inverted Pendulum

sys = InvPendulum(); % State Space

sim('InvPendSimPID');

%%
figure(1);hold on;grid on;
plot(position,'.-');shg;

figure(2);hold on;grid on;
plot(theta,'.-');shg;

%%
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
