%% Optimal Controller Inverted Pendulum (LQR)

clear;close all;warning('off');

s = tf('s');

%% State Space Inverted Pendulum

sys = InvPendulum(); % State Space

%% Design Parameters

R = 1;
Q = sys.C'*sys.C;

%% Controller Design

K_lqr   = lqr(sys,Q,R);
aug_A   = sys.A - sys.B*K_lqr;
aug_sys = ss(aug_A,sys.B,sys.C,sys.D);

G = tf(aug_sys);
position = G(1,:);

step(position,'.-');
return

%% Closed Loop Performance
figure(1);hold on;
step(aug_sys);title('Step Response');grid;shg
figure(2);hold on;
impulse(aug_sys);title('Impulse Response');grid;shg

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
