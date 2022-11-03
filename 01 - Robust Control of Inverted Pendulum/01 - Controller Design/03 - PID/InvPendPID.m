%% Hinf Controller Inverted Pendulum (Loop Shaping)

clear;close all;
warning('off');
set(groot,'defaultlineLineWidth',2);

s = tf('s');

%% State Space Inverted Pendulum

sys = InvPendulum(); % State Space
G   = tf(sys);       % Transfer Functions

position = G(1,:);   % Transfer Function for Position
theta    = G(2,:);   % Transfer Function for Angle

%% Controller Design (Angle)

Kth = pidtune(theta,'PID',35);

%% Closed Loop Performance

L = theta*Kth;      % Nominal Loop Transfer Function
S = inv(1+L); % Nominal Sensitivity Function
T = 1-S;      % Nominal Complementary Sensitivity Function

[GM,PM,W180,Wc] = margin(L);
MS = norm(S,'inf');
MT = norm(T,'inf');

fprintf('\n%6s%6s%6s%6s%6s\n','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         MS,MT,GM,PM,Wc);

figure(2);hold on;
step(T);title('Step Response');grid;shg

figure(3);hold on;
impulse(T);title('Impulse Response');grid;shg

%% Controller Design (Position)

aug_sys  = feedback(Kth*sys,1,1,2);
position = aug_sys(1,:);

Kx = pidtune(position,'PID');

%% Closed Loop Performance

L = position*Kx;      % Nominal Loop Transfer Function
S = inv(1+L); % Nominal Sensitivity Function
T = 1-S;      % Nominal Complementary Sensitivity Function

[GM,PM,W180,Wc] = margin(L);
MS = norm(S,'inf');
MT = norm(T,'inf');

fprintf('\n%6s%6s%6s%6s%6s\n','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         MS,MT,GM,PM,Wc);

figure(5);hold on;
step(T);title('Step Response');grid;shg

figure(6);hold on;
impulse(T);title('Impulse Response');grid;shg

%%
save ControllerPID.mat Kx Kth

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
