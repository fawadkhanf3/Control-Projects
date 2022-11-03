%% Hinf Controller Inverted Pendulum (Mixed Sensitivity)

clear;close all;
warning('off');
set(groot,'defaultlineLineWidth',2);

[w,Ha,Hb,s] = bode_opts();

%% State Space Inverted Pendulum

sys = InvPendulum(); % State Space
G   = tf(sys);       % Transfer Functions

position = G(1,:);   % Transfer Function for Position
theta    = G(2,:);   % Transfer Function for Angle

%% Design Parameters (Angle)
 
Wp.M  = 1.5;
Wp.A  = 0.1;
Wp.wB = 15;

Wt.tau  = 0.04;
Wt.r0   = 2;
Wt.rinf = 2;

Wu = [];

[Wp,Wt] = perf_weights(Wp,Wt);

%% Controller Design (Angle)

[Kth,CL,GAM,INFO] = mixsyn(theta,Wp,Wu,Wt,'Display','On');

%% Closed Loop Performance

L = theta*Kth;      % Nominal Loop Transfer Function
S = inv(1+L); % Nominal Sensitivity Function
T = 1-S;      % Nominal Complementary Sensitivity Function

[GM,PM,W180,Wc] = margin(L);
MS = norm(S,'inf');
MT = norm(T,'inf');

fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);
  
figure(1);hold on;
subplot(211),bodemag(S,'b',1/Wp,'r',Ha);legend('S','1/Wp');grid;xlim([1e-1 1e3]);
subplot(212),bodemag(Wp*S,Hb);title('Wp*S');grid;xlim([1e-1 1e3]);

figure(2);hold on;
subplot(211);bodemag(T,'b',1/Wt,'r',Ha);legend('T','1/Wt');grid;xlim([1e-1 1e3]);
subplot(212);bodemag(Wt*T,Hb);title('Wt*T');grid;xlim([1e-1 1e3]);

figure(3);hold on;
step(T);title('Step Response');grid;shg

figure(4);hold on;
impulse(T);title('Impulse Response');grid;shg


%% Design Parameters (Position)
 
clear Wp Wt Wu 

Wp.M  = 1.5;
Wp.A  = 0.1;
Wp.wB = 5;

Wt.tau  = 0.04;
Wt.r0   = 2;
Wt.rinf = 2;

Wu = [];

[Wp,Wt] = perf_weights(Wp,Wt);

%% Controller Design (Position)

aug_sys = feedback(Kth*sys,1,1,2);
position = aug_sys(1,:);

[Kx,CL,GAM,INFO] = mixsyn(position,Wp,Wu,Wt,'Display','On');

%% Closed Loop Performance

L = position*Kx;  % Nominal Loop Transfer Function
S = inv(1+L); % Nominal Sensitivity Function
T = 1-S;      % Nominal Complementary Sensitivity Function

[GM,PM,W180,Wc] = margin(L);
MS = norm(S,'inf');
MT = norm(T,'inf');

fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);
  
figure(5);hold on;
subplot(211),bodemag(S,'b',1/Wp,'r',Ha);legend('S','1/Wp');grid;xlim([1e-1 1e3]);
subplot(212),bodemag(Wp*S,Hb);title('Wp*S');grid;xlim([1e-1 1e3]);

figure(6);hold on;
subplot(211);bodemag(T,'b',1/Wt,'r',Ha);legend('T','1/Wt');grid;xlim([1e-1 1e3]);
subplot(212);bodemag(Wt*T,Hb);title('Wt*T');grid;xlim([1e-1 1e3]);

figure(7);hold on;
step(T);title('Step Response');grid;shg

figure(8);hold on;
impulse(T);title('Impulse Response');grid;shg

%%
save ControllerMS.mat Kx Kth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w,Ha,Hb,s] = bode_opts()

w = logspace(-3,3,1e3);

Ha = bodeoptions;
Ha.MagScale = 'log';
Ha.MagUnits = 'abs';

Hb = bodeoptions;
Hb.MagScale = 'linear';
Hb.MagUnits = 'abs';

s = tf('s');

end

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

function [Wp,Wt] = perf_weights(Wp,Wt)

M  = Wp.M;
A  = Wp.A;
wB = Wp.wB;

tau  = Wt.tau;
r0   = Wt.r0;
rinf = Wt.rinf;

Wp = tf([1/M wB],[1 wB*A]);
Wt = tf([tau r0],[tau/rinf 1]);

end
