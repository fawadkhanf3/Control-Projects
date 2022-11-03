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

%% Design Parameters

wc = 20;
Gd = wc/s;

%% Controller Design (Angle)

[Kth,CL,GAM,INFO] = loopsyn(theta,Gd);

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

figure(1);
sigma(L,'r',Gd*GAM,'k-.',Gd/GAM,'k-.',{.1,100})  % plot result
legend('L','Gd*GAM','Gd/GAM')

figure(2);hold on;
step(T);title('Step Response');grid;shg

figure(3);hold on;
impulse(T);title('Impulse Response');grid;shg

%% Controller Design (Position)

Gd = 5/s;
aug_sys  = feedback(Kth*sys,1,1,2);
position = aug_sys(1,:);

[Kx,CL,GAM,INFO] = loopsyn(position,Gd);

%% Closed Loop Performance

L = position*Kx;      % Nominal Loop Transfer Function
S = inv(1+L); % Nominal Sensitivity Function
T = 1-S;      % Nominal Complementary Sensitivity Function

[GM,PM,W180,Wc] = margin(L);
MS = norm(S,'inf');
MT = norm(T,'inf');

fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);

figure(4);
sigma(L,'r',Gd*GAM,'k-.',Gd/GAM,'k-.',{.1,100})  % plot result
legend('L','Gd*GAM','Gd/GAM')

figure(5);hold on;
step(T);title('Step Response');grid;shg

figure(6);hold on;
impulse(T);title('Impulse Response');grid;shg

%%
save ControllerLS.mat Kx Kth

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
