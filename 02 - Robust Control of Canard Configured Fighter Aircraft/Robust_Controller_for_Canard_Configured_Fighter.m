%% Robust Control of a Canard Configured Fighter Aircraft

clear;close all;clc
set(0, 'DefaultLineLineWidth', 2);
format long g

% 
s = tf('s');
w = logspace(-3,3,1e3);

Ha = bodeoptions;
Ha.MagScale = 'log';
Ha.MagUnits = 'abs';

Hb = bodeoptions;
Hb.MagScale = 'linear';
Hb.MagUnits = 'abs';

plot_LS  = 0;
plot_MS  = 0;
plot_LQG = 0;
plot_PI  = 0;

% Aircraft Linear Model
sys = aircraft_model();

% Inner Loop (Stability Augmentation)
K_alpha = 1.72;
K_q = 1.6;

sys_f1 = feedback(sys,K_q,1,3);
sys_f2 = feedback(sys_f1,K_alpha,1,2);

G     = tf(sys_f2);
theta = G(4,:);

%% Loop Shaping
Gd             = theta*1.5/(s+0.35);
[K_LS,~,GAM,~] = loopsyn(theta,Gd);

% Closed Loop Performance

L_LS = theta*K_LS;  % Nominal Loop Transfer Function
S_LS = 1/(1+L_LS);  % Nominal Sensitivity Function
T_LS = 1-S_LS;      % Nominal Complementary Sensitivity Function

[GM,PM,~,Wc] = margin(L_LS);
MS = norm(S_LS,'inf');
MT = norm(T_LS,'inf');

disp('Loop Shaping:');
fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);

if (plot_LS == 1)
figure(1);hold on;grid on;box on
sigma(L_LS,'r',Gd*GAM,'k-.',Gd/GAM,'k-.',{.1,100})
legend('L','Gd*GAM','Gd/GAM')
xlabel('t [sec]');ylabel('\theta [deg]');
figure(2);hold on;grid on;box on
step(T_LS,20);title('Step Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
figure(3);hold on;grid on;box on
impulse(T_LS,20);title('Impulse Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
end
%% Mixed Sensitivity

wp.M  = 2;
wp.A  = 1e-4;
wp.wB = 0.6;

Wp = tf([1/wp.M wp.wB],[1 wp.wB*wp.A]);
Wu = tf(0.1);
Wt = [];

[K_MS,~,GAM,~] = mixsyn(theta,Wp,Wu,Wt,'Display','Off');

L_MS = theta*K_MS;  % Nominal Loop Transfer Function
S_MS = 1/(1+L_MS);  % Nominal Sensitivity Function
T_MS = 1-S_MS;      % Nominal Complementary Sensitivity Function

[GM,PM,~,Wc] = margin(L_MS);
MS = norm(S_MS,'inf');
MT = norm(T_MS,'inf');

disp('Mixed Sensitivity:');
fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);

if (plot_MS == 1)     
figure(4);hold on;grid on;box on
subplot(211),bodemag(S_MS,'b',1/Wp,'r',Ha);legend('S','1/Wp');grid;xlim([1e-1 1e3]);
subplot(212),bodemag(Wp*S_MS,Hb);title('Wp*S');xlim([1e-1 1e3]);
xlabel('t [sec]');ylabel('\theta [deg]');
figure(5);hold on;grid on;box on
step(T_MS,20);title('Step Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
figure(6);hold on;grid on;box on
impulse(T_MS,20);title('Impulse Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
end
%% LQG

Q  = 500*eye(4);
Qn = 0.9*eye(4);
R  = 10;
Rn = 0.90;
QI = 10;

QXW = 0.9*eye(5);
QXU = blkdiag(Q,R);
QWV = blkdiag(Qn,Rn);

K_LQG = lqg(sys_f2(4,:),QXU,QWV,QI,'1dof');

L_LQG = K_LQG*theta;  % Nominal Loop Transfer Function
S_LQG = 1/(1+L_LQG);  % Nominal Sensitivity Function
T_LQG = 1-S_LQG;      % Nominal Complementary Sensitivity Function

% T_LQG = feedback(K_LQG*sys_f2,1,1,4);
% T_LQG = T_LQG(4,:);

[GM,PM,~,Wc] = margin(L_LQG);
MS = norm(S_LQG,'inf');
MT = norm(T_LQG,'inf');

disp('LQG:');
fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);

if (plot_LQG == 1)
figure(7);hold on;grid on;box on;
step(T_LQG,20);title('Step Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
figure(8);hold on;grid on;box on;
impulse(T_LQG,20);title('Impulse Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
end

%%
Kp = 0.593346649788978;
Ki = 0.00362527940530768;
K_PI = Kp + Ki/s;

L_PI = theta*K_PI;  % Nominal Loop Transfer Function
S_PI = 1/(1+L_PI);  % Nominal Sensitivity Function
T_PI = 1-S_PI;      % Nominal Complementary Sensitivity Function

[GM,PM,~,Wc] = margin(L_PI);
MS = norm(S_MS,'inf');
MT = norm(T_MS,'inf');

disp('Proportional-Integral:');
fprintf('\n%6s%6s%6s%6s%6s%6s\n','gamma','Ms','Mt','GM','PM','Wc');
fprintf('%6.2f%6.2f%6.2f%6.2f%6.2f%6.2f\n',...
         GAM,MS,MT,GM,PM,Wc);

if (plot_PI == 1)
figure(9);hold on;grid on;box on;
step(T_PI,20);title('Step Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
figure(10);hold on;grid on;box on;
impulse(T_PI,20);title('Impulse Response');shg
xlabel('t [sec]');ylabel('\theta [deg]');
end

%%
figure(11);hold on;
step(T_LS,'r.-',T_MS,'b.-',T_LQG,'k.-',T_PI,'g.-',20);grid;
xlabel('t [sec]');ylabel('\theta [deg]');
legend('Loop Shaping','Mixed Sensitivity','LQG','PID','Location','Best');shg
figure(12);hold on;
impulse(T_LS,'r.-',T_MS,'b.-',T_LQG,'k.-',T_PI,'g.-',20);grid;
xlabel('t [sec]');ylabel('\theta [deg]');
legend('Loop Shaping','Mixed Sensitivity','LQG','PID','Location','Best');shg

figure(13);hold on;
t = 0:0.1:30;
u = zeros(length(t),1);
u(t>=1) = 5;
y1 = lsim(T_LS,u,t);
y2 = lsim(T_MS,u,t);
y3 = lsim(T_LQG,u,t);
y4 = lsim(T_PI,u,t);
plot(t,y1,'r.-',t,y2,'b.-',t,y3,'k.-',t,y4,'g.-');grid
xlabel('t [sec]');ylabel('\theta [deg]');
legend('Loop Shaping','Mixed Sensitivity','LQG','PID','Location','Best');shg