%% Linear Controller for a Scalar first order 1D example
% Example 03 from "Numerical Approximation of the H8 norm for Nonlinear
% Systems"

clearvars;clc;close all;
format short g;warning off;
set(0,'DefaultLineLineWidth',2);

SYM = @(x)(x + x');
mag = @(v) (sqrt(sum(v.^2,2)));

tsim = 1.0;

%% System

f = @(x,w,u) -5.*x.*(1+(sin(x)).^2) + u + w;

k = 10;
g = k/sqrt(2*5-1);

x0 = 2;

% Linear System

A = -5;
G = 1;
B = 1;
C = [k;0];
D = [0;1];

nx = size(A,1);
nw = size(G,2);
nu = size(B,2);
nz = size(C,1);

%% LMI

Q = sdpvar(nx,nx,'symmetric');
Y = sdpvar(nu,nx,'full');

eps = 1e-5;

Constr1 = Q >= eps*eye(nx);

Constr2 = [SYM(A*Q+B*Y),G,(Q*C.'+Y.'*D.');
           G.',-g*eye(nw),zeros(nw,nz);
           (Q*C.'+Y.'*D.').',zeros(nz,nw),-g*eye(nz)];
       
LMIs = [Constr1;Constr2 <= -eps*eye(nx+nw+nz)];       

options = sdpsettings('solver','sedumi','verbose',0);
optimize(LMIs,[],options);

K_Linear = value(Y)/value(Q);

%% Simulation

out = sim('A_Linear_Sim',tsim);

figure;hold on;grid on;box on;
plot(out.simout.time,mag(out.simout.data(:,[2,3])),'b.-');
xlabel('t [sec]');ylabel('z');
