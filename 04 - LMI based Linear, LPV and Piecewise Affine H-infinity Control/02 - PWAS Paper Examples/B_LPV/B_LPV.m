%% LPV Controler for a Scalar first order 1D example
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

% LPV System
b = [-pi,-pi/2,0,pi/2,pi];

A(:,:,1) = -5;
A(:,:,2) = -10;
A(:,:,3) = -5;
A(:,:,4) = -10;
A(:,:,5) = -5;

G = 1;
B = 1;
C = [k;0];
D = [0;1];

np = size(A,3);
nx = size(A,1);
nw = size(G,2);
nu = size(B,2);
nz = size(C,1);

%% LMI

g2 = g^2;

Q = sdpvar(nx,nx,'symmetric');

LMIs = Q >= eps*eye(nx);

for i = 1:np
    
    Y{i} = sdpvar(nu,nx,'full');
    
    F = [SYM(A(:,:,i)*Q+B*Y{i}),G,(Q*C.'+Y{i}.'*D.');
           G.',-g2*eye(nw),zeros(nw,nz);
           (Q*C.'+Y{i}.'*D.').',zeros(nz,nw),-eye(nz)];
    
    LMIs = [LMIs; F <= -eps*eye(nx+nw+nz)];
    
end

options = sdpsettings('solver','sedumi','verbose',0);
optimize(LMIs,[],options);

K = [];
for i = 1:np
    K = [K;value(Y{i})/value(Q)];
end

b_LPV = b(:);
K_LPV = K(:);

%% Simulation

out = sim('B_LPV_Sim',tsim);

figure;hold on;grid on;box on;
plot(out.simout.time,mag(out.simout.data(:,[2,3])),'r.-');
xlabel('t [sec]');ylabel('z');