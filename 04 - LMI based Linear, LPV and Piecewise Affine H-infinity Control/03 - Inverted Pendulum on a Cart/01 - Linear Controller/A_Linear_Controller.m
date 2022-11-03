%% Linear Control of an Inverted Pendulum with a Cart

clearvars;clc;close all;
format short g
set(0,'DefaultLineLineWidth',2);

s = tf('s');

% Linear Model for Inverted Pendulum on a Cart
% taken from the Book "Introduction to Feedback Control - Using Design
% Examples by Randal & Timothy

g  = 9.80665; % Acceleration due to gravity [m/s^2]
Mc = 1.0;     % Mass of Cart [kg]
Mr = 0.25;    % Mass of rod [kg]
l  = 1.0;     % Length of rod [m]
b  = 0.05;    % Viscous Force Coefficient

% State Space Representation
% \dot{x} = Ax + Bu
% y = Cx + Du

A(:,:,1) = [0 0 1 0;
    0 0 0 1;
    0 (-g*Mr/Mc) -b/Mc 0;
    0 2*g*(Mc+Mr)/(Mc*l) 2*b/(Mc*l) 0];

B1 = zeros(4,1);
B2 = [0;0;1/Mc;-2/(Mc*l)];
C1 = [0.2*eye(2), zeros(2);zeros(1,4)];
D1 = zeros(3,1);
D2 = [zeros(2,1);0.6];

np = size(A,3);
nx = size(A,1);
nw = size(B1,2);
nu = size(B2,2);
nz = size(C1,1);

for k=1:np
    P(:,:,k) = ss(A(:,:,k),[B1 B2],[C1;eye(nx)],[D1,D2;zeros(nx,nw+nu)]);
end
P.InputGroup.w  = 1:nw;
P.InputGroup.u  = (nw+1):(nw+nu);
P.OutputGroup.z = 1:nz;
P.OutputGroup.x = (nz+1):(nz+nx);

%%
X = sdpvar(nx,nx);

for i=1:np
    Z{i} = sdpvar(nu,nx);
end

g = sdpvar(1,1);
eps = 1;

LMI{1} = X >= eps*eye(nx);
SYM = @(x)(x + x');

for k=1:np
temp = C1*X + D2*Z{k};
LMI{k+1} = [SYM(A(:,:,k)*X + B2*Z{k}),     B1,     temp'   ;...
                        B1',          -g*eye(nw),   D1'    ;...
                       temp,              D1,   -g*eye(nz)] <= -eps*eye(nx+nw+nz);       
end
LMIs = [LMI{:}];

options = sdpsettings('solver','sedumi','verbose',0);
optimize(LMIs,g,options)
g = value(g);

K = value(Z{i})/value(X);

save LinearController.mat K