%% LPV Control of an Inverted Pendulum with a Stationary Pivot Point

clearvars;clc;close all;
format short g
set(0,'DefaultLineLineWidth',2);

s = tf('s');

% Dynamics (from "On Model Based Synthesis of Embedded Control Software")

g = 9.81;  % gravitational acceleration
l = g;     % length of pendulum
m = 1/l^2; % mass of pendulum

pi_r  = -pi:0.01:pi;
b_pif = @(theta) (g/l).*sin(theta)./theta;
b_pi  = b_pif(pi_r);

b_min = min(b_pi);
b_max = max(b_pi);

A(:,:,1) = [0,1;
            b_min 0];
        
A(:,:,2) = [0,1;
            b_max 0];        
 
B1 = [0;0];
B2 = [0;1/(m*l^2)];
C1 = [0.2,0;0,0.2;0,0];
D1 = [0;0;0];
D2 = [0;0;0.6];

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
g = value(g)

K = [];
for i=1:np
    K = [K;value(Z{i})/value(X)];
end

b = [b_min;b_max];

save LPV_Controller.mat b K

%%
for k=1:np
Gcl(:,:,k) = lft(P(:,:,k),K(:,k));
end

fprintf('GAMMA: %g\t MATLAB: [%g,%g]\n',g,norm(Gcl,inf));
fprintf('Max. real(eig(Gcl)): %g\n', max(max(real(eig(Gcl)))));