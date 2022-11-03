%% Quadratic Stability & Control of Piecewise Linear Systems
% A. Hassibi & S. Boyd
clear; close all; %clc

% Only Mass-Spring Part
M1 = 1;   M2 = 1;
k1 = 1;   k2 = 1;
b1 = 0.1; b2 = 0.1;

M = blkdiag(M1,M2);

C = [b1+b2, -b2;
       -b2,  b2];

K = [k1+k2, -k2;
       -k2,  k2];
   
A0 = [zeros(2), eye(2); -M\K, -M\C];
B0 = [zeros(2); M\eye(2)];

clearvars M C K M1 M2 b1 b2 k1 k2

% Actuator
tau = 1/10;

%% Piecewise Linear Systems Model
% for x \in R{i}
% xd = A{i}*x + b{i} + B1{i}*w + B2{i}*u
%  z = C1{i}*x       + D1{i}*w + D2{i}*u

% aa = 10e3; 

%% R{1} -aa <= x5 <= -1

A{1}  = blkdiag(A0,-1/tau);
b{1}  = [0;0;0;-1;0];
B1{1} = [B0(:,1);0];
B2{1} = [zeros(4,1);1/tau];
C1{1} = [1 zeros(1,4); zeros(1,5)];
D1{1} = [0;0];
D2{1} = [0;0.1];
% E{1} = 2*[zeros(1,4),1]/(-1 + aa);
% f{1} = -(-1-aa)/(-1+aa);

%% R{2}   -1 <= x5 <= 1

A{2}  = [A0, B0(:,2);zeros(1,4), -1/tau];
b{2}  = zeros(5,1);
B1{2} = [B0(:,1);0];
B2{2} = [zeros(4,1);1/tau];
C1{2} = [1 zeros(1,4); zeros(1,5)];
D1{2} = [0;0];
D2{2} = [0;0.1];
% E{2} = 2*[zeros(1,4),1]/(1+1);
% f{2} = 0;

% R{3}    1 <= x5 <= aa
A{3}  = blkdiag(A0,-1/tau);
b{3}  = [0;0;0;1;0];
B1{3} = [B0(:,1);0];
B2{3} = [zeros(4,1);1/tau];
C1{3} = [1 zeros(1,4); zeros(1,5)];
D1{3} = [0;0];
D2{3} = [0;0.1];
% E{3} = 2*[zeros(1,4),1]/(-1 + aa);
% f{3} = (-1-aa)/(-1+aa);

aa = 10e1; 
% 
E{1} = 2*[zeros(1,4),1]/(-1 + aa);
f{1} = -(-1-aa)/(-1+aa);

E{2} = 2*[zeros(1,4),1]/(1+1);
f{2} = 0;

E{3} = 2*[zeros(1,4),1]/(-1 + aa);
f{3} = (-1-aa)/(-1+aa);

% E{1} = 2*blkdiag(0,0,0,0,1)/(-1 + aa);
% f{1} = [0;0;0;0;-(-1-aa)/(-1+aa)];
% 
% E{2} = 2*blkdiag(0,0,0,0,1)/(1+1);
% f{2} = [0;0;0;0;0];
% 
% E{3} = 2*blkdiag(0,0,0,0,1)/(-1 + aa);
% f{3} = [0;0;0;0;(-1-aa)/(-1+aa)];

%% LMI for Synthesis [sec 7.2]
gam2 = 7.0^2;
eta = sdpvar;
eps  = 1e-5;

M = length(A);
n = size(A{1},1);
m = size(B2{1},2);

Q = sdpvar(n,n,'symmetric');
Y = cell(M,1); 
mu = cell(M,1);

LMIs = Q >= eps*eye(n);

for i = 1:M
    Y{i} = sdpvar(m,n,'full');
    mu{i} = sdpvar(1);
    
    nf  = size(f{i},1);
    nD1 = size(D1{i},1);
    temp = A{i}*Q + B2{i}*Y{i};
    F11 = temp + temp.' + mu{i}*b{i}*b{i}.' + eta*B1{i}*B1{i}.';
    F12 = mu{i}*b{i}*f{i}.' + Q*E{i}.';
    F13 = B1{i}*D1{i}.' + Q*C1{i}.' + Y{i}.'*D2{i}.';
    
    F21 = F12.';
    F22 = -mu{i}*(eye(nf) -(f{i}*f{i}.'));
    F23 = zeros(nf,nD1);
    
    F31 = F13.';
    F32 = F23.';
    F33 = -(1)*eye(nD1) + D1{i}*D1{i}.';
    
    if(i==2) %% Origin Lies in RE{2} (so 2nd block row & column removed)
        F = [F11  F13;...
            F31 F33];
        nF = n + nD1;
    else
        F = [F11  F12  F13;...
            F21 F22  F23;...
            F31 F32 F33];
        nF = n + nf + nD1;
    end
    
    LMIs = [LMIs;
            F <= -eps*eye(nF);
            mu{i} <= -eps;
            eta >= eps;
            ]; %#ok<*AGROW>
    
end

options = sdpsettings('solver','sedumi','verbose',true);
% options = sdpsettings('solver','lmilab','verbose',true);


optimize(LMIs,-eta,options)

for i=1:length(A)
    K{i} = value(Y{i})/value(Q);
    [value(Y{i})/value(Q)]
end
