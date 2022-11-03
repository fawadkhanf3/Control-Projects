%% Quadratic Stability & Control of Piecewise Linear Systems
% A. Hassibi & S. Boyd
clear; close all; %clc

% Only Mass-Spring Part
M = [1,0;0,1]; C = [0.1+0.1,-0.1;-0.1,0.1]; K = [1+1,-1;-1,1];
A0 = [zeros(2), eye(2); -M\K, -M\C];
B0 = [zeros(2); M\eye(2)];
clearvars M C K
% Actuator
tau = 1/10;
%% Piecewise Linear Systems Model
% for x \in R{i}
% xd = A{i}*x + b{i} + B1{i}*w + B2{i}*u
%  z = C1{i}*x + D1{i}*w + D2{i}*u

aa = 100; % Only works for aa <= 139, I dont know why? there must be some issue in code or some typo in paper

% R{1} -aa <= x5 <= -1
R{1}  = struct('h',[zeros(4,2);1 -1],'g',[-1;aa]);  % Polytopic
RE{1} = struct('E',2*[zeros(1,4),1]/(-1 + aa),'f',-(-1-aa)/(-1+aa)); % Degenrate Ellipsoid
A{1}  = [A0, zeros(4,1);zeros(1,4), -1/tau];
b{1}  = [0;0;0;-1;0];
B1{1} = [B0(:,1);0];
B2{1} = [zeros(4,1);1/tau];
C1{1} = [1 zeros(1,4); zeros(1,5)];
D1{1} = [0;0];
D2{1} = [0;0.1];

% R{2}   -1 <= x5 <= 1
R{2}  = struct('h',[zeros(4,2);1 -1],'g',[1;1]); % Polytopic
RE{2} = struct('E',2*[zeros(1,4),1]/(1 + 1),'f',0); % Degenrate Ellipsoid
A{2}  = [A0, B0(:,2);zeros(1,4), -1/tau];
b{2}  = zeros(5,1);
B1{2} = [B0(:,1);0];
B2{2} = [zeros(4,1);1/tau];
C1{2} = [1 zeros(1,4); zeros(1,5)];
D1{2} = [0;0];
D2{2} = [0;0.1];

% R{3}    1 <= x5 <= aa
R{3}  = struct('h',[zeros(4,2);1 -1],'g',[aa;-1]); % Polytopic
RE{3} = struct('E',2*[zeros(1,4),1]/(-1 + aa),'f',(-1-aa)/(-1+aa)); % Degenrate Ellipsoid
A{3}  = [A0, zeros(4,1);zeros(1,4), -1/tau];
b{3}  = [0;0;0;1;0];
B1{3} = [B0(:,1);0];
B2{3} = [zeros(4,1);1/tau];
C1{3} = [1 zeros(1,4); zeros(1,5)];
D1{3} = [0;0];
D2{3} = [0;0.1];

%% LMI for Synthesis [sec 7.2]
gam2 = 7.0^2;
eps = 1e-5;

M = length(A);
n = size(A{1},1);
m = size(B2{1},2);

Q = sdpvar(n,n);
Y = cell(M,1); mu = cell(M,1);
for i=1:M
    Y{i} = sdpvar(m,n,'full');
    mu{i} = sdpvar(1);
end

LMI1 = Q >= eps*eye(n);

LMI2 = [];
for i=1:M
    LMI2 = [LMI2, mu{i} <= -eps];
end

LMI3 = [];
for i=1:M
    nf  = size(RE{i}.f,1);
    nD1 = size(D1{i},1);
    temp = A{i}*Q + B2{i}*Y{i};
    F11 = temp + temp.' + mu{i}*b{i}*b{i}.' + B1{i}*B1{i}.';
    
    F12 = mu{i}*b{i}*RE{i}.f.' + Q*RE{i}.E.';
    F13 = B1{i}*D1{i}' + Q*C1{i}.' + Y{i}.'*D2{i}.';
    
    F22 = -mu{i}*(eye(nf) -(RE{i}.f*RE{i}.f.'));
    F23 = zeros(nf,nD1);
    F33 = -(gam2)*eye(nD1) + D1{i}*D1{i}.';
    if(i==2) %% Origin Lies in RE{2} (so 2nd block row & column removed)
        F = [F11  F13;...
            F13' F33];
        nF = n + nD1;
    else
        F = [F11  F12  F13;...
            F12' F22  F23;...
            F13' F23' F33];
        nF = n + nf + nD1;
    end
    
    LMI3 = [LMI3, F <= -eps*eye(nF)];
end

LMIs = [LMI1 LMI2 LMI3];

options = sdpsettings('solver','sedumi','verbose',0);

optimize(LMIs,[],options)

for i=1:length(A)
    K{i} = value(Y{i})*inv(value(Q));
end
