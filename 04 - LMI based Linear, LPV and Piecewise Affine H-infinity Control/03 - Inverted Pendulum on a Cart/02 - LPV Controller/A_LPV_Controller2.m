%% LPV Control of an Inverted Pendulum with a Cart

clearvars;clc;close all;
format short g
set(0,'DefaultLineLineWidth',2);

load LPVModels

A_LPV = A;
B_LPV = B;

%% LMI Synthesis

bounds{1,1} = [-1.3963,-1.0472];
bounds{2,1} = [-1.0472,1.0472];
bounds{3,1} = [1.0472,1.3963];

A = cell(3,1);
B = cell(3,1);

for i = 1:3
    
    A{i,1}  = A_LPV(:,:,i)
           
    b{i,1}  = zeros(4,1);
    B1{i,1} = zeros(4,1);
    
    B2{i,1} = B_LPV(:,:,i);
    
    C1{i,1} = [0.2*eye(2), zeros(2);zeros(1,4)];
    D1{i,1} = zeros(3,1);
    D2{i,1} = [zeros(2,1);0.6];
    
    lb = bounds{i,1}(:,1);
    ub = bounds{i,1}(:,2);
    
    E{i,1} = [0,2/(ub(1)-lb(1)),0,0];
    f{i,1} = -(ub(1)+lb(1))/(ub(1)-lb(1));
    
end

gam2 = 3.0^2;%sdpvar(1);
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
    F11 = temp + temp.' + mu{i}*b{i}*b{i}.' + B1{i}*B1{i}.';
    F12 = mu{i}*b{i}*f{i}.' + Q*E{i}.';
    F13 = B1{i}*D1{i}' + Q*C1{i}.' + Y{i}.'*D2{i}.';
    
    F21 = F12.';
    F22 = -mu{i}*(eye(nf) -(f{i}*f{i}.'));
    F23 = zeros(nf,nD1);
    
    F31 = F13.';
    F32 = F23.';
    F33 = -(gam2)*eye(nD1) + D1{i}*D1{i}.';
    
    if (i==2)
        F = [F11  F13;
            F31 F33];
        nF = n + nD1;
    else
        F = [F11  F12  F13;
            F21 F22  F23;
            F31 F32 F33];
        nF = n + nf + nD1;
    end
    
    LMIs = [LMIs;
        F <= -eps*eye(nF);
        mu{i} <= -eps]; %#ok<*AGROW>
    
end

% LMIs = [LMIs;gam2>=eps];

options = sdpsettings('solver','sedumi','verbose',true);

optimize(LMIs,[],options);

K = [];

for i = 1:length(A)
    K = [K;value(Y{i})/value(Q)];
end
    
save LPVController2.mat K