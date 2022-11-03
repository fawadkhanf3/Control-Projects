% Linear program for Affine Abstraction in Theorem 1
function [error,planes] = abstract(bound,param)

func  = param.func;

m    = param.m;
algo = param.algo;

lambda  = param.lambda;
lambda4 = param.lambda4;
lambda5 = param.lambda5;

xp = bound(1,1);
xq = bound(1,2);

yp = bound(2,1);
yq = bound(2,2);

%% Affine abstraction
% Overapproximation (Eqs. 16(a) and 16(b))

x = sdpvar(7,1);

Xo = sample(m,bound); % Creates a Mesh for v and Fw within their bounds
Yo = zeros(size(Xo,1),1);

for i = 1:size(Xo,1)
        Yo(i) = func(Xo(i,1),Xo(i,2)); % Computes vdot against sample mesh
end

A = zeros(2*m*m+5,7);

A(1:m*m,1:3) = -Xo;
A(m*m+1:2*m*m,4:6) = Xo;

% Distance between the two points on the corner of the domain
% is less than theta (Eq. 16(c))

A(2*m*m+1,:) = [xp yp 1 -[xp yp 1] -1];
A(2*m*m+2,:) = [xp yq 1 -[xp yq 1] -1];
A(2*m*m+3,:) = [xq yp 1 -[xq yp 1] -1];
A(2*m*m+4,:) = [xq yq 1 -[xq yq 1] -1];

A(2*m*m+5,:) = [zeros(1,6) -1]; % theta >= 0

b = [-Yo;Yo;zeros(5,1)];

%%
% Separation calculations

delta = sqrt((bound(1,2)-bound(1,1))^2 +...
             (bound(2,2)-bound(2,1))^2)/(m-1);
         
Rs = delta*sqrt(1/3);

if (algo == 1)
    sep = 2*lambda*Rs;
elseif (algo == 2) 
    sep = lambda*Rs;
elseif (algo == 3)
    sep = lambda4*Rs;
elseif (algo==4)
    sep = lambda5*Rs^2/2;
else
    error('Undefined algo')
end

%%     
% solver 

% Objective function
obj = x(7); % Minimize theta

options = sdpsettings('verbose',0);
optimize(A*x<=b,obj,options);

Au = [x(1) x(2)];
Ab = [x(4) x(5)];
hu = x(3);
hb = x(6);

% Account for approximation error corresponding to Proposition 1
hu = hu+sep;
hb = hb-sep;

%% Outputs
% Account for approximation error corresponding to Proposition 1
error = value(x(7))+2*sep;

planes = cell(5,1);

planes{1} = value(Au);
planes{2} = value(hu);
planes{3} = value(Ab);
planes{4} = value(hb);
planes{5} = bound;

end

