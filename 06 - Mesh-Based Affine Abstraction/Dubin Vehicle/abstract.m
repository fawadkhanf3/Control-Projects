% Linear program for Affine Abstraction in Theorem 1
function [err,planes] = abstract(bound,param)

lambda5 = param{7};
lambda4 = param{6};
algo = param{5};
lambda = param{2};
func = param{3};
m = param{4};

xp = bound(1,1);
xq = bound(1,2);

%% Affine abstraction
% Overapproximation (Eqs. 16(a) and 16(b))

x = sdpvar(5,1);

Xo = sample(m,bound);
Yo = zeros(size(Xo,1),1);

for i = 1:size(Xo,1)
        Yo(i) = func(Xo(i,1));
end

A = zeros(m*m*2+5,5);

A(1:m*m,1:2) = -Xo;
A(m*m+1:2*m*m,3:4) = Xo;

% Distance between the two points on the corner of the domain
% is less than theta (Eq. 16(c))

A(2*m*m+1,:) = [xp 1 -[xp 1] -1];
A(2*m*m+2,:) = [xp 1 -[xp 1] -1];
A(2*m*m+3,:) = [xq 1 -[xq 1] -1];
A(2*m*m+4,:) = [xq 1 -[xq 1] -1];

A(2*m*m+5,:) = [zeros(1,4) -1]; % theta >= 0

b = [-Yo;Yo;zeros(5,1)];

%%
% Separation calculations

delta = (bound(1,2)-bound(1,1))/(m-1);
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
obj = x(5); % Minimize theta

options = sdpsettings('verbose',0);
optimize(A*x<=b,obj,options);

Au = x(1);
Ab = x(3);
hu = x(2);
hb = x(4);

% Account for approximation error corresponding to Proposition 1
hu = hu+sep;
hb = hb-sep;

%% Outputs
% Account for approximation error corresponding to Proposition 1
err = value(x(5))+2*sep;

planes = cell(5,1);

planes{1} = value(Au);
planes{2} = value(hu);
planes{3} = value(Ab);
planes{4} = value(hb);
planes{5} = bound;

end

