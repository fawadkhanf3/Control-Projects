%% IGM (SLV 2D - Multi Stage)  - (Without Coasting)

% S-function (dynam)
%=======================================================%

function [sys,x0,str,ts] = InvPend_nonlinearM(t,x,u,flag,xi)

%=======================================================%

switch flag
    
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0
        [sys,x0,str,ts] = mdlInitializeSizes(xi);
        
        %%%%%%%%%%%%%%%
        % Derivatives %
        %%%%%%%%%%%%%%%
    case 1
        sys = mdlDerivatives(t,x,u);
        
        %%%%%%%%%%%%%%
        % Do Nothing %
        %%%%%%%%%%%%%%
    case {2,9}
        sys = [];
        
        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3
        sys = mdlOutputs(t,x,u);
        
        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        error(['Unhandled Flag = ',num2str(flag)]);
end

%end dynam

%=======================================================%

% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.

function [sys,x0,str,ts] = mdlInitializeSizes(xi)

sizes = simsizes; % Calls simsizes for a sizes structure

sizes.NumContStates  = 4; % Number of Continuous States
sizes.NumDiscStates  = 0; % Number of Discrete States
sizes.NumOutputs     = 4; % Number of Outputs
sizes.NumInputs      = 1; % Number of Inputs
sizes.DirFeedthrough = 1;
sizes.NumSampleTimes = 1; % Number of Sample Times

sys = simsizes(sizes); % Converts it to a sizes array.

x0 = xi; % Initializes the initial conditions
str = []; % State ordering strings which is generally specified as []
ts = [0 0]; % Initializes the array of Sample Times

%end mdlInitializeSizes

%=======================================================%

% mdlDerivatives
% Return the derivatives for the continuous states.

function sys = mdlDerivatives(~,x,u)

Mc = 1.0;
Mr = 0.25;
l  = 0.5;
b  = 0.05;
g0 = 9.80665;

y  = x(1);
dy = x(2);

theta  = x(3);
dtheta = x(4);

F = u;

M = [Mc+Mr               0.5*Mr*l*cos(theta);
     0.5*Mr*l*cos(theta) Mr*(l^2)/4];

c = [0.5*Mr*l*dtheta^2*sin(theta) + F - b*dy;
     0.5*Mr*g0*l*sin(theta)];

tmp = M\c;

ddy     = tmp(1);
ddtheta = tmp(2);

sys = [dy;ddy;dtheta;ddtheta];

%end mdlDerivatives

%=======================================================%

% mdlOutputs
% Return the block outputs.

function sys = mdlOutputs(~,x,~)

sys = x;

%end mdlOutputs

%=======================================================%