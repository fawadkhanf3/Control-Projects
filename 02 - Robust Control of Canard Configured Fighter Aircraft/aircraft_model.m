function sys = aircraft_model()

% Dimensionless Derivatives

Xu  = 0.050;
Xw  = 0.260;
dXw = 0;
Xq  = 0;
Xde = 0;

Zu  = -1.2;
Zw  = -2.8;
dZw = -0.7;
Zq  = -1.2;
Zde = -0.04;

Mu  = 0.003;
Mw  = 0.280; 
dMw = 0.380;
Mq  = -0.5;
Mde = 0.16;

% Flight Conditions

h       = 0;
gamma   = 0;
alpha_e = 0;
V0      = 123.7;
m       = 12500;
Iy      = 105592;
rho     = 1.225;
S       = 50.0;
c       = 5.7;
g0      = 9.81;

%
theta_e = alpha_e;
Ue = V0*cos(theta_e);
We = V0*sin(theta_e);

XU   = 1/2*rho*V0*S*Xu;
XW   = 1/2*rho*V0*S*Xw;
dXW  = 1/2*rho*S*c*dXw;
XQ   = 1/2*rho*V0*S*c*Xq;
XDE  = 1/2*rho*V0^2*S*Xde;

ZU   = 1/2*rho*V0*S*Zu;
ZW   = 1/2*rho*V0*S*Zw;
dZW  = 1/2*rho*S*c*dZw;
ZQ   = 1/2*rho*V0*S*c*Zq;
ZDE  = 1/2*rho*V0^2*S*Zde;

MU   = 1/2*rho*V0*S*c*Mu;
MW   = 1/2*rho*V0*S*c*Mw;
dMW  = 1/2*rho*S*c^2*dMw;
MQ   = 1/2*rho*V0*S*c^2*Mq;
MDE  = 1/2*rho*V0^2*S*c*Mde;

M = [m -dXW 0 0;
     0 m-dZW 0 0;
     0 -dMW Iy 0;
     0 0 0 1];
 
Abar = [XU XW XQ-m*We -m*g0*cos(theta_e);
        ZU ZW ZQ+m*Ue -m*g0*sin(theta_e);
        MU MW MQ 0;
        0 0 1 0];
    
Bbar = [XDE;ZDE;MDE;0];

A = M\Abar;
B = M\Bbar;

C = eye(4);
C(2,2) = C(2,2)/V0;

D = zeros(4,1);

sys = ss(A,B,C,D);