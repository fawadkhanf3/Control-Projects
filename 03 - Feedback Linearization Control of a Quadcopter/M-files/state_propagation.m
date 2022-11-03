function X = state_propagation(~,X,control,d,fp)

% Assigning States
x     = X(1);
y     = X(2);
z     = X(3);
psi   = X(4);
theta = X(5);
phi   = X(6);
u     = X(7);
v     = X(8);
w     = X(9);
zeta  = X(10);
xi    = X(11);
p     = X(12);
q     = X(13);
r     = X(14);

% Control Input
u1  = zeta;
u1_ = control{1};
u2  = control{2};
u3  = control{3};
u4  = control{4};

% Differential Equations
xdot     = u;
ydot     = v;
zdot     = w;

psidot   = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);
thetadot = q*cos(phi) - r*sin(phi);
phidot   = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);

udot     = d(1)/fp.m - u1*(1/fp.m)*(cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi));
vdot     = d(2)/fp.m - u1*(1/fp.m)*(sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi));
wdot     = d(3)/fp.m + fp.g - u1*(1/fp.m)*cos(theta)*cos(phi);

zetadot  = xi;
xidot    = u1_;

pdot     = (fp.Iyy-fp.Izz)/fp.Ixx*q*r + u2*fp.l/fp.Ixx;
qdot     = (fp.Izz-fp.Ixx)/fp.Iyy*p*r + u3*fp.l/fp.Iyy;
rdot     = (fp.Ixx-fp.Iyy)/fp.Izz*p*q + u4/fp.Izz;


% Forward Euler Integration
x_next = x+fp.T*xdot;
y_next = y+fp.T*ydot;
z_next = z+fp.T*zdot;

psi_next   = psi+fp.T*psidot;
theta_next = theta+fp.T*thetadot;
phi_next   = phi+fp.T*phidot;

u_next = u+fp.T*udot;
v_next = v+fp.T*vdot;
w_next = w+fp.T*wdot;

zeta_next = zeta+fp.T*zetadot;
xi_next   = xi+fp.T*xidot;

p_next = p+fp.T*pdot;
q_next = q+fp.T*qdot;
r_next = r+fp.T*rdot;

position_next     = [x_next;y_next;z_next];
velocity_next     = [u_next;v_next;w_next];
euler_angles_next = [psi_next;theta_next;phi_next];
rates_next        = [p_next;q_next;r_next];

X = [position_next;euler_angles_next;velocity_next;zeta_next;xi_next;rates_next];

end