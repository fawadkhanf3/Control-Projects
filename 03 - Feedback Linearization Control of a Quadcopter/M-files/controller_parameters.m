function fp = controller_parameters(fp)

fp.c0  = 625;
fp.c1  = 500;
fp.c2  = 150;
fp.c3  = 20;
fp.c4  = 4;
fp.c5  = 4;

% syms x y z psi theta phi u v w zeta xi p q r
x     = sym('x');
y     = sym('y');
z     = sym('z');
psi   = sym('psi');
theta = sym('theta');
phi   = sym('phi');
u     = sym('u');
v     = sym('v');
w     = sym('w');
zeta  = sym('zeta');
xi    = sym('xi');
p     = sym('p');
q     = sym('q');
r     = sym('r');

d = [0;0;0];

X = [x;y;z;psi;theta;phi;u;v;w;zeta;xi;p;q;r];
h = [x;y;z;psi];

f = [u;
     v;
     w;
     q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);
     q*cos(phi) - r*sin(phi);
     p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
     d(1)/fp.m - 1/fp.m*(cos(psi)*sin(theta)*cos(phi) + sin(psi)*sin(phi))*zeta;
     d(2)/fp.m - 1/fp.m*(sin(psi)*sin(theta)*cos(phi) - cos(psi)*sin(phi))*zeta;
     d(3)/fp.m + fp.g - 1/fp.m*cos(theta)*cos(phi)*zeta;
     xi;
     0;
     (fp.Iyy-fp.Izz)/fp.Ixx*q*r;
     (fp.Izz-fp.Ixx)/fp.Iyy*p*r;
     (fp.Ixx-fp.Iyy)/fp.Izz*p*q];
 
g1 = [zeros(10,1);1;zeros(3,1)];
g2 = [zeros(11,1);fp.l/fp.Ixx;zeros(2,1)];
g3 = [zeros(12,1);fp.l/fp.Iyy;zeros(1,1)];
g4 = [zeros(13,1);1/fp.Izz];

Lfh  = simplify(jacobian(h,X)*f);
L2fh = simplify(jacobian(Lfh,X)*f);
L3fh = simplify(jacobian(L2fh,X)*f);
L4fh = simplify(jacobian(L3fh,X)*f);

Lg1r3 = simplify(jacobian(L3fh,X)*g1);
Lg2r3 = simplify(jacobian(L3fh,X)*g2);
Lg3r3 = simplify(jacobian(L3fh,X)*g3);
Lg4r3 = simplify(jacobian(L3fh,X)*g4);

deltaXr3 = [Lg1r3 Lg2r3 Lg3r3 Lg4r3]; 

Lg1r1 = simplify(jacobian(Lfh,X)*g1);
Lg2r1 = simplify(jacobian(Lfh,X)*g2);
Lg3r1 = simplify(jacobian(Lfh,X)*g3);
Lg4r1 = simplify(jacobian(Lfh,X)*g4);

deltaXr1 = [Lg1r1 Lg2r1 Lg3r1 Lg4r1];

deltaX = [deltaXr3(1:3,:);deltaXr1(4,:)];

fp.deltaX = matlabFunction(deltaX,'Vars',{X});

b = [L4fh(1:3);L2fh(4)];

fp.b = matlabFunction(b,'Vars',{X});

xyzpsidot    = Lfh;
fp.xyzpsidot = matlabFunction(xyzpsidot,'Vars',{X});

xyz2dot    = L2fh(1:3);
fp.xyz2dot = matlabFunction(xyz2dot,'Vars',{X});

xyz3dot    = L3fh(1:3);
fp.xyz3dot = matlabFunction(xyz3dot,'Vars',{X});

end