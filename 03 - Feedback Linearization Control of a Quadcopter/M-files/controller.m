function [u1,u2,u3,u4] = controller(t,X,ref_traj,fp)

xref   = ref_traj{1};
yref   = ref_traj{2};
zref   = ref_traj{3};
psiref = ref_traj{4};

xyzpsidot = fp.xyzpsidot(X);
xyz2dot   = fp.xyz2dot(X);
xyz3dot   = fp.xyz3dot(X);

v1 = fp.xref4dot(t) - ...
     fp.c3*(xyz3dot(1)   - xref(4)) - ...
     fp.c2*(xyz2dot(1)   - xref(3)) - ...
     fp.c1*(xyzpsidot(1) - xref(2)) - ...
     fp.c0*(X(1)         - xref(1));
 
v2 = fp.yref4dot(t) - ...
     fp.c3*(xyz3dot(2)   - yref(4)) - ...
     fp.c2*(xyz2dot(2)   - yref(3)) - ...
     fp.c1*(xyzpsidot(2) - yref(2)) - ...
     fp.c0*(X(2)         - yref(1));
 
v3 = fp.zref4dot(t) - ...
     fp.c3*(xyz3dot(3)   - zref(4)) - ...
     fp.c2*(xyz2dot(3)   - zref(3)) - ...
     fp.c1*(xyzpsidot(3) - zref(2)) - ...
     fp.c0*(X(3)         - zref(1));
 
v4 = fp.psiref2dot - ...
     fp.c5*(xyzpsidot(4) - psiref(2)) - ...
     fp.c4*(X(4)         - psiref(1));
 
v = [v1;v2;v3;v4];

deltaX = fp.deltaX(X);

b = fp.b(X);

alpha = -inv(deltaX)*b;
beta  = inv(deltaX);

ubar = alpha + beta*v; %#ok<MINV>

u1 = ubar(1);
u2 = ubar(2);
u3 = ubar(3);
u4 = ubar(4);

end