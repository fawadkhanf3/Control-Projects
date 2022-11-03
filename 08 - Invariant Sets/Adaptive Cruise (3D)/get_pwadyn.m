function pwadyn = get_pwadyn(con,VH)

dT = con.dT;

%% State Space Matrices (Continuous)

A = [-con.f1_bar/con.mass 0 0;-1 0 1;0 0 0];
B = [1/con.mass;0;0];
E = [0;0;1];
K = [-con.f0_bar/con.mass;0;0];

[~,n]  = size(A);

%% Discretization (Forward Euler)

Ad = eye(n)+dT*A;
Bd = dT*B;
B_scale = max(abs(Bd));
Bd = Bd/B_scale;
Ed = dT*E;
Kd = dT*K;

%% Constraints on Control Input

Huu = [0 0 0 1;0 0 0 -1];
huu = [con.u_max_bar; -con.u_min_bar]*B_scale;
XU = Polyhedron(Huu,huu);

%% Piecewise Dynamics

upper_cutoff = con.vL_max - con.aL_max*con.dT;
lower_cutoff = con.vL_min - con.aL_min*con.dT;

reg_upper  = intersect(VH,Polyhedron('A',[0 0 -1],'b',-upper_cutoff));
reg_middle = intersect(VH,Polyhedron('A',[0 0 1;0 0 -1],'b',[upper_cutoff;-lower_cutoff]));
reg_lower  = intersect(VH,Polyhedron('A',[0 0 1],'b',lower_cutoff));

Hdx_minus_upper = [0 0 0];
Hdx_plus_upper  = [0 0 -1/con.dT];

Hdu_minus_upper = 0;
Hdu_plus_upper  = 0;

hd_minus_upper  = con.aL_min;
hd_plus_upper   = con.vL_max/con.dT;

Hdx_minus_middle = [0 0 0];
Hdx_plus_middle  = [0 0 0];

Hdu_minus_middle = 0;
Hdu_plus_middle  = 0;

hd_minus_middle  = con.aL_min;
hd_plus_middle   = con.aL_max;

Hdx_minus_lower = [0 0 -1/con.dT];
Hdx_plus_lower  = [0 0 0];

Hdu_minus_lower = 0;
Hdu_plus_lower  = 0;

hd_minus_lower  = con.vL_min/con.dT;
hd_plus_lower   = con.aL_max;

XW_V_upper = {[Hdx_plus_upper hd_plus_upper Hdu_plus_upper],...
              [Hdx_minus_upper hd_minus_upper Hdu_minus_upper]};

dyn_upper = Dyn(Ad, Kd, Bd, XU, [], [], [], [], [], [], [], [], Ed, XW_V_upper);

XW_V_middle = {[Hdx_plus_middle hd_plus_middle Hdu_plus_middle],...
              [Hdx_minus_middle hd_minus_middle Hdu_minus_middle]};

dyn_middle = Dyn(Ad, Kd, Bd, XU, [], [], [], [], [], [], [], [], Ed, XW_V_middle);

XW_V_lower = {[Hdx_plus_lower hd_plus_lower Hdu_plus_lower],...
              [Hdx_minus_lower hd_minus_lower Hdu_minus_lower]};

dyn_lower = Dyn(Ad, Kd, Bd, XU, [], [], [], [], [], [], [], [], Ed, XW_V_lower);

reg_list = {reg_upper,reg_middle,reg_lower};
dyn_list = {dyn_upper,dyn_middle,dyn_lower};

pwadyn = PwDyn(VH, reg_list, dyn_list);
