function d = get_dyn(con)

dT = con.dT;

%% State Space Matrices (Continuous)
A = [-con.f1_bar/con.mass 0;-1 0];
B = [1/con.mass;0];
E = [0;0];
K = [-con.f0_bar/con.mass;con.v_L];

[~,n]  = size(A);

%% Discretization (Forward Euler Approximation)

Ad = eye(n)+dT*A;
Bd = dT*B;
B_scale = max(abs(Bd));
Bd = Bd/B_scale;
Ed = dT*E;
Kd = dT*K;

%% Constraints (Control Input)

Huu = [0 0 1;0 0 -1];
huu = [con.u_max_bar; -con.u_min_bar]*B_scale;

XU = Polyhedron(Huu,huu);

%% Disturbances (State-Dependent)

Hdx_minus = [0 0];
Hdx_plus  = [0 0];

Hdu_minus = 0;
Hdu_plus  = 0;

hd_minus  = 0;
hd_plus   = 0;

XW_V = {[Hdx_plus hd_plus Hdu_plus], [Hdx_minus hd_minus Hdu_minus]};

d = Dyn(Ad, Kd, Bd, XU, [], [], [], [], [], [], [], [], Ed, XW_V);