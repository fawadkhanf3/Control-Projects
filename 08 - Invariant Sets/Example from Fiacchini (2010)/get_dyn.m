function dyn = get_dyn(con,hyperplane1,hyperplane2)

dT = con.dT;

au_1 = hyperplane1{1}(1);
bu_1 = hyperplane1{1}(2);
cu_1 = hyperplane1{1}(3);
hu_1 = real(hyperplane1{2});

ab_1 = hyperplane1{3}(1);
bb_1 = hyperplane1{3}(2);
cb_1 = hyperplane1{3}(3);
hb_1 = real(hyperplane1{4});

a_1 = 1/2*(au_1+ab_1);
b_1 = 1/2*(bu_1+bb_1);
c_1 = 1/2*(cu_1+cb_1);
h_1 = 1/2*(hu_1+hb_1);

au_2 = hyperplane2{1}(1);
bu_2 = hyperplane2{1}(2);
cu_2 = hyperplane2{1}(3);
hu_2 = real(hyperplane2{2});

ab_2 = hyperplane2{3}(1);
bb_2 = hyperplane2{3}(2);
cb_2 = hyperplane2{3}(3);
hb_2 = real(hyperplane2{4});

a_2 = 1/2*(au_2+ab_2);
b_2 = 1/2*(bu_2+bb_2);
c_2 = 1/2*(cu_2+cb_2);
h_2 = 1/2*(hu_2+hb_2);

%% State space
A = [a_1 b_1;a_2 b_2];
B = [c_1;c_2];
E = [1 0;0 1];
K = [h_1;h_2];

[~,n]  = size(A);

%% Discretization (Forward Euler Approximation)

Ad = eye(n)+dT*A;
Bd = dT*B;
Ed = dT*E;
Kd = dT*K;

%% Constraints (Control Input)
Huu = [0 0 1;0 0 -1];
huu = [con.u_max; -con.u_min];

XU = Polyhedron(Huu,huu);

%% Disturbances (State-Dependent)
Hdx_minus = [1/2*(ab_1-au_1) 1/2*(bb_1-bu_1);1/2*(ab_2-au_2) 1/2*(bb_2-bu_2)];
Hdx_plus  = [1/2*(-ab_1+au_1) 1/2*(-bb_1+bu_1);1/2*(-ab_2+au_2) 1/2*(-bb_2+bu_2)];

Hdu_minus = [1/2*(cb_1-cu_1);1/2*(cb_2-cu_2)];
Hdu_plus  = [1/2*(-cb_1+cu_1);1/2*(-cb_2+cu_2)];

hd_minus  = [1/2*(hb_1-hu_1);1/2*(hb_2-hu_2)];
hd_plus   = [1/2*(-hb_1+hu_1);1/2*(-hb_2+hu_2)];

XW_V = {[Hdx_plus hd_plus Hdu_plus],[Hdx_minus hd_minus Hdu_minus]};

dyn = Dyn(Ad, Kd, Bd, XU, [], [], [], [], [], [], [], [], Ed, XW_V);

end