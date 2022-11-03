%% ACC 2014 (Nilsson Dynamics)

clear
close all
format short g

%% Constants & Data
con = get_constants();

%% State Space and Target Set/Goal Set

% Define V-H State space (Fig. 2)
VH = Polyhedron([1 0; -eye(2);0 1], [con.v_max; -con.v_min; -con.h_min; con.h_max]);

% Distant Constraint
S1  = Polyhedron('A',[con.tau_min -1;0 -1],'b',[0;-con.h_safe]);

% Safe Set
SafeSet = intersect(VH,S1);

% Goal Set
G2 = Polyhedron('A',[con.tau_des  -1; 1 0],'b',[0; con.v_des]);
TargetSet = intersect(SafeSet,G2);

%% Dynamics

dyn = get_dyn(con);

%% Invariant Set
Xinv = win_always(dyn, TargetSet, 0.00, 1, 1);
