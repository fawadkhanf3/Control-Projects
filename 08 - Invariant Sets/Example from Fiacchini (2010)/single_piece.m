%% ACC 2014 (Nilsson Dynamics)

clear
close all
format short g

%% Constants & Data
con = get_constants();

load eq1_6
Reg1 = final{1};
load eq2_6
Reg2 = final{1};

%% State Space and Target Set

TargetSet = Polyhedron([1 0; -eye(2);0 1], [con.x1_max; -con.x1_min; -con.x2_min; con.x2_max]);

%% Dynamics

dyn = get_dyn(con,Reg1,Reg2);

%% Invariant Set
Xinv = win_always(dyn, TargetSet, 0.00, 1, 1);
