%% ACC 2014 (Nilsson Dynamics)

clear
close all
format short g

%% Constants & Data
con = get_constants();

load hyperplane1
Reg1a = final{1};

load hyperplane2
Reg1b = final{1};

%% Dynamics

dyn = get_dyn(con,Reg1a,Reg1b);

%% State Space and Target Set/Goal Set

Targetset = Polyhedron([1 0 0; -1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1],...
    [con.x_max;-con.x_min;con.y_max;-con.y_min;con.theta_max;-con.theta_min]);

%% Invariant Set

Inv1 = win_always(dyn,Targetset,0.00000,1,1);