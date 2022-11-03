%% Fiacchini

clear
close all
format short g

%% Constants & Data
con = get_constants();

load eq1_6
Reg_1 = {final{1},final{2},final{3},final{5}};

load eq2_6
Reg_2 = {final{1},final{2},final{3},final{5}};

reg_list = cell(1,length(Reg_1));
for i = 1:length(Reg_1)
    reg_list{i} = regions(con,Reg_1{i});
end

%% Target Set/Goal Set
TargetSet = Polyhedron([1 0;-1 0;0 1;0 -1],...
    [con.x1_max;-con.x1_min;con.x2_max;-con.x2_min]);

%% Dynamics

dyn_list = cell(1,length(Reg_1));
for i = 1:length(Reg_1)
    dyn_list{i} = get_dyn(con,Reg_1{i},Reg_2{i});
end

pwd = PwDyn(TargetSet, reg_list, dyn_list);

%% Invariant Set
Xinv = win_always2(pwd, TargetSet, 0.00, 1, 1);
clf;
Xinv(end).plot;