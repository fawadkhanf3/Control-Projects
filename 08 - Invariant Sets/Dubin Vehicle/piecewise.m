%% ACC 2014 (Nilsson Dynamics)

clear
close all
format short g

%% Constants & Data
con = get_constants();

pieces = 4;%4;%8;

if pieces == 1
    load hyperplane1
    Reg_1 = {final{1}};
    load hyperplane2
    Reg_2 = {final{1}};
    
elseif pieces == 4
    load hyperplane1_4
    Reg_1 = {final{1},final{2},final{3},final{4}};
    load hyperplane2_4
    Reg_2 = {final{1},final{2},final{3},final{4}};
    
elseif pieces == 8
    load hyperplane1_8
    Reg_1 = {final{1},final{2},final{3},final{4},...
        final{5},final{6},final{7},final{8}};
    load hyperplane2_8
    Reg_2 = {final{1},final{2},final{3},final{4},...
        final{5},final{6},final{7},final{8}};
    
end
    
reg_list = cell(1,length(Reg_1));
for i = 1:length(Reg_1)
    reg_list{i} = regions(con,Reg_1{i});
end

%% Target Set/Goal Set
TargetSet = Polyhedron([1 0 0; -1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1],...
    [con.x_max;-con.x_min;con.y_max;-con.y_min;con.theta_max;-con.theta_min]);

%% Dynamics

dyn_list = cell(1,length(Reg_1));
for i = 1:length(Reg_1)
    dyn_list{i} = get_dyn(con,Reg_1{i},Reg_2{i});
end

pwd = PwDyn(TargetSet, reg_list, dyn_list);

% dd = dyn_list{1};
% rr = reg_list{1};
% Inv1 = win_always(dd,rr,0.00000,1,1);
% 
% return
%% Invariant Set
Xinv = win_always2(pwd, TargetSet, 0.00, 1, 1);
clf;
Xinv(end).plot;