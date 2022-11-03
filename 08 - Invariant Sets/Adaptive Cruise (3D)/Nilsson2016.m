%% ACC 2016

clear
close all
format short g

%% Constants & Data
con = get_constants;

%% V-H-VL Space

A = [1 0 0;-1 0 0;0 -1 0;0 0 1;0 0 -1];
b = [con.v_max; -con.v_min; -con.h_min;con.vL_max;-con.vL_min];
VH = Polyhedron(A,b);

% SafeSet
S2 = Polyhedron('H', [con.tau_min -1 0 0]);

%% Dynamics

dyn    = get_dyn(con);
pwadyn = get_pwadyn(con,VH);

%% Specification Sets (Target/Goal Sets)

M11 = Polyhedron('H', [0 -1 0 -con.v_des*con.tau_des]);
M1  = intersect(VH,M11);
% 
G11 = Polyhedron('H', [1 0 0 con.v_des]);
G1  = intersect(VH,G11);
% 
M22 = Polyhedron('H', [0 1 0 con.v_des*con.tau_des]);
M2  = intersect(VH,M22);
% 
G22 = Polyhedron('H', [con.tau_des -1 0 0;1 0 0 con.v_des]);
G2  = intersect(VH,G22);

%% Set-Valued Mappings

C1 = M1;
C2 = M2;

ALL = merge_in(C1, C2);

%% Target Set
Goal1 = intersect(S2,intersect(G1, ALL));
Goal2 = intersect(S2,intersect(G1, ALL));

%% Dynamics

Xinv1 = win_always(dyn, Goal1, 0.00, 1, 1);
% Xinv1_expanded = expand(pwadyn,PolyUnion(Goal1),PolyUnion(Xinv1), 0.00,'plot_stuff');

figure(2)
C1_1 = expand(pwadyn,PolyUnion(intersect(S2, ALL)),PolyUnion(Xinv1), 0.00);

C1_1_M1 = [];
for i = 1:length(C1_1{1})
C1_1_M1 = [C1_1_M1 intersect1(C1_1{1}(i), M1)];
end

C1 = C1_1_M1;

%%

Xinv2 = win_always(dyn, Goal2, 0.00, 1, 1);

% Xinv2_expanded = expand(pwadyn,PolyUnion(Goal2),PolyUnion(Xinv2), 0.00);

figure(4)
C2_1 = expand(pwadyn,PolyUnion(intersect(S2, ALL)),PolyUnion(Xinv2), 0.00);

C2_1_M2 = [];
for i = 1:length(C2_1{1})
C2_1_M2 = [C2_1_M2 intersect1(C1_1{1}(i), M2)];
end

C2 = C2_1_M2;

%% Final Figure

AA = [];
for i = 1:length(C1)
AA = [AA intersect1(C1(i), Polyhedron('H', [0 1 0 200]))];
end

BB = [];
for i = 1:length(AA)
BB = [BB  intersect1(AA(i),Polyhedron('H', [1 0 0 30]))];
end
    
figure(101); clf; hold on
plot(BB, 'color', 'red','alpha',0.3);
plot(C2, 'color', 'black','alpha',0.7)