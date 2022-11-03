%% ACC 2016

clear
% close all
format short g


Ad = [0.5 0;
    1 -0.5];

Bd = [0;0];

Ed = [1;1];

Kd = 0;%[0.5;0.5];

dyn.n = 2;
dyn.A = Ad;
dyn.B = Bd;
dyn.E = Ed;
dyn.K = Kd;

dyn.Hu = [0 0 0;0 0 0];
dyn.hu = [0;0];

dyn.Hdx_minus = [-1 0];
dyn.Hdu_minus = -1;
dyn.Hdx_plus  = [1 0];
dyn.Hdu_plus  = +1;

A = [1 0;-1 0; 0 1;0 -1];
b = [10;10;10;10];
Target = Polyhedron(A,b);

figure(1);
plot(Target, 'alpha', 0.5, 'color', 'blue')
hold on
Inv1 = invariant_set(dyn,Target);
